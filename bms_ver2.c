// bms2d.c : Berlekamp-Massey-Sakata(BMS)アルゴリズムによる2次元シンドロームの復号(GF(16) 上)
//
//   gcc -O2 -o bms2d bms2d.c
//   ./bms2d demo 8 4          # D=8、誤り4個: 符号語・受信語・シンドローム・基底・訂正結果を表示
//   ./bms2d check 8 300       # BMS と、ガウス消去による基準版(sakata2d.c と同じ方法)の基底を突き合わせる
//   ./bms2d 8 10 300          # BMS で復号したときの、誤りの個数ごとの成功率の表
//
// 設定: 体 GF(16)(x^4+x+1、alpha = 2)、符号長 15 x 15 の2次元巡回符号。
//   位置 (a,b)、値 e の誤りに対し、シンドローム S[i][j] = sum e * alpha^(a*i + b*j)、既知の領域は i+j < D。
//
// BMS の構成(単項式順序は、全次数、同じ次数では (d,0) < (d-1,1) < ... < (0,d)):
//   F : 現在の極小多項式の集合(Delta の補集合の「角」の単項式ごとに 1 本)
//   G : 補助多項式(過去に不一致を起こした多項式 g と、その添字 q、不一致の値 d_g)。
//       最初は、添字が負の仮想的な失敗 (1, q=(-1,0)) と (1, q=(0,-1)) を入れておく。
//   各点 p を順序どおりに処理する:
//     1. F の各 f(先頭項 s <= p のもの)について、不一致 d_f = sum f_k S[k + p - s] を計算する。
//     2. すべて 0 なら何もしない。
//     3. そうでなければ、単項式 r のうち「修正できないもの」を新しい Delta に加える。
//        r が修正できるのは、(a) s(f) <= r で d_f = 0 の f がある、(b) r が p 以下でない、
//        (c) G に q - s(g) >= p - r を満たす g がある、のいずれか。
//     4. 新しい Delta の補集合の角 s' ごとに、新しい多項式を作る:
//          f' = x^a f - (d_f / d_g) x^b g   (a = s' - s(f), b = s' - (p - q) - s(g) >= 0)
//        (d_f = 0 か s' が p 以下でないときは f' = x^a f)
//     5. 不一致を起こした古い f を G に加える。
//   最終的な F の共通零点(15x15 の全点を試す)が誤り位置、続いて連立方程式で誤りの値を求める。
//
// 注意: 動作確認用。乱数は rand()。基準版との突き合わせで、同じ基底が出ることを確かめている。

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>

#define NN 15
#define MAXE 24 // 指数の上限
#define MAXD 40 // Delta の大きさの上限(誤りの個数)
#define MAXF 64
#define MAXG 6000
#define MAXS 24 // S の添字の上限

// ---------------------------------------------------------------- GF(16)
static int EXPT[2 * NN], LOGT[16];
static void gf_init(void)
{
    int x = 1;
    for (int i = 0; i < NN; i++)
    {
        EXPT[i] = x;
        LOGT[x] = i;
        x <<= 1;
        if (x & 16)
            x ^= 0x13;
    }
    for (int i = NN; i < 2 * NN; i++)
        EXPT[i] = EXPT[i - NN];
}
static int gmul(int a, int b) { return (a && b) ? EXPT[LOGT[a] + LOGT[b]] : 0; }
static int ginv(int a) { return EXPT[(NN - LOGT[a]) % NN]; }
static int apow(int e)
{
    e %= NN;
    if (e < 0)
        e += NN;
    return EXPT[e];
}

// 連立方程式(加法は XOR)。A は R x (C+1) の拡大係数行列。整合していれば 1
static int solve(int A[][MAXD + 1], int R, int C, int *x)
{
    int piv[MAXD + 1], rank = 0;
    for (int c = 0; c < C && rank < R; c++)
    {
        int p = -1;
        for (int i = rank; i < R; i++)
            if (A[i][c])
            {
                p = i;
                break;
            }
        if (p < 0)
            continue;
        if (p != rank)
            for (int j = 0; j <= C; j++)
            {
                int t = A[p][j];
                A[p][j] = A[rank][j];
                A[rank][j] = t;
            }
        int iv = ginv(A[rank][c]);
        for (int j = 0; j <= C; j++)
            A[rank][j] = gmul(A[rank][j], iv);
        for (int i = 0; i < R; i++)
            if (i != rank && A[i][c])
            {
                int f = A[i][c];
                for (int j = 0; j <= C; j++)
                    A[i][j] ^= gmul(f, A[rank][j]);
            }
        piv[rank++] = c;
    }
    for (int i = rank; i < R; i++)
        if (A[i][C])
            return 0;
    for (int j = 0; j < C; j++)
        x[j] = 0;
    for (int i = 0; i < rank; i++)
        x[piv[i]] = A[i][C];
    return 1;
}

// ---------------------------------------------------------------- 多項式(係数の2次元配列)
typedef struct
{
    uint8_t c[MAXE][MAXE];
    int s1, s2;
} poly; // s = 先頭項の指数
static int D;
static int S[MAXS][MAXS];

static void poly_one(poly *f)
{
    memset(f, 0, sizeof *f);
    f->c[0][0] = 1;
    f->s1 = f->s2 = 0;
}
static int leq(int a1, int a2, int b1, int b2) { return a1 <= b1 && a2 <= b2; }

static int disc(const poly *f, int p1, int p2) // 前提: s <= p
{
    int acc = 0;
    for (int i = 0; i < MAXE; i++)
        for (int j = 0; j < MAXE; j++)
            if (f->c[i][j])
                acc ^= gmul(f->c[i][j], S[i + p1 - f->s1][j + p2 - f->s2]);
    return acc;
}
// dst ^= coef * x^(e1,e2) * src
static void add_shifted(poly *dst, const poly *src, int e1, int e2, int coef)
{
    for (int i = 0; i < MAXE; i++)
        for (int j = 0; j < MAXE; j++)
            if (src->c[i][j])
            {
                int a = i + e1, b = j + e2;
                if (a < 0 || b < 0 || a >= MAXE || b >= MAXE)
                {
                    fprintf(stderr, "指数が範囲外\n");
                    exit(1);
                }
                dst->c[a][b] ^= gmul(src->c[i][j], coef);
            }
}

// ---------------------------------------------------------------- BMS
static poly Fs[MAXF], NF[MAXF];
static int nF;
static int delta[MAXE][MAXE], dn;
static poly Gp[MAXG];
static int Gq1[MAXG], Gq2[MAXG], Gd[MAXG], nG;

static int bms(void)
{
    nF = 1;
    poly_one(&Fs[0]);
    nG = 2;
    poly_one(&Gp[0]);
    Gq1[0] = -1;
    Gq2[0] = 0;
    Gd[0] = 1;
    poly_one(&Gp[1]);
    Gq1[1] = 0;
    Gq2[1] = -1;
    Gd[1] = 1;
    memset(delta, 0, sizeof delta);
    for (int dg = 0; dg < D; dg++)
        for (int pi = dg; pi >= 0; pi--)
        {
            int p1 = pi, p2 = dg - pi;
            int ds[MAXF], any = 0;
            for (int k = 0; k < nF; k++)
            {
                ds[k] = leq(Fs[k].s1, Fs[k].s2, p1, p2) ? disc(&Fs[k], p1, p2) : 0;
                if (ds[k])
                    any = 1;
            }
            if (!any)
                continue;
            // 1. 新しい Delta
            static int nd[MAXE][MAXE];
            memcpy(nd, delta, sizeof nd);
            for (int d2 = 0; d2 < D; d2++)
                for (int i2 = d2; i2 >= 0; i2--)
                {
                    int r1 = i2, r2 = d2 - i2;
                    if (delta[r1][r2])
                        continue;
                    int zero = 0;
                    for (int k = 0; k < nF; k++)
                        if (leq(Fs[k].s1, Fs[k].s2, r1, r2) && ds[k] == 0)
                        {
                            zero = 1;
                            break;
                        }
                    if (zero)
                        continue;
                    if (!leq(r1, r2, p1, p2))
                        continue;
                    int ok = 0;
                    for (int g = 0; g < nG; g++)
                    {
                        poly *gg = &Gp[g];
                        if (Gq1[g] - gg->s1 >= p1 - r1 && Gq2[g] - gg->s2 >= p2 - r2)
                        {
                            ok = 1;
                            break;
                        }
                    }
                    if (ok)
                        continue;
                    nd[r1][r2] = 1;
                }
            // 2. 角ごとに新しい多項式
            int nn = 0;
            for (int d2 = 0; d2 < D; d2++)
                for (int i2 = d2; i2 >= 0; i2--)
                {
                    int s1 = i2, s2 = d2 - i2;
                    if (nd[s1][s2])
                        continue;
                    if (!((s1 == 0 || nd[s1 - 1][s2]) && (s2 == 0 || nd[s1][s2 - 1])))
                        continue; // 角ではない
                    if (nn >= MAXF)
                        return 0;
                    poly *out = &NF[nn];
                    int fi = -1;
                    for (int k = 0; k < nF; k++)
                        if (Fs[k].s1 == s1 && Fs[k].s2 == s2)
                        {
                            fi = k;
                            break;
                        }
                    if (fi >= 0)
                    { // 古い角
                        if (ds[fi] == 0)
                        {
                            *out = Fs[fi];
                            nn++;
                            continue;
                        }
                    }
                    else
                    { // 新しい角: s(f) <= s' の f を探す
                        int zi = -1, ci = -1;
                        for (int k = 0; k < nF; k++)
                            if (leq(Fs[k].s1, Fs[k].s2, s1, s2))
                            {
                                if (ci < 0)
                                    ci = k;
                                if (ds[k] == 0 && zi < 0)
                                    zi = k;
                            }
                        if (ci < 0)
                            return 0;
                        if (zi >= 0 || !leq(s1, s2, p1, p2))
                        { // 修正不要
                            int k = zi >= 0 ? zi : ci;
                            memset(out, 0, sizeof *out);
                            add_shifted(out, &Fs[k], s1 - Fs[k].s1, s2 - Fs[k].s2, 1);
                            out->s1 = s1;
                            out->s2 = s2;
                            nn++;
                            continue;
                        }
                        fi = ci;
                    }
                    // g を使って打ち消す
                    int gi = -1, b1 = 0, b2 = 0;
                    for (int g = 0; g < nG; g++)
                    {
                        b1 = s1 - (p1 - Gq1[g]) - Gp[g].s1;
                        b2 = s2 - (p2 - Gq2[g]) - Gp[g].s2;
                        if (b1 >= 0 && b2 >= 0)
                        {
                            gi = g;
                            break;
                        }
                    }
                    if (gi < 0)
                        return 0;
                    int coef = gmul(ds[fi], ginv(Gd[gi]));
                    memset(out, 0, sizeof *out);
                    add_shifted(out, &Fs[fi], s1 - Fs[fi].s1, s2 - Fs[fi].s2, 1);
                    add_shifted(out, &Gp[gi], b1, b2, coef);
                    out->s1 = s1;
                    out->s2 = s2;
                    nn++;
                }
            // 3. 不一致を起こした古い f を G に加える
            for (int k = 0; k < nF; k++)
                if (ds[k])
                {
                    if (nG >= MAXG)
                        return 0;
                    Gp[nG] = Fs[k];
                    Gq1[nG] = p1;
                    Gq2[nG] = p2;
                    Gd[nG] = ds[k];
                    nG++;
                }
            for (int k = 0; k < nn; k++)
                Fs[k] = NF[k];
            nF = nn;
            memcpy(delta, nd, sizeof delta);
        }
    dn = 0;
    for (int i = 0; i < MAXE; i++)
        for (int j = 0; j < MAXE; j++)
            dn += delta[i][j];
    return 1;
}

// ---------------------------------------------------------------- 基準: ガウス消去版(sakata2d.c と同じ方法)
static poly OB[MAXF];
static int nOB;
static int odelta[MAXE][MAXE], odn;
static int oracle(void)
{
    static int A[MAXS * MAXS][MAXD + 1];
    int dl1[MAXD], dl2[MAXD];
    nOB = 0;
    odn = 0;
    memset(odelta, 0, sizeof odelta);
    for (int d = 0; d <= D - 1; d++)
    {
        int ms[MAXS * MAXS][2], nm = 0;
        for (int m1 = 0; m1 <= D - 1 - d; m1++)
            for (int m2 = 0; m1 + m2 <= D - 1 - d; m2++)
            {
                ms[nm][0] = m1;
                ms[nm][1] = m2;
                nm++;
            }
        for (int i = d; i >= 0; i--)
        {
            int r1 = i, r2 = d - i, div = 0;
            for (int b = 0; b < nOB; b++)
                if (leq(OB[b].s1, OB[b].s2, r1, r2))
                    div = 1;
            if (div)
                continue;
            if (odn >= MAXD)
                return 0;
            for (int row = 0; row < nm; row++)
            {
                for (int t = 0; t < odn; t++)
                    A[row][t] = S[dl1[t] + ms[row][0]][dl2[t] + ms[row][1]];
                A[row][odn] = S[r1 + ms[row][0]][r2 + ms[row][1]];
            }
            int x[MAXD + 1];
            if (solve(A, nm, odn, x))
            {
                poly *f = &OB[nOB++];
                memset(f, 0, sizeof *f);
                f->c[r1][r2] = 1;
                f->s1 = r1;
                f->s2 = r2;
                for (int t = 0; t < odn; t++)
                    f->c[dl1[t]][dl2[t]] = (uint8_t)x[t];
            }
            else
            {
                dl1[odn] = r1;
                dl2[odn] = r2;
                odelta[r1][r2] = 1;
                odn++;
            }
        }
    }
    return 1;
}

// 基底を正規形(先頭係数 1、先頭項以外は他の多項式の先頭項で割れない形)にする
static void normal_form(poly *B, int n)
{
    for (int i = 0; i < n; i++)
    {
        int iv = ginv(B[i].c[B[i].s1][B[i].s2]);
        for (int a = 0; a < MAXE; a++)
            for (int b = 0; b < MAXE; b++)
                B[i].c[a][b] = (uint8_t)gmul(B[i].c[a][b], iv);
    }
    for (int i = 0; i < n; i++)
    {
        int changed = 1;
        while (changed)
        {
            changed = 0;
            for (int d = 2 * MAXE; d >= 0 && !changed; d--)
                for (int a = 0; a <= d && !changed; a++)
                {
                    int k1 = a, k2 = d - a; // 次数 d の単項式を大きい順(i が小さいほど大きい)に調べる
                    if (k1 >= MAXE || k2 >= MAXE || !B[i].c[k1][k2])
                        continue;
                    if (k1 == B[i].s1 && k2 == B[i].s2)
                        continue;
                    for (int j = 0; j < n; j++)
                        if (j != i && leq(B[j].s1, B[j].s2, k1, k2))
                        {
                            add_shifted(&B[i], &B[j], k1 - B[j].s1, k2 - B[j].s2, B[i].c[k1][k2]);
                            changed = 1;
                            break;
                        }
                }
        }
    }
}
static int key_lt(const poly *a, const poly *b)
{
    int da = a->s1 + a->s2, db = b->s1 + b->s2;
    return da != db ? da < db : a->s1 > b->s1;
}
static int same_basis(void)
{
    if (nF != nOB || dn != odn)
        return 0;
    static poly A[MAXF], B[MAXF];
    memcpy(A, Fs, sizeof(poly) * nF);
    memcpy(B, OB, sizeof(poly) * nOB);
    normal_form(A, nF);
    normal_form(B, nOB);
    for (int i = 0; i < nF; i++)
        for (int j = i + 1; j < nF; j++)
        {
            if (key_lt(&A[j], &A[i]))
            {
                poly t = A[i];
                A[i] = A[j];
                A[j] = t;
            }
            if (key_lt(&B[j], &B[i]))
            {
                poly t = B[i];
                B[i] = B[j];
                B[j] = t;
            }
        }
    for (int i = 0; i < nF; i++)
        if (memcmp(A[i].c, B[i].c, sizeof A[i].c))
            return 0;
    for (int i = 0; i < MAXE; i++)
        for (int j = 0; j < MAXE; j++)
            if (delta[i][j] != odelta[i][j])
                return 0;
    return 1;
}

// ---------------------------------------------------------------- 復号
typedef struct
{
    int a, b, v;
} err;
static int decode(err *out)
{
    int allzero = 1;
    for (int i = 0; i < D; i++)
        for (int j = 0; i + j < D; j++)
            if (S[i][j])
                allzero = 0;
    if (allzero)
        return 0;
    if (!bms() || dn == 0)
        return -1;
    int pa[NN * NN], pb[NN * NN], zn = 0;
    for (int a = 0; a < NN; a++)
        for (int b = 0; b < NN; b++)
        {
            int ok = 1;
            for (int k = 0; k < nF && ok; k++)
            {
                int v = 0;
                for (int i = 0; i < MAXE; i++)
                    for (int j = 0; j < MAXE; j++)
                        if (Fs[k].c[i][j])
                            v ^= gmul(Fs[k].c[i][j], apow(a * i + b * j));
                if (v)
                    ok = 0;
            }
            if (ok)
            {
                if (zn >= MAXD)
                    return -1;
                pa[zn] = a;
                pb[zn] = b;
                zn++;
            }
        }
    if (zn != dn)
        return -1;
    static int A[MAXS * MAXS][MAXD + 1];
    int R = 0;
    for (int i = 0; i < D; i++)
        for (int j = 0; i + j < D; j++)
        {
            for (int k = 0; k < zn; k++)
                A[R][k] = apow(pa[k] * i + pb[k] * j);
            A[R][zn] = S[i][j];
            R++;
        }
    int x[MAXD + 1];
    if (!solve(A, R, zn, x))
        return -1;
    for (int k = 0; k < zn; k++)
    {
        if (!x[k])
            return -1;
        out[k].a = pa[k];
        out[k].b = pb[k];
        out[k].v = x[k];
    }
    return zn;
}
static void make_syndromes(const err *e, int w)
{
    memset(S, 0, sizeof S);
    for (int i = 0; i < D; i++)
        for (int j = 0; i + j < D; j++)
        {
            int s = 0;
            for (int k = 0; k < w; k++)
                s ^= gmul(e[k].v, apow(e[k].a * i + e[k].b * j));
            S[i][j] = s;
        }
}
static void random_errors(err *e, int w)
{
    int used[NN][NN];
    memset(used, 0, sizeof used);
    for (int k = 0; k < w; k++)
    {
        int a, b;
        do
        {
            a = rand() % NN;
            b = rand() % NN;
        } while (used[a][b]);
        used[a][b] = 1;
        e[k].a = a;
        e[k].b = b;
        e[k].v = 1 + rand() % 15;
    }
}
static int same_set(const err *x, int nx, const err *y, int ny)
{
    if (nx != ny)
        return 0;
    for (int i = 0; i < nx; i++)
    {
        int f = 0;
        for (int j = 0; j < ny; j++)
            if (x[i].a == y[j].a && x[i].b == y[j].b && x[i].v == y[j].v)
                f = 1;
        if (!f)
            return 0;
    }
    return 1;
}
static void print_poly(const poly *f)
{
    printf("x^%d y^%d", f->s1, f->s2);
    for (int d = 2 * MAXE; d >= 0; d--)
        for (int a = 0; a <= d; a++)
        {
            int i = a, j = d - a;
            if (i >= MAXE || j >= MAXE || !f->c[i][j] || (i == f->s1 && j == f->s2))
                continue;
            printf(" + %d*x^%d y^%d", f->c[i][j], i, j);
        }
}

// ---------------------------------------------------------------- 符号語の生成と表示(デモ用)
// 符号 C = { c (15x15、GF(16) 上) : c(alpha^i, alpha^j) = 0  (i+j < D の全ての (i,j)) }
// 次元は 225 - D(D+1)/2。ランダムな語 u を、H x = H u の解 x で補正して、符号語 c = u + x を作る。
#define MAXT 260
static void gen_codeword(int c[NN][NN])
{
    static int M[MAXT][NN * NN + 1];
    int u[NN][NN], T = 0;
    for (int a = 0; a < NN; a++)
        for (int b = 0; b < NN; b++)
            u[a][b] = rand() % 16;
    for (int i = 0; i < D; i++)
        for (int j = 0; i + j < D; j++)
        {
            int s = 0;
            for (int a = 0; a < NN; a++)
                for (int b = 0; b < NN; b++)
                {
                    int h = apow(a * i + b * j);
                    M[T][a * NN + b] = h;
                    s ^= gmul(u[a][b], h);
                }
            M[T][NN * NN] = s;
            T++;
        }
    int piv[MAXT], rank = 0;
    for (int col = 0; col < NN * NN && rank < T; col++)
    {
        int p = -1;
        for (int r = rank; r < T; r++)
            if (M[r][col])
            {
                p = r;
                break;
            }
        if (p < 0)
            continue;
        if (p != rank)
            for (int j = 0; j <= NN * NN; j++)
            {
                int t = M[p][j];
                M[p][j] = M[rank][j];
                M[rank][j] = t;
            }
        int iv = ginv(M[rank][col]);
        for (int j = 0; j <= NN * NN; j++)
            M[rank][j] = gmul(M[rank][j], iv);
        for (int r = 0; r < T; r++)
            if (r != rank && M[r][col])
            {
                int f = M[r][col];
                for (int j = 0; j <= NN * NN; j++)
                    M[r][j] ^= gmul(f, M[rank][j]);
            }
        piv[rank++] = col;
    }
    int x[NN * NN];
    memset(x, 0, sizeof x);
    for (int r = 0; r < rank; r++)
        x[piv[r]] = M[r][NN * NN];
    for (int a = 0; a < NN; a++)
        for (int b = 0; b < NN; b++)
            c[a][b] = u[a][b] ^ x[a * NN + b];
}
static void syndromes_of(const int w[NN][NN])
{
    memset(S, 0, sizeof S);
    for (int i = 0; i < D; i++)
        for (int j = 0; i + j < D; j++)
        {
            int s = 0;
            for (int a = 0; a < NN; a++)
                for (int b = 0; b < NN; b++)
                    if (w[a][b])
                        s ^= gmul(w[a][b], apow(a * i + b * j));
            S[i][j] = s;
        }
}
static void print_grid(const char *title, const int g[NN][NN], const err *mark, int nm)
{
    printf("%s(行 a、列 b、値は 16 進。[ ] は誤りの位置)\n     ", title);
    for (int b = 0; b < NN; b++)
        printf(" %X ", b);
    printf("\n");
    for (int a = 0; a < NN; a++)
    {
        printf("  %X  ", a);
        for (int b = 0; b < NN; b++)
        {
            int m = 0;
            for (int k = 0; k < nm; k++)
                if (mark[k].a == a && mark[k].b == b)
                    m = 1;
            if (m)
                printf("[%X]", g[a][b]);
            else
                printf(" %X ", g[a][b]);
        }
        printf("\n");
    }
}
static void print_syndrome_table(void)
{
    printf("シンドローム S[i][j](i+j < %d、値は 16 進):\n     ", D);
    for (int j = 0; j < D; j++)
        printf(" j=%X", j);
    printf("\n");
    for (int i = 0; i < D; i++)
    {
        printf("  i=%X", i);
        for (int j = 0; j < D; j++)
            if (i + j < D)
                printf("  %X  ", S[i][j]);
            else
                printf("     ");
        printf("\n");
    }
}

int main(int argc, char **argv)
{
    gf_init();
    srand((unsigned)time(NULL));
    int argi = 1;
    const char *mode = "table";
    if (argc > 1 && (!strcmp(argv[1], "demo") || !strcmp(argv[1], "check")))
    {
        mode = argv[1];
        argi = 2;
    }
    D = argc > argi ? atoi(argv[argi]) : 8;
    if (D < 2 || D > MAXS - 2)
    {
        fprintf(stderr, "D は 2..%d\n", MAXS - 2);
        return 2;
    }

 if (!strcmp(mode, "demo"))
{
    int w = argc > argi + 1 ? atoi(argv[argi + 1]) : 4;
    if (w < 1 || w > MAXD) {
        fprintf(stderr, "誤り数は 1..%d を指定してください\n", MAXD);
        return 2;
    }

    err e[MAXD], f[MAXD];
    int c[NN][NN], r[NN][NN], fixed[NN][NN];

    gen_codeword(c);
    random_errors(e, w);

    memcpy(r, c, sizeof c);
    for (int k = 0; k < w; k++)
        r[e[k].a][e[k].b] ^= e[k].v;

    // 符号語 → シンドローム → エラーの順で表示
    syndromes_of(r);

    printf("2次元巡回符号(GF(16) 上、15x15、零点 i+j<%d の %d 個): "
           "次元 %d、誤り %d 個\n\n",
           D, D * (D + 1) / 2,
           NN * NN - D * (D + 1) / 2, w);

    print_grid("送信した符号語", c, NULL, 0);
    printf("\n");

    print_syndrome_table();
    printf("\n");

    printf("入れたエラー:\n");
    for (int k = 0; k < w; k++)
        printf("  位置 (%2d,%2d)  値 %2d\n", e[k].a, e[k].b, e[k].v);

    printf("\n");
    print_grid("受信語", r, e, w);
    printf("\n");
    fflush(stdout);

    int n = decode(f);  // 復号は一度だけ

    printf("BMS が求めた極小多項式(%d 本、補助多項式を %d 個ためた):\n",
           nF, nG);
    for (int k = 0; k < nF; k++) {
        printf("  ");
        print_poly(&Fs[k]);
        printf("\n");
    }

    printf("Delta(誤りの個数 %d):", dn);
    for (int d = 0; d < D; d++)
        for (int i = d; i >= 0; i--)
            if (delta[i][d - i])
                printf(" x^%d y^%d", i, d - i);
    printf("\n\n復号結果:");

    if (n < 0) {
        printf(" 失敗(訂正できない)\n");
        return 0;
    }

    printf(" %d 個の誤りを訂正\n", n);
    for (int k = 0; k < n; k++)
        printf("  位置 (%2d,%2d)  値 %2d\n",
               f[k].a, f[k].b, f[k].v);

    memcpy(fixed, r, sizeof r);
    for (int k = 0; k < n; k++)
        fixed[f[k].a][f[k].b] ^= f[k].v;

    printf("\n");
    print_grid("訂正後の語", fixed, f, n);
    printf("\n見つけた誤りは入れた誤りと一致: %s\n",
           same_set(e, w, f, n) ? "はい" : "いいえ");
    printf("訂正後の語は送信した符号語と一致: %s\n",
           memcmp(fixed, c, sizeof c) ? "いいえ" : "はい");

    oracle();
    printf("基準版(ガウス消去)と同じ基底: %s\n",
           same_basis() ? "はい" : "いいえ");
    return 0;
}

    if (!strcmp(mode, "check"))
    {
        int trials = argc > argi + 1 ? atoi(argv[argi + 1]) : 300;
        int bad = 0, fail = 0;
        double tb = 0, to = 0;
        for (int t = 0; t < trials; t++)
        {
            int w = 1 + rand() % 10;
            err e[MAXD];
            random_errors(e, w);
            make_syndromes(e, w);
            clock_t c0 = clock();
            int ok1 = bms();
            clock_t c1 = clock();
            int ok2 = oracle();
            clock_t c2 = clock();
            tb += (double)(c1 - c0) / CLOCKS_PER_SEC;
            to += (double)(c2 - c1) / CLOCKS_PER_SEC;
            if (!ok1 || !ok2)
            {
                fail++;
                continue;
            }
            if (!same_basis())
            {
                bad++;
                if (bad <= 3)
                    printf("不一致: 試行 %d、誤り %d 個\n", t, w);
            }
        }
        printf("D=%d、%d 回(誤り 1..10 個):基準版と基底が異なる %d 回、どちらかが範囲超過 %d 回\n", D, trials, bad, fail);
        printf("計算時間の合計: BMS %.3f 秒、ガウス消去版 %.3f 秒\n", tb, to);
        return bad != 0;
    }
    int wmax = argc > argi + 1 ? atoi(argv[argi + 1]) : 10;
    int trials = argc > argi + 2 ? atoi(argv[argi + 2]) : 500;
    printf("BMS による復号: GF(16) 上、15x15、既知のシンドローム i+j < %d の %d 個\n\n", D, D * (D + 1) / 2);
    printf("誤りの個数 | 訂正成功 | 検出失敗 | 誤訂正\n");
    for (int w = 1; w <= wmax; w++)
    {
        int succ = 0, fail = 0, mis = 0;
        for (int t = 0; t < trials; t++)
        {
            err e[MAXD], f[MAXD];
            random_errors(e, w);
            make_syndromes(e, w);
            int n = decode(f);
            if (n < 0)
                fail++;
            else if (same_set(e, w, f, n))
                succ++;
            else
                mis++;
        }
        printf("    %3d     | %6.1f%%  | %6.1f%%  | %6.1f%%\n", w, 100.0 * succ / trials, 100.0 * fail / trials, 100.0 * mis / trials);
    }
    return 0;
}