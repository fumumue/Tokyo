#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>

#define Q 16
#define NN 15
#define MAXD 40
#define MAXB 64
#define MAXS 40

// ---------------------------------------------------------------- GF(16)
static int EXPT[2 * NN], LOGT[Q];
static void gf_init(void)
{
    int x = 1;
    for (int i = 0; i < NN; i++) { EXPT[i] = x; LOGT[x] = i; x <<= 1; if (x & 16) x ^= 0x13; }
    for (int i = NN; i < 2 * NN; i++) EXPT[i] = EXPT[i - NN];
}
static int gmul(int a, int b) { return (a && b) ? EXPT[LOGT[a] + LOGT[b]] : 0; }
static int ginv(int a) { return EXPT[(NN - LOGT[a]) % NN]; }
static int apow(int e) { e %= NN; if (e < 0) e += NN; return EXPT[e]; }

// ---------------------------------------------------------------- 連立方程式 (GF(16) 上、加法は XOR)
static int solve(int A[][MAXD + 1], int R, int C, int *x)
{
    int piv[MAXD + 1], rank = 0;
    for (int c = 0; c < C && rank < R; c++) {
        int p = -1;
        for (int i = rank; i < R; i++) if (A[i][c]) { p = i; break; }
        if (p < 0) continue;
        if (p != rank) for (int j = 0; j <= C; j++) { int t = A[p][j]; A[p][j] = A[rank][j]; A[rank][j] = t; }
        int iv = ginv(A[rank][c]);
        for (int j = 0; j <= C; j++) A[rank][j] = gmul(A[rank][j], iv);
        for (int i = 0; i < R; i++)
            if (i != rank && A[i][c]) { int f = A[i][c]; for (int j = 0; j <= C; j++) A[i][j] ^= gmul(f, A[rank][j]); }
        piv[rank++] = c;
    }
    for (int i = rank; i < R; i++) if (A[i][C]) return 0;
    for (int j = 0; j < C; j++) x[j] = 0;
    for (int i = 0; i < rank; i++) x[piv[i]] = A[i][C];
    return 1;
}

// ----------------------------------------------------------------
// 2次元 Sakata/BMS 風の誤り位置イデアル作成
//
// これは、元の sakata2d.c を「BMS 風」に整理した実験版です。
// つまり、理論的には BM-Sakata に近い差分更新の視点を使いつつ、
// 既存の線形漸化式判定を簡潔に保つために、依存判定を Gauss 消去で行います。
//
// 完全な最適化版の BM-Sakata では、線形系を毎回全体解くのではなく、
// discrepancy の更新を逐次適用して O(n^2) 近似で進めます。
// ここでは、実験コードとして「BMS の考え方」を反映しつつ、
// 既存のデコード・検証フレームワークを壊さないようにしています。
// ----------------------------------------------------------------
typedef struct { int r1, r2; int c[MAXD]; } poly;
static int S[MAXS][MAXS];
static int D;
static int dn, dl1[MAXD], dl2[MAXD];
static poly basis[MAXB];
static int nb;

static int divisible(int r1, int r2)
{
    for (int b = 0; b < nb; b++) if (r1 >= basis[b].r1 && r2 >= basis[b].r2) return 1;
    return 0;
}

// 2D Sakata 風の基底構築。
// ここでの列挙順は、単項式を total degree で昇順に並べ、同一 degree では
// x^r y^(d-r) の順に処理する。
static int build_basis_bms(void)
{
    static int A[MAXS * MAXS][MAXD + 1];
    nb = 0; dn = 0;

    for (int d = 0; d <= D - 1; d++) {
        int ms[MAXS * MAXS][2], nm = 0;
        for (int m1 = 0; m1 <= D - 1 - d; m1++)
            for (int m2 = 0; m1 + m2 <= D - 1 - d; m2++) {
                ms[nm][0] = m1; ms[nm][1] = m2; nm++;
            }

        for (int i = d; i >= 0; i--) {
            int r1 = i, r2 = d - i;
            if (divisible(r1, r2)) continue;
            if (dn >= MAXD) return 0;

            for (int row = 0; row < nm; row++) {
                for (int t = 0; t < dn; t++)
                    A[row][t] = S[dl1[t] + ms[row][0]][dl2[t] + ms[row][1]];
                A[row][dn] = S[r1 + ms[row][0]][r2 + ms[row][1]];
            }

            int x[MAXD + 1];
            if (solve(A, nm, dn, x)) {
                poly *f = &basis[nb++];
                f->r1 = r1; f->r2 = r2;
                for (int t = 0; t < MAXD; t++) f->c[t] = t < dn ? x[t] : 0;
            } else {
                dl1[dn] = r1; dl2[dn] = r2; dn++;
            }
            if (nb >= MAXB) return 0;
        }
    }
    return 1;
}

// ---------------------------------------------------------------- 復号
typedef struct { int a, b, v; } err;
static int decode(err *out)
{
    int allzero = 1;
    for (int i = 0; i < D; i++) for (int j = 0; i + j < D; j++) if (S[i][j]) allzero = 0;
    if (allzero) return 0;

    if (!build_basis_bms()) return -1;
    if (dn == 0) return -1;

    int pa[NN * NN], pb[NN * NN], zn = 0;
    for (int a = 0; a < NN; a++) for (int b = 0; b < NN; b++) {
        int ok = 1;
        for (int k = 0; k < nb && ok; k++) {
            int v = apow(a * basis[k].r1 + b * basis[k].r2);
            for (int t = 0; t < dn; t++) if (basis[k].c[t])
                v ^= gmul(basis[k].c[t], apow(a * dl1[t] + b * dl2[t]));
            if (v) ok = 0;
        }
        if (ok) { if (zn >= MAXD) return -1; pa[zn] = a; pb[zn] = b; zn++; }
    }
    if (zn != dn) return -1;

    static int A[MAXS * MAXS][MAXD + 1];
    int R = 0;
    for (int i = 0; i < D; i++) for (int j = 0; i + j < D; j++) {
        for (int k = 0; k < zn; k++) A[R][k] = apow(pa[k] * i + pb[k] * j);
        A[R][zn] = S[i][j]; R++;
    }
    int x[MAXD + 1];
    if (!solve(A, R, zn, x)) return -1;

    for (int k = 0; k < zn; k++) {
        if (!x[k]) return -1;
        out[k].a = pa[k]; out[k].b = pb[k]; out[k].v = x[k];
    }
    return zn;
}

static void make_syndromes(const err *e, int w)
{
    memset(S, 0, sizeof S);
    for (int i = 0; i < D; i++) for (int j = 0; i + j < D; j++) {
        int s = 0;
        for (int k = 0; k < w; k++) s ^= gmul(e[k].v, apow(e[k].a * i + e[k].b * j));
        S[i][j] = s;
    }
}

static void random_errors(err *e, int w)
{
    int used[NN][NN]; memset(used, 0, sizeof used);
    for (int k = 0; k < w; k++) {
        int a, b; do { a = rand() % NN; b = rand() % NN; } while (used[a][b]);
        used[a][b] = 1; e[k].a = a; e[k].b = b; e[k].v = 1 + rand() % 15;
    }
}

static int same_set(const err *x, int nx, const err *y, int ny)
{
    if (nx != ny) return 0;
    for (int i = 0; i < nx; i++) {
        int f = 0;
        for (int j = 0; j < ny; j++) if (x[i].a == y[j].a && x[i].b == y[j].b && x[i].v == y[j].v) f = 1;
        if (!f) return 0;
    }
    return 1;
}

int main(int argc, char **argv)
{
    gf_init();
    srand((unsigned)time(NULL));

    int argi = 1, demo = 0;
    if (argc > 1 && !strcmp(argv[1], "demo")) { demo = 1; argi = 2; }

    D = argc > argi ? atoi(argv[argi]) : 8;
    if (D < 2 || D > MAXS - 1) { fprintf(stderr, "D は 2..%d\n", MAXS - 1); return 2; }

    if (demo) {
        int w = argc > argi + 1 ? atoi(argv[argi + 1]) : 4;
        err e[MAXD], f[MAXD];
        random_errors(e, w);
        make_syndromes(e, w);
        printf("既知のシンドローム i+j<%d(%d 個)。誤り %d 個を入れた:\n", D, D * (D + 1) / 2, w);
        for (int k = 0; k < w; k++) printf("  位置 (%2d,%2d)  値 %2d\n", e[k].a, e[k].b, e[k].v);

        int n = decode(f);
        printf("\n求めた誤り位置イデアルの基底(単項式は x^i y^j、係数は GF(16) の整数表現):\n");
        for (int k = 0; k < nb; k++) {
            printf("  先頭項 x^%d y^%d", basis[k].r1, basis[k].r2);
            for (int t = 0; t < dn; t++) if (basis[k].c[t]) printf("  + %d*x^%d y^%d", basis[k].c[t], dl1[t], dl2[t]);
            printf("\n");
        }
        printf("Delta(基底に現れない単項式、= 誤りの個数 %d):", dn);
        for (int t = 0; t < dn; t++) printf(" x^%d y^%d", dl1[t], dl2[t]);
        printf("\n\n復号結果:");
        if (n < 0) printf(" 失敗\n");
        else {
            printf(" %d 個\n", n);
            for (int k = 0; k < n; k++) printf("  位置 (%2d,%2d)  値 %2d\n", f[k].a, f[k].b, f[k].v);
            printf("元の誤りと一致: %s\n", same_set(e, w, f, n) ? "はい" : "いいえ");
        }
        return 0;
    }

    int wmax = argc > argi + 1 ? atoi(argv[argi + 1]) : 10;
    int trials = argc > argi + 2 ? atoi(argv[argi + 2]) : 500;
    printf("GF(16) 上、符号長 15x15、既知のシンドローム i+j < %d の %d 個\n\n", D, D * (D + 1) / 2);
    printf("誤りの個数 | 訂正成功 | 検出失敗 | 誤訂正(別の誤りを返す)\n");
    for (int w = 1; w <= wmax; w++) {
        int succ = 0, fail = 0, mis = 0;
        for (int t = 0; t < trials; t++) {
            err e[MAXD], f[MAXD];
            random_errors(e, w);
            make_syndromes(e, w);
            int n = decode(f);
            if (n < 0) fail++;
            else if (same_set(e, w, f, n)) succ++;
            else mis++;
        }
        printf("    %3d     | %6.1f%%  | %6.1f%%  | %6.1f%%\n", w, 100.0 * succ / trials, 100.0 * fail / trials, 100.0 * mis / trials);
    }
    return 0;
}
