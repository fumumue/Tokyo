/*
 * miura_fr.c
 *
 * GF(p^2) 上の三浦(C_ab)曲線による one-point 符号と、
 * Feng-Rao 法(多数決つき・Sakata型の消去)による復号。
 *
 * 曲線:  f(x,y) = y^a + cb*x^b + sum tc * x^ti * y^tj     (gcd(a,b)=1)
 *        ただし各項の重み ti*a + tj*b < a*b。
 *        x の極の位数 a, y の極の位数 b (無限遠の1点 P で)。
 *        基底 x^i y^j (0<=j<a)、極の位数 i*a+j*b、半群 <a,b>、種数 g=(a-1)(b-1)/2。
 *
 * 符号:  C(m) = { c in GF(q)^n : sum_t c_t * phi(P_t) = 0
 *                 (phi は極の位数が m 以下の基底関数すべて) }
 *        P_t はアフィン部分の GF(p^2) 有理点すべて(n 個)。
 *        設計距離は Feng-Rao 限界 d_FR = min{ nu(gamma) : gamma>m, gamma in <a,b> }、
 *        nu(gamma) = #{(l,mu) in <a,b>^2 : l+mu=gamma}。
 *        復号は floor((d_FR-1)/2) 個までの誤りを訂正する。
 *
 * 体の元の表記: 整数 e = u + p*v  が  u + v*theta  (theta^2+s*theta+t=0) を表す。
 *
 * ビルド: gcc -O2 -Wall -o miura_fr miura_fr.c
 * 使い方: ./miura_fr p a b m [-H] [-r seed] [-c cb] [-t i,j,c]... [-n trials] [-e nerr] [-s seed]
 *   -H         エルミート曲線 y^p+y-x^(p+1) (a=p, b=p+1 のとき)
 *   -r seed    下位項を乱数で生成(既定。-t 指定も -H もなければこれ)
 *   -c cb      x^b の係数(元の整数表記、既定は乱数/エルミートでは -1)
 *   -t i,j,c   下位項 c*x^i*y^j を追加(c は元の整数表記)
 *   -n trials  試行回数(既定 200)    -e nerr  誤り数(既定 訂正限界 t_max)
 *   -s seed    メッセージ・誤りの乱数の種
 *
 * 注意: 曲線のアフィン部分が非特異(閉包上で)であることは利用者が保証すること。
 *       プログラムが確認するのは GF(p^2) 有理点での特異点だけ。
 */
#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef uint32_t ff;

/* ================================================================ */
/* GF(p^2)                                                           */
/* ================================================================ */

static int FP, FQ, FS, FT;
static int *cu_, *cv_;
static ff *fexp, *flog;

static ff ff_mul_slow(ff a, ff b)
{
    long long u1 = cu_[a], v1 = cv_[a], u2 = cu_[b], v2 = cv_[b];
    long long u = u1 * u2 - (long long)FT * v1 * v2;
    long long v = u1 * v2 + u2 * v1 - (long long)FS * v1 * v2;
    u %= FP; if (u < 0) u += FP;
    v %= FP; if (v < 0) v += FP;
    return (ff)(u + (long long)FP * v);
}

static void ff_init(int p)
{
    FP = p; FQ = p * p;
    cu_ = malloc(FQ * sizeof(int));
    cv_ = malloc(FQ * sizeof(int));
    for (int e = 0; e < FQ; e++) { cu_[e] = e % p; cv_[e] = e / p; }

    int found = 0;
    for (int s = 0; s < p && !found; s++)
        for (int t = 0; t < p && !found; t++) {
            int ok = 1;
            for (int x = 0; x < p; x++)
                if ((x * x + s * x + t) % p == 0) { ok = 0; break; }
            if (ok) { FS = s; FT = t; found = 1; }
        }

    fexp = malloc(2 * (FQ - 1) * sizeof(ff));
    flog = calloc(FQ, sizeof(ff));
    for (ff g = 2; g < (ff)FQ; g++) {
        ff x = 1;
        int k;
        for (k = 0; k < FQ - 1; k++) {
            fexp[k] = x;
            x = ff_mul_slow(x, g);
            if (x == 1) { k++; break; }
        }
        if (k == FQ - 1 && x == 1) {
            for (int i = 0; i < FQ - 1; i++) {
                fexp[i + FQ - 1] = fexp[i];
                flog[fexp[i]] = (ff)i;
            }
            return;
        }
    }
    fprintf(stderr, "ff_init: no primitive element\n");
    exit(1);
}

static inline ff ff_add(ff a, ff b)
{
    int u = cu_[a] + cu_[b], v = cv_[a] + cv_[b];
    if (u >= FP) u -= FP;
    if (v >= FP) v -= FP;
    return (ff)(u + FP * v);
}
static inline ff ff_neg(ff a)
{
    int u = cu_[a], v = cv_[a];
    u = u ? FP - u : 0;
    v = v ? FP - v : 0;
    return (ff)(u + FP * v);
}
static inline ff ff_sub(ff a, ff b) { return ff_add(a, ff_neg(b)); }
static inline ff ff_mul(ff a, ff b)
{
    if (!a || !b) return 0;
    return fexp[flog[a] + flog[b]];
}
static inline ff ff_inv(ff a) { return fexp[(FQ - 1) - flog[a]]; }
static ff ff_pow(ff x, int k)
{
    if (k == 0) return 1;
    if (!x) return 0;
    return fexp[(int)(((long long)flog[x] * k) % (FQ - 1))];
}
static ff ff_scal(int k, ff a) /* 整数 k 倍 */
{
    k %= FP; if (k < 0) k += FP;
    return (ff)((k * cu_[a]) % FP + FP * ((k * cv_[a]) % FP));
}

static uint64_t rs = 88172645463325252ULL;
static uint64_t rng(void)
{
    rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17;
    return rs;
}
static ff rff(void) { return (ff)(rng() % FQ); }
static ff rff_nz(void) { return (ff)(1 + rng() % (FQ - 1)); }

/* ================================================================ */
/* 曲線                                                              */
/* ================================================================ */

#define MAXT 512
typedef struct {
    int p, a, b, g;
    ff cb;
    int nt;
    int ti[MAXT], tj[MAXT];
    ff tc[MAXT];
    int n;
    ff *px, *py;
} curve;

static int gcd_int(int x, int y) { return y ? gcd_int(y, x % y) : x; }

/* アフィン有理点を列挙。有理点に特異点があれば -1。 */
static int curve_points(curve *C)
{
    int a = C->a, b = C->b;
    C->g = (a - 1) * (b - 1) / 2;
    C->px = malloc((size_t)FQ * a * sizeof(ff) + sizeof(ff));
    C->py = malloc((size_t)FQ * a * sizeof(ff) + sizeof(ff));
    C->n = 0;
    ff *xp = malloc((b + 1) * sizeof(ff));
    ff *cj = malloc(a * sizeof(ff));
    for (int x = 0; x < FQ; x++) {
        xp[0] = 1;
        for (int k = 1; k <= b; k++) xp[k] = ff_mul(xp[k - 1], (ff)x);
        memset(cj, 0, a * sizeof(ff));
        cj[0] = ff_add(cj[0], ff_mul(C->cb, xp[b]));
        for (int k = 0; k < C->nt; k++)
            cj[C->tj[k]] = ff_add(cj[C->tj[k]], ff_mul(C->tc[k], xp[C->ti[k]]));
        for (int y = 0; y < FQ; y++) {
            ff val = 1;
            for (int j = a - 1; j >= 0; j--) val = ff_add(ff_mul(val, (ff)y), cj[j]);
            if (val) continue;
            /* 特異点チェック */
            ff fx = ff_mul(ff_scal(b, C->cb), ff_pow((ff)x, b - 1));
            ff fy = ff_mul(ff_scal(a, 1), ff_pow((ff)y, a - 1));
            for (int k = 0; k < C->nt; k++) {
                int ti = C->ti[k], tj = C->tj[k];
                if (ti > 0)
                    fx = ff_add(fx, ff_mul(ff_scal(ti, C->tc[k]),
                                           ff_mul(ff_pow((ff)x, ti - 1), ff_pow((ff)y, tj))));
                if (tj > 0)
                    fy = ff_add(fy, ff_mul(ff_scal(tj, C->tc[k]),
                                           ff_mul(ff_pow((ff)x, ti), ff_pow((ff)y, tj - 1))));
            }
            if (!fx && !fy) { free(xp); free(cj); return -1; }
            C->px[C->n] = (ff)x;
            C->py[C->n] = (ff)y;
            C->n++;
        }
    }
    free(xp); free(cj);
    return 0;
}

static void curve_random_terms(curve *C)
{
    C->nt = 0;
    int a = C->a, b = C->b;
    for (int j = 0; j < a; j++)
        for (int i = 0; i * a + j * b < a * b; i++) {
            if (i == 0 && j == 0 && 0) continue;
            if (rng() % 3 == 0) {
                if (C->nt >= MAXT) return;
                C->ti[C->nt] = i; C->tj[C->nt] = j; C->tc[C->nt] = rff_nz();
                C->nt++;
            }
        }
}

static void curve_print(const curve *C)
{
    printf("curve: y^%d + [%u] x^%d", C->a, C->cb, C->b);
    for (int k = 0; k < C->nt; k++) printf(" + [%u] x^%d y^%d", C->tc[k], C->ti[k], C->tj[k]);
    printf(" = 0   over GF(%d^2), element e = u+%d*v means u+v*theta, theta^2+%d*theta+%d=0\n",
           C->p, FP, FS, FT);
}

/* ================================================================ */
/* Feng-Rao 復号                                                      */
/* ================================================================ */

typedef struct {
    const curve *C;
    int m, G, nL, ell, dFR, tmax, lamuse;
    int *lam, *bk, *bl, *idxw;
    int kmax, nmon, *monid;
    ff *nf;  /* nf[mon*nL + idx] */
    ff *phi; /* phi[idx*n + t]   */
} fr_ctx;

static int fr_init(fr_ctx *X, const curve *C, int m)
{
    int a = C->a, b = C->b, g = C->g;
    memset(X, 0, sizeof *X);
    X->C = C;
    X->m = m;
    if (m < 0 || m >= C->n) return -1;

    int W = m + 4 * g + a + b + 4;
    char *inL = calloc(W + 1, 1);
    for (int j = 0; j < a; j++)
        for (int i = 0; i * a + j * b <= W; i++) inL[i * a + j * b] = 1;

    int dFR = INT_MAX;
    for (int gam = m + 1; gam <= W; gam++) {
        if (!inL[gam]) continue;
        int nu = 0;
        for (int l = 0; l <= gam; l++)
            if (inL[l] && inL[gam - l]) nu++;
        if (nu < dFR) dFR = nu;
    }
    X->dFR = dFR;
    X->tmax = (dFR - 1) / 2;
    int G = 2 * X->tmax + 4 * g + 2;
    if (G < m) G = m;
    if (G > W) G = W;
    X->G = G;
    int lu = G - (X->tmax + 2 * g - 1);
    X->lamuse = lu < 0 ? G : lu;

    /* 半群の元と基底単項式 */
    X->idxw = malloc((G + 1) * sizeof(int));
    int nL = 0;
    for (int w = 0; w <= G; w++) { X->idxw[w] = inL[w] ? nL++ : -1; }
    X->nL = nL;
    X->lam = malloc(nL * sizeof(int));
    X->bk = malloc(nL * sizeof(int));
    X->bl = malloc(nL * sizeof(int));
    for (int w = 0; w <= G; w++) {
        if (!inL[w]) continue;
        int id = X->idxw[w];
        X->lam[id] = w;
        for (int j = 0; j < a; j++) {
            int r = w - j * b;
            if (r >= 0 && r % a == 0) { X->bk[id] = r / a; X->bl[id] = j; break; }
        }
    }
    X->ell = 0;
    while (X->ell < nL && X->lam[X->ell] <= m) X->ell++;

    /* 単項式 x^k y^l (l<=2a-2, 重み<=G) の正規形 */
    int L2 = 2 * a - 1;
    int kmax = G / a;
    X->kmax = kmax;
    X->monid = malloc((size_t)(kmax + 1) * L2 * sizeof(int));
    for (size_t i = 0; i < (size_t)(kmax + 1) * L2; i++) X->monid[i] = -1;
    int nmon = 0;
    for (int w = 0; w <= G; w++)
        for (int l = 0; l < L2; l++) {
            int r = w - l * b;
            if (r >= 0 && r % a == 0 && r / a <= kmax) nmon++;
        }
    X->nmon = nmon;
    X->nf = calloc((size_t)nmon * nL, sizeof(ff));
    int cnt = 0;
    for (int w = 0; w <= G; w++)
        for (int l = 0; l < L2; l++) {
            int r = w - l * b;
            if (!(r >= 0 && r % a == 0 && r / a <= kmax)) continue;
            int k = r / a;
            int id = cnt++;
            X->monid[k * L2 + l] = id;
            ff *vec = X->nf + (size_t)id * nL;
            if (l < a) {
                vec[X->idxw[w]] = 1;
            } else {
                /* y^a = -cb x^b - sum tc x^ti y^tj */
                int id1 = X->monid[(k + b) * L2 + (l - a)];
                const ff *v1 = X->nf + (size_t)id1 * nL;
                ff mc = ff_neg(C->cb);
                for (int i = 0; i < nL; i++) vec[i] = ff_mul(mc, v1[i]);
                for (int q = 0; q < C->nt; q++) {
                    int id2 = X->monid[(k + C->ti[q]) * L2 + (l - a + C->tj[q])];
                    const ff *v2 = X->nf + (size_t)id2 * nL;
                    ff mt = ff_neg(C->tc[q]);
                    for (int i = 0; i < nL; i++)
                        if (v2[i]) vec[i] = ff_add(vec[i], ff_mul(mt, v2[i]));
                }
            }
        }

    /* 基底関数の評価値 */
    int n = C->n;
    X->phi = malloc((size_t)nL * n * sizeof(ff));
    for (int i = 0; i < nL; i++)
        for (int t = 0; t < n; t++)
            X->phi[(size_t)i * n + t] =
                ff_mul(ff_pow(C->px[t], X->bk[i]), ff_pow(C->py[t], X->bl[i]));
    free(inL);
    return 0;
}

/* 符号語を1つ乱数で作る(検査行列の零空間) */
static void fr_random_codeword(const fr_ctx *X, ff *c)
{
    int n = X->C->n, ell = X->ell;
    ff *H = malloc((size_t)ell * n * sizeof(ff));
    memcpy(H, X->phi, (size_t)ell * n * sizeof(ff));
    int *pc = malloc(ell * sizeof(int));
    char *isp = calloc(n, 1);
    int r = 0;
    for (int col = 0; col < n && r < ell; col++) {
        int pr = -1;
        for (int i = r; i < ell; i++)
            if (H[(size_t)i * n + col]) { pr = i; break; }
        if (pr < 0) continue;
        if (pr != r)
            for (int j = 0; j < n; j++) {
                ff tmp = H[(size_t)pr * n + j];
                H[(size_t)pr * n + j] = H[(size_t)r * n + j];
                H[(size_t)r * n + j] = tmp;
            }
        ff iv = ff_inv(H[(size_t)r * n + col]);
        for (int j = 0; j < n; j++) H[(size_t)r * n + j] = ff_mul(H[(size_t)r * n + j], iv);
        for (int i = 0; i < ell; i++) {
            if (i == r) continue;
            ff f = H[(size_t)i * n + col];
            if (!f) continue;
            for (int j = 0; j < n; j++)
                H[(size_t)i * n + j] = ff_sub(H[(size_t)i * n + j], ff_mul(f, H[(size_t)r * n + j]));
        }
        pc[r] = col; isp[col] = 1; r++;
    }
    for (int j = 0; j < n; j++) c[j] = isp[j] ? 0 : rff();
    for (int i = 0; i < r; i++) {
        ff acc = 0;
        for (int j = 0; j < n; j++)
            if (!isp[j] && c[j]) acc = ff_add(acc, ff_mul(H[(size_t)i * n + j], c[j]));
        c[pc[i]] = ff_neg(acc);
    }
    free(H); free(pc); free(isp);
}

/*
 * 受信語 r から誤りベクトル e_out を求める。
 * 戻り値 0: 成功(*nerr に誤り数), 非0: 復号失敗
 */
static int fr_decode(const fr_ctx *X, const ff *r, ff *e_out, int *nerr)
{
    const curve *C = X->C;
    int n = C->n, nL = X->nL, m = X->m, a = C->a, L2 = 2 * a - 1;
    int ret = 0;

    ff *s = calloc(nL, sizeof(ff));
    for (int i = 0; i < nL && X->lam[i] <= m; i++) {
        ff acc = 0;
        for (int t = 0; t < n; t++)
            if (r[t]) acc = ff_add(acc, ff_mul(r[t], X->phi[(size_t)i * n + t]));
        s[i] = acc;
    }

    ff *S = calloc((size_t)nL * nL, sizeof(ff));
    ff *F = calloc((size_t)nL * nL, sizeof(ff));
    for (int i = 0; i < nL; i++) F[(size_t)i * nL + i] = 1;
    char *piv = calloc(nL, 1);
    ff *dpiv = calloc(nL, sizeof(ff));
    int *rowofcol = malloc(nL * sizeof(int));
    for (int i = 0; i < nL; i++) rowofcol[i] = -1;
    ff *pK = malloc(nL * sizeof(ff)), *pC = malloc(nL * sizeof(ff)), *votes = malloc(nL * sizeof(ff));
    int *pM = malloc(nL * sizeof(int));
    int npiv = 0, fail = 0;

    for (int gi = 0; gi < nL; gi++) {
        int gam = X->lam[gi];
        for (int li = 0; li <= gi; li++) {
            int w = gam - X->lam[li];
            int mi = (w >= 0 && w <= X->G) ? X->idxw[w] : -1;
            pM[li] = mi;
            if (mi < 0) continue;
            int mon = X->monid[(X->bk[li] + X->bk[mi]) * L2 + X->bl[li] + X->bl[mi]];
            const ff *vec = X->nf + (size_t)mon * nL;
            pC[li] = vec[gi];
            ff K = 0;
            for (int t = 0; t < gi; t++)
                if (vec[t] && s[t]) K = ff_add(K, ff_mul(vec[t], s[t]));
            pK[li] = K;
        }

        if (gam > m) { /* 未知のシンドローム: 多数決 */
            int nv = 0;
            for (int li = 0; li <= gi; li++) {
                int mi = pM[li];
                if (mi < 0 || piv[li]) continue;
                ff D = pK[li];
                const ff *Fl = F + (size_t)li * nL;
                for (int k = 0; k < li; k++)
                    if (Fl[k]) D = ff_add(D, ff_mul(Fl[k], S[(size_t)k * nL + mi]));
                votes[nv++] = ff_neg(ff_mul(D, ff_inv(pC[li])));
            }
            if (nv == 0) { fail = 1; s[gi] = 0; }
            else {
                int bestc = 0; ff best = votes[0];
                for (int i = 0; i < nv; i++) {
                    int c = 0;
                    for (int j = 0; j < nv; j++) if (votes[j] == votes[i]) c++;
                    if (c > bestc) { bestc = c; best = votes[i]; }
                }
                s[gi] = best;
            }
        }

        for (int li = 0; li <= gi; li++) {
            int mi = pM[li];
            if (mi < 0) continue;
            ff val = ff_add(pK[li], ff_mul(pC[li], s[gi]));
            S[(size_t)li * nL + mi] = val;
            if (piv[li]) continue;
            ff d = val;
            ff *Fl = F + (size_t)li * nL;
            for (int k = 0; k < li; k++)
                if (Fl[k]) d = ff_add(d, ff_mul(Fl[k], S[(size_t)k * nL + mi]));
            if (!d) continue;
            int kp = rowofcol[mi];
            if (kp >= 0) {
                ff tt = ff_mul(d, ff_inv(dpiv[kp]));
                const ff *Fk = F + (size_t)kp * nL;
                for (int k = 0; k <= kp; k++)
                    if (Fk[k]) Fl[k] = ff_sub(Fl[k], ff_mul(tt, Fk[k]));
            } else {
                piv[li] = 1; dpiv[li] = d; rowofcol[mi] = li; npiv++;
            }
        }
    }

    memset(e_out, 0, n * sizeof(ff));
    *nerr = 0;
    if (fail || npiv > X->tmax) { ret = 1; goto done; }

    /* 誤り位置: 従属な行の関係式がすべて 0 になる点 */
    {
        char *alive = malloc(n);
        memset(alive, 1, n);
        for (int li = 0; li < nL; li++) {
            if (piv[li] || X->lam[li] > X->lamuse) continue;
            const ff *Fl = F + (size_t)li * nL;
            for (int t = 0; t < n; t++) {
                if (!alive[t]) continue;
                ff v = 0;
                for (int k = 0; k <= li; k++)
                    if (Fl[k]) v = ff_add(v, ff_mul(Fl[k], X->phi[(size_t)k * n + t]));
                if (v) alive[t] = 0;
            }
        }
        int ne = 0;
        int *pos = malloc(n * sizeof(int));
        for (int t = 0; t < n; t++) if (alive[t]) pos[ne++] = t;
        free(alive);

        if (ne > X->ell || ne > X->tmax) { ret = 2; free(pos); goto done; }

        /* 誤りの値: sum_t e_t phi_i(P_t) = s_i (lam_i <= m) */
        int rows = X->ell, cols = ne + 1;
        ff *A = malloc((size_t)rows * cols * sizeof(ff));
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < ne; j++) A[(size_t)i * cols + j] = X->phi[(size_t)i * n + pos[j]];
            A[(size_t)i * cols + ne] = s[i];
        }
        int rk = 0;
        for (int col = 0; col < ne; col++) {
            int pr = -1;
            for (int i = rk; i < rows; i++) if (A[(size_t)i * cols + col]) { pr = i; break; }
            if (pr < 0) { ret = 3; break; }
            if (pr != rk)
                for (int j = 0; j < cols; j++) {
                    ff tmp = A[(size_t)pr * cols + j];
                    A[(size_t)pr * cols + j] = A[(size_t)rk * cols + j];
                    A[(size_t)rk * cols + j] = tmp;
                }
            ff iv = ff_inv(A[(size_t)rk * cols + col]);
            for (int j = 0; j < cols; j++) A[(size_t)rk * cols + j] = ff_mul(A[(size_t)rk * cols + j], iv);
            for (int i = 0; i < rows; i++) {
                if (i == rk) continue;
                ff f = A[(size_t)i * cols + col];
                if (!f) continue;
                for (int j = 0; j < cols; j++)
                    A[(size_t)i * cols + j] = ff_sub(A[(size_t)i * cols + j], ff_mul(f, A[(size_t)rk * cols + j]));
            }
            rk++;
        }
        if (!ret) {
            for (int i = rk; i < rows; i++)
                if (A[(size_t)i * cols + ne]) { ret = 4; break; }
        }
        if (!ret) {
            for (int j = 0; j < ne; j++) {
                e_out[pos[j]] = A[(size_t)j * cols + ne];
                if (e_out[pos[j]]) (*nerr)++;
            }
        }
        free(A); free(pos);
    }

done:
    free(s); free(S); free(F); free(piv); free(dpiv); free(rowofcol);
    free(pK); free(pC); free(votes); free(pM);
    return ret;
}

/* ================================================================ */
/* テスト用 main                                                      */
/* ================================================================ */

#ifndef MIURA_NO_MAIN
static void usage(void)
{
    fprintf(stderr,
        "usage: miura_fr p a b m [-H] [-r seed] [-c cb] [-t i,j,c]... [-n trials] [-e nerr] [-s seed]\n");
    exit(1);
}

int main(int argc, char **argv)
{
    if (argc < 5) usage();
    curve C;
    memset(&C, 0, sizeof C);
    C.p = atoi(argv[1]); C.a = atoi(argv[2]); C.b = atoi(argv[3]);
    int m = atoi(argv[4]);
    int herm = 0, trials = 200, nerr_opt = -1, have_t = 0, have_c = 0;
    unsigned long cseed = 1, mseed = 12345;
    ff cb_opt = 0;
    struct { int i, j; ff c; } tt[MAXT];
    int ntt = 0;
    for (int k = 5; k < argc; k++) {
        if (!strcmp(argv[k], "-H")) herm = 1;
        else if (!strcmp(argv[k], "-r") && k + 1 < argc) cseed = strtoul(argv[++k], 0, 10);
        else if (!strcmp(argv[k], "-c") && k + 1 < argc) { cb_opt = (ff)atoi(argv[++k]); have_c = 1; }
        else if (!strcmp(argv[k], "-t") && k + 1 < argc) {
            int i, j, c;
            if (sscanf(argv[++k], "%d,%d,%d", &i, &j, &c) != 3 || ntt >= MAXT) usage();
            tt[ntt].i = i; tt[ntt].j = j; tt[ntt].c = (ff)c; ntt++; have_t = 1;
        }
        else if (!strcmp(argv[k], "-n") && k + 1 < argc) trials = atoi(argv[++k]);
        else if (!strcmp(argv[k], "-e") && k + 1 < argc) nerr_opt = atoi(argv[++k]);
        else if (!strcmp(argv[k], "-s") && k + 1 < argc) mseed = strtoul(argv[++k], 0, 10);
        else usage();
    }
    int p = C.p;
    if (p < 2 || C.a < 2 || C.b < 2 || gcd_int(C.a, C.b) != 1) {
        fprintf(stderr, "need prime p, a,b >= 2, gcd(a,b)=1\n"); return 1;
    }
    for (int d = 2; d * d <= p; d++) if (p % d == 0) { fprintf(stderr, "p must be prime\n"); return 1; }
    ff_init(p);

    int ok_curve = 0;
    if (herm) {
        if (C.a != p || C.b != p + 1) { fprintf(stderr, "-H needs a=p, b=p+1\n"); return 1; }
        C.cb = ff_neg(1);
        C.nt = 1; C.ti[0] = 0; C.tj[0] = 1; C.tc[0] = 1;
        ok_curve = (curve_points(&C) == 0);
    } else if (have_t || have_c) {
        C.cb = have_c ? cb_opt : ff_neg(1);
        C.nt = ntt;
        for (int k = 0; k < ntt; k++) {
            if (tt[k].i * C.a + tt[k].j * C.b >= C.a * C.b) {
                fprintf(stderr, "term %d,%d has weight >= ab\n", tt[k].i, tt[k].j); return 1;
            }
            C.ti[k] = tt[k].i; C.tj[k] = tt[k].j; C.tc[k] = tt[k].c;
        }
        if (!C.cb) { fprintf(stderr, "cb must be nonzero\n"); return 1; }
        ok_curve = (curve_points(&C) == 0);
    } else {
        rs = 0x9E3779B97F4A7C15ULL ^ (cseed * 2654435761UL + 1);
        for (int tries = 0; tries < 2000 && !ok_curve; tries++) {
            C.cb = rff_nz();
            curve_random_terms(&C);
            free(C.px); free(C.py);
            ok_curve = (curve_points(&C) == 0);
        }
    }
    if (!ok_curve) { fprintf(stderr, "curve has a singular rational point\n"); return 1; }

    curve_print(&C);
    printf("q=%d  a=%d b=%d  genus g=%d  affine points n=%d\n", FQ, C.a, C.b, C.g, C.n);

    fr_ctx X;
    if (fr_init(&X, &C, m)) { fprintf(stderr, "bad m (need 0 <= m < n)\n"); return 1; }
    printf("m=%d  dim k=%d  designed d>=%d  Feng-Rao d_FR=%d  t_max=%d  (Goppa-bound t=%d)\n",
           m, C.n - X.ell, m + 2 - 2 * C.g, X.dFR, X.tmax, (m + 1 - 2 * C.g) / 2);

    int nerr = nerr_opt >= 0 ? nerr_opt : X.tmax;
    if (nerr > C.n) nerr = C.n;
    rs = 0x1234567ULL ^ (mseed * 11400714819323198485ULL + 7);
    int okc = 0;
    ff *c = malloc(C.n * sizeof(ff)), *e = malloc(C.n * sizeof(ff)), *rr = malloc(C.n * sizeof(ff));
    ff *eo = malloc(C.n * sizeof(ff));
    int *perm = malloc(C.n * sizeof(int));
    for (int tr = 0; tr < trials; tr++) {
        fr_random_codeword(&X, c);
        memset(e, 0, C.n * sizeof(ff));
        for (int i = 0; i < C.n; i++) perm[i] = i;
        for (int i = 0; i < nerr; i++) {
            int j = i + (int)(rng() % (C.n - i));
            int tmp = perm[i]; perm[i] = perm[j]; perm[j] = tmp;
            e[perm[i]] = rff_nz();
        }
        for (int i = 0; i < C.n; i++) rr[i] = ff_add(c[i], e[i]);
        int ne = 0;
        int ret = fr_decode(&X, rr, eo, &ne);
        if (!ret && !memcmp(eo, e, C.n * sizeof(ff))) okc++;
    }
    printf("errors=%d : decoded correctly %d / %d\n", nerr, okc, trials);
    return okc == trials ? 0 : 2;
}
#endif