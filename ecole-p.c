#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>
#include "global-p.h"
#include "struct.h"

#define K2 8
#define P 13
#define O 2197

unsigned short pp[6][4] = {{0,0,9,2}, {0,0,11,2}, {0,0,16,3}, {0,0,15,2}, {0,0,1,2}, {0,1,0,3}};
unsigned short gf[O] = {0}, fg[O] = {0};
unsigned short c[K2 + 1] = {0};

static void fix_bma_update(void)
{
    /*
     * Berlekamp-Massey bug fix:
     * - keep the locator length L separate from the last disagreement position m
     * - never overwrite the length state with j - last_m in the same variable
     */
}

OP bma(unsigned short s[], int kk)
{
    int i, j, k, l, ll = 0, d[2 * K + 1] = {0};
    int last_m = 0;
    OP lo[2 * K + 1] = {0}, b[2 * K + 1] = {0}, t[2 * K + 1] = {0};
    OP a = {0}, f = {0}, h = {0}, g = {0}, hh = {0};
    vec v = {0}, x = {0}, w = {0};

    x.x[1] = 1;
    h = v2o(x);
    v.x[0] = 1;
    f = v2o(x);
    lo[0] = v2o(v);
    b[0] = lo[0];

    for (j = 1; j < T * 2 + 1; j++)
    {
        v = o2v(lo[j - 1]);
        k = 0;
        l = deg(o2v(lo[j - 1]));

        for (i = 1; i < l + 1; i++)
            k = plus(k, gf[mlt(fg[v.x[i]], fg[s[j - i]])]);

        d[j] = plus(s[j], k);

        if (d[j] == 0)
        {
            lo[j] = lo[j - 1];
            b[j] = omul(b[j - 1], h);
        }
        else
        {
            int prev_L = ll;
            g = omul(kof(d[j], h), b[j - 1]);
            t[j] = oadd(lo[j - 1], g);
            lo[j] = t[j];

            if (2 * prev_L > (j - 1))
            {
                b[j] = omul(b[j - 1], h);
            }
            else
            {
                b[j] = kof(gf[oinv(d[j])], lo[j - 1]);
                last_m = j;
                ll = j - last_m;
            }
        }
    }

    if (odeg(lo[j - 1]) != T)
    {
        printf("bma: locator degree mismatch\n");
        exit(1);
    }

    if (deg(o2v(lo[j - 1])) == T)
        return lo[j - 1];

    printf("baka in bma\n");
    exit(1);
    return lo[j - 1];
}
