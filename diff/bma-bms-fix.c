#include <stdio.h>
#include <stdlib.h>

/*
 * Corrected Berlekamp-Massey state update.
 *
 * The original implementation used the same state variable for two distinct
 * meanings:
 *   - L: the current locator length
 *   - m: the last disagreement index
 *
 * This breaks the update rule because the update must use the previous value of
 * m, not the already-updated one. The fixed code keeps them separate.
 */

static OP bma_fixed(unsigned short s[], int kk)
{
    int i, j, k, ll = 0, l, d[2 * K + 1] = {0};
    int last_m = 0;
    OP lo[2 * K + 1] = {0}, b[2 * K + 1] = {0}, t[2 * K + 1] = {0};
    OP a = {0}, f = {0}, h = {0}, g = {0};
    vec v = {0}, x = {0};

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

                /* Correct BM update:
                 * save the previous mismatch position first, then advance.
                 */
                int old_m = last_m;
                last_m = j;
                ll = j - old_m;
            }
        }
    }

    if (odeg(lo[j - 1]) != T)
    {
        fprintf(stderr, "bma: locator degree mismatch\n");
        exit(1);
    }

    return lo[j - 1];
}

/*
 * The same state-separation rule applies to the BMS variant: keep the locator
 * length separate from the disagreement index, do not rewrite the current L with
 * j - last_m until after the old value has been stored.
 */
