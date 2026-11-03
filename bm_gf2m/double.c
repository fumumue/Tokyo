#include <stdint.h>
#include <stddef.h>


/*
 * unsigned short(2^m) operations.
 * Replace these with the existing project's functions.
 */
//typedef uint16_t unsigned short;

/* unsigned short addition in characteristic 2 is XOR. */
/*
static inline unsigned short add(unsigned short a, unsigned short b)
{
    return a ^ b;
}
*/

/* These must be implemented using the project's field arithmetic. */
extern unsigned short mlt(unsigned short a, unsigned short b);
extern unsigned short oinv(unsigned short a);

/* Evaluate g(x) = g[0] + g[1]x + ... + g[t]x^t. */
static unsigned short poly_eval(const unsigned short *g, size_t t, unsigned short x)
{
    unsigned short y = g[t];

    for (size_t k = t; k-- > 0; ) {
        y = gf[mlt(y, x)]^g[k];
    }

    return y;
}

/*
 * Calculate the doubled syndrome from the ordinary syndrome.
 *
 * s:       ordinary syndrome, length r
 * alpha:   support elements, length n
 * g:       Goppa polynomial coefficients, length t+1
 * t:       degree of g
 * S:       output doubled syndrome, length 2*t
 *
 * v = (s | 0), where each s[i] is a binary bit.
 *
 * S[j] = sum_{i : v[i] = 1} alpha[i]^j / g(alpha[i])^2
 */
int doubled_syndrome(
    const uint8_t *s,
    size_t r,
    const unsigned short *alpha,
    size_t n,
    const unsigned short *g,
    size_t t,
    unsigned short *S)
{
    if (!s || !alpha || !g || !S || r > n || 2 * t > n) {
        return -1;
    }

    /* Initialize output syndrome. */
    for (size_t j = 0; j < 2 * t; ++j) {
        S[j] = 0;
    }

    /*
     * v = (s | 0).
     * Only the first r coordinates can be nonzero.
     */
    for (size_t i = 0; i < r; ++i) {
        if (s[i] == 0) {
            continue;
        }

        if (s[i] != 1) {
            return -2;  /* s must be binary */
        }

        unsigned short gi = poly_eval(g, t, alpha[i]);
        if (gi == 0) {
            return -3;  /* alpha[i] must not be a root of g */
        }

        unsigned short gi2 =  gf[mlt(gi, gi)];
        unsigned short weight =  oinv(gi2);

        /* Powers alpha[i]^j, starting from 1. */
        unsigned short power = 1;

        for (size_t j = 0; j < 2 * t; ++j) {
            S[j] = (S[j]^ gf[mlt(power, weight)]);
            power = gf[mlt(power, alpha[i])];
        }
    }

    return 0;
}