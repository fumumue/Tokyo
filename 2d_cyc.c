/*
 * Demonstration of a two-dimensional cyclic product code over F_17.
 *
 * Each axis uses a length-16, dimension-11 Reed-Solomon cyclic code with
 * generator roots alpha^1,...,alpha^5 (alpha=3).  The product code has
 * parameters [256,121,36].  Encoding is systematic in both dimensions.
 *
 * Decoding is iterative row/column bounded-distance decoding.  Each line
 * decoder searches for error patterns of weight at most two and checks all
 * five syndromes.  This is deliberately NOT the Berlekamp-Massey-Sakata
 * (BMS) algorithm and does NOT guarantee correction of every pattern of 16
 * errors.  It succeeds when the alternating row/column passes can resolve
 * the error pattern; an uncorrectable/ambiguous result is reported.
 *
 * Build: cc -std=c11 -O2 -Wall -Wextra f17_2d_product.c -o f17_2d_product
 * Run:   ./f17_2d_product [number-of-random-errors] [seed]
 */
#include <errno.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define Q 17
#define N 16
#define K 11
#define NSYM 5
#define MAX_ERRORS 16

typedef struct { int v[N][N]; } Grid;

static int mod(int x)
{
    x %= Q;
    return x < 0 ? x + Q : x;
}

static int power(int a, int n)
{
    int r = 1;
    while (n-- > 0) r = mod(r * a);
    return r;
}

static int inverse(int a)
{
    if (mod(a) == 0) return 0;
    return power(mod(a), Q - 2);
}

static void generator(int g[NSYM + 1])
{
    int tmp[NSYM + 1] = {1};
    int degree = 0;
    for (int root = 1; root <= NSYM; ++root) {
        int a = power(3, root);
        memset(&tmp[degree + 1], 0, (size_t)(NSYM - degree) * sizeof(tmp[0]));
        for (int i = degree; i >= 0; --i) {
            tmp[i + 1] = mod(tmp[i + 1] + tmp[i]);
            tmp[i] = mod(tmp[i] - a * tmp[i]);
        }
        ++degree;
    }
    memcpy(g, tmp, sizeof(tmp));
}

/* Systematic cyclic RS encoding: 11 input symbols -> 16 code symbols. */
static void encode_line(const int in[K], int out[N])
{
    int g[NSYM + 1];
    generator(g);
    memset(out, 0, N * sizeof(out[0]));
    for (int i = 0; i < K; ++i) out[i + NSYM] = mod(in[i]);

    for (int d = N - 1; d >= NSYM; --d) {
        int factor = out[d];
        if (factor == 0) continue;
        for (int j = 0; j <= NSYM; ++j)
            out[d - NSYM + j] = mod(out[d - NSYM + j] - factor * g[j]);
    }
    /* Restore the systematic message portion after remainder computation. */
    for (int i = 0; i < K; ++i) out[i + NSYM] = mod(in[i]);
}

static int syndrome(const int line[N], int s[NSYM])
{
    int nonzero = 0;
    for (int j = 1; j <= NSYM; ++j) {
        int sum = 0;
        for (int i = 0; i < N; ++i)
            sum = mod(sum + line[i] * power(power(3, j), i));
        s[j - 1] = sum;
        nonzero |= (sum != 0);
    }
    return nonzero;
}

static int verify_candidate(const int s[NSYM], int p, int e, int q, int f,
                            int count)
{
    for (int j = 1; j <= NSYM; ++j) {
        int predicted = mod(e * power(power(3, j), p));
        if (count == 2) predicted = mod(predicted + f * power(power(3, j), q));
        if (predicted != s[j - 1]) return 0;
    }
    return 1;
}

/* Correct one line if it has a unique error pattern of weight <= 2. */
static int decode_line(int line[N])
{
    int s[NSYM];
    if (!syndrome(line, s)) return 1;

    /* Try a single error. */
    for (int p = 0; p < N; ++p) {
        int a = power(3, p);
        int e = mod(s[0] * inverse(a));
        if (e && verify_candidate(s, p, e, 0, 0, 1)) {
            line[p] = mod(line[p] - e);
            return 1;
        }
    }

    /* Try two errors, solving the first two syndrome equations. */
    for (int p = 0; p < N; ++p) {
        int a = power(3, p), a2 = mod(a * a);
        for (int q = p + 1; q < N; ++q) {
            int b = power(3, q), b2 = mod(b * b);
            int det = mod(a * b2 - b * a2);
            if (det == 0) continue;
            int e = mod((s[0] * b2 - b * s[1]) * inverse(det));
            int f = mod((a * s[1] - s[0] * a2) * inverse(det));
            if (e && f && verify_candidate(s, p, e, q, f, 2)) {
                line[p] = mod(line[p] - e);
                line[q] = mod(line[q] - f);
                return 1;
            }
        }
    }
    return 0;
}

static int line_is_codeword(const int line[N])
{
    int s[NSYM];
    return !syndrome(line, s);
}

static void encode_grid(const int message[K][K], Grid *code)
{
    int temp[K][N];
    for (int r = 0; r < K; ++r) encode_line(message[r], temp[r]);
    for (int c = 0; c < N; ++c) {
        int in[K], out[N];
        for (int r = 0; r < K; ++r) in[r] = temp[r][c];
        encode_line(in, out);
        for (int r = 0; r < N; ++r) code->v[r][c] = out[r];
    }
}

/* Iterate bounded-distance row and column decoding until stable. */
static int decode_grid(Grid *received)
{
    int changed;
    for (int pass = 0; pass < 100; ++pass) {
        changed = 0;
        for (int r = 0; r < N; ++r) {
            int before[N];
            memcpy(before, received->v[r], sizeof(before));
            if (!decode_line(received->v[r])) continue;
            if (memcmp(before, received->v[r], sizeof(before)) != 0) changed = 1;
        }
        for (int c = 0; c < N; ++c) {
            int line[N], before[N];
            for (int r = 0; r < N; ++r) line[r] = received->v[r][c];
            memcpy(before, line, sizeof(before));
            if (!decode_line(line)) continue;
            if (memcmp(before, line, sizeof(before)) != 0) {
                changed = 1;
                for (int r = 0; r < N; ++r) received->v[r][c] = line[r];
            }
        }
        if (!changed) break;
    }

    for (int r = 0; r < N; ++r)
        if (!line_is_codeword(received->v[r])) return 0;
    for (int c = 0; c < N; ++c) {
        int line[N];
        for (int r = 0; r < N; ++r) line[r] = received->v[r][c];
        if (!line_is_codeword(line)) return 0;
    }
    return 1;
}

static int parse_int(const char *text, int min, int max, int *out)
{
    char *end;
    errno = 0;
    long value = strtol(text, &end, 10);
    if (errno || *text == '\0' || *end != '\0' || value < min || value > max)
        return 0;
    *out = (int)value;
    return 1;
}

int main(int argc, char **argv)
{
    int errors = MAX_ERRORS;
    unsigned seed = 1;
    if (argc > 3 || (argc >= 2 && !parse_int(argv[1], 0, MAX_ERRORS, &errors)) ||
        (argc == 3 && !parse_int(argv[2], 0, 2147483647, (int *)&seed))) {
        fprintf(stderr, "Usage: %s [errors:0..16] [seed]\n", argv[0]);
        return 2;
    }
    srand(seed);

    int message[K][K];
    for (int r = 0; r < K; ++r)
        for (int c = 0; c < K; ++c) message[r][c] = rand() % Q;

    Grid code, received;
    encode_grid(message, &code);
    received = code;

    for (int e = 0; e < errors; ++e) {
        int r = rand() % N, c = rand() % N, delta;
        do { delta = rand() % Q; } while (delta == 0);
        received.v[r][c] = mod(received.v[r][c] + delta);
    }

    printf("F_17 2D cyclic product-code demo: 16x16 symbols, %d injected errors\n", errors);
    if (!decode_grid(&received)) {
        puts("Decode failed: residual syndrome remains (this pattern is outside the decoder's guaranteed region).");
        return 1;
    }
    if (memcmp(&received, &code, sizeof(code)) != 0) {
        puts("Decoder reached a valid codeword, but it is not the transmitted codeword.");
        return 1;
    }
    puts("Decoded successfully; recovered codeword matches the transmitted one.");
    return 0;
}

