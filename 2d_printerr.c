/*
 * F_17 two-dimensional cyclic product-code demonstration.
 *
 * A length-16, dimension-11 cyclic Reed-Solomon code (roots alpha^1..alpha^5,
 * alpha=3) is applied along both axes, giving a [256,121,36] product code.
 * This demo uses iterative row/column bounded-distance decoding up to two
 * errors per line. It is not the Berlekamp-Massey-Sakata algorithm and does
 * not guarantee correction of every pattern of 16 errors.
 *
 * Prints the original codeword, received word before correction, error vector
 * (received - original, over F_17), decoded word, and residual error vector.
 * Build: cc -std=c11 -O2 -Wall -Wextra f17_2d_product.c -o f17_2d_product
 * Run:   ./f17_2d_product [errors 0..16] [seed]
 */
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define Q 17
#define N 16
#define K 11
#define NSYM 5
#define MAX_ERRORS 16

typedef struct { int v[N][N]; } Grid;

static int mod(int x) { x %= Q; return x < 0 ? x + Q : x; }
static int pow17(int a, int n)
{
    int r = 1;
    while (n-- > 0) r = mod(r * a);
    return r;
}
static int inv17(int a) { return a ? pow17(mod(a), Q - 2) : 0; }

/* Coefficients are stored in ascending powers. */
static void make_generator(int g[NSYM + 1])
{
    int poly[NSYM + 1] = {1};
    int degree = 0;
    for (int root = 1; root <= NSYM; ++root) {
        int a = pow17(3, root);
        for (int i = degree; i >= 0; --i) {
            poly[i + 1] = mod(poly[i + 1] + poly[i]);
            poly[i] = mod(poly[i] - a * poly[i]);
        }
        ++degree;
    }
    memcpy(g, poly, sizeof(poly));
}

static void encode_line(const int in[K], int out[N])
{
    int g[NSYM + 1];
    make_generator(g);
    memset(out, 0, N * sizeof(*out));
    for (int i = 0; i < K; ++i) out[i + NSYM] = mod(in[i]);
    for (int d = N - 1; d >= NSYM; --d) {
        int factor = out[d];
        for (int j = 0; j <= NSYM; ++j)
            out[d - NSYM + j] = mod(out[d - NSYM + j] - factor * g[j]);
    }
    for (int i = 0; i < K; ++i) out[i + NSYM] = mod(in[i]);
}

static int get_syndromes(const int line[N], int s[NSYM])
{
    int any = 0;
    for (int j = 1; j <= NSYM; ++j) {
        int sum = 0, x = pow17(3, j), term = 1;
        for (int i = 0; i < N; ++i) {
            sum = mod(sum + line[i] * term);
            term = mod(term * x);
        }
        s[j - 1] = sum;
        any |= sum != 0;
    }
    return any;
}

static int matches(const int s[NSYM], int p, int e, int q, int f, int count)
{
    for (int j = 1; j <= NSYM; ++j) {
        int x = pow17(3, j);
        int predicted = mod(e * pow17(x, p));
        if (count == 2) predicted = mod(predicted + f * pow17(x, q));
        if (predicted != s[j - 1]) return 0;
    }
    return 1;
}

/* Decode one line when there is a unique correctable pattern of weight <= 2. */
static int decode_line(int line[N])
{
    int s[NSYM];
    if (!get_syndromes(line, s)) return 1;

    for (int p = 0; p < N; ++p) {
        int a = pow17(3, p);
        int e = mod(s[0] * inv17(a));
        if (e && matches(s, p, e, 0, 0, 1)) {
            line[p] = mod(line[p] - e);
            return 1;
        }
    }
    for (int p = 0; p < N; ++p) {
        int a = pow17(3, p), a2 = mod(a * a);
        for (int q = p + 1; q < N; ++q) {
            int b = pow17(3, q), b2 = mod(b * b);
            int det = mod(a * b2 - b * a2);
            if (!det) continue;
            int e = mod((s[0] * b2 - b * s[1]) * inv17(det));
            int f = mod((a * s[1] - s[0] * a2) * inv17(det));
            if (e && f && matches(s, p, e, q, f, 2)) {
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
    return !get_syndromes(line, s);
}

static void encode_grid(const int msg[K][K], Grid *out)
{
    int temp[K][N];
    for (int r = 0; r < K; ++r) encode_line(msg[r], temp[r]);
    for (int c = 0; c < N; ++c) {
        int in[K], encoded[N];
        for (int r = 0; r < K; ++r) in[r] = temp[r][c];
        encode_line(in, encoded);
        for (int r = 0; r < N; ++r) out->v[r][c] = encoded[r];
    }
}

static int decode_grid(Grid *g)
{
    for (int pass = 0; pass < 100; ++pass) {
        int changed = 0;
        for (int r = 0; r < N; ++r) {
            int before[N];
            memcpy(before, g->v[r], sizeof(before));
            (void)decode_line(g->v[r]);
            if (memcmp(before, g->v[r], sizeof(before))) changed = 1;
        }
        for (int c = 0; c < N; ++c) {
            int line[N], before[N];
            for (int r = 0; r < N; ++r) line[r] = g->v[r][c];
            memcpy(before, line, sizeof(before));
            (void)decode_line(line);
            if (memcmp(before, line, sizeof(before))) {
                changed = 1;
                for (int r = 0; r < N; ++r) g->v[r][c] = line[r];
            }
        }
        if (!changed) break;
    }
    for (int r = 0; r < N; ++r)
        if (!line_is_codeword(g->v[r])) return 0;
    for (int c = 0; c < N; ++c) {
        int line[N];
        for (int r = 0; r < N; ++r) line[r] = g->v[r][c];
        if (!line_is_codeword(line)) return 0;
    }
    return 1;
}

static void print_grid(const char *title, const Grid *g)
{
    printf("\n%s\n", title);
    for (int r = 0; r < N; ++r) {
        for (int c = 0; c < N; ++c) printf("%2d%c", g->v[r][c], c == N - 1 ? '\n' : ' ');
    }
}

static int parse_int(const char *s, long min, long max, long *out)
{
    char *end;
    errno = 0;
    long v = strtol(s, &end, 10);
    if (errno || *s == '\0' || *end != '\0' || v < min || v > max) return 0;
    *out = v;
    return 1;
}

int main(int argc, char **argv)
{
    long errors = MAX_ERRORS, seed = 1;
    if (argc > 3 || (argc > 1 && !parse_int(argv[1], 0, MAX_ERRORS, &errors)) ||
        (argc > 2 && !parse_int(argv[2], 0, 2147483647L, &seed))) {
        fprintf(stderr, "Usage: %s [errors:0..16] [seed]\n", argv[0]);
        return 2;
    }
    srand((unsigned)seed);

    int message[K][K];
    for (int r = 0; r < K; ++r)
        for (int c = 0; c < K; ++c) message[r][c] = rand() % Q;

    Grid original, received, corrected, error, residual;
    encode_grid(message, &original);
    received = original;

    /* Inject exactly the requested number of errors at distinct positions. */
    int used[N][N] = {{0}};
    for (int i = 0; i < errors; ++i) {
        int r, c, delta;
        do { r = rand() % N; c = rand() % N; } while (used[r][c]);
        used[r][c] = 1;
        do { delta = rand() % Q; } while (!delta);
        received.v[r][c] = mod(received.v[r][c] + delta);
    }

    for (int r = 0; r < N; ++r)
        for (int c = 0; c < N; ++c)
            error.v[r][c] = mod(received.v[r][c] - original.v[r][c]);

    corrected = received;
    int success = decode_grid(&corrected);
    for (int r = 0; r < N; ++r)
        for (int c = 0; c < N; ++c)
            residual.v[r][c] = mod(corrected.v[r][c] - original.v[r][c]);

    printf("F_17 two-dimensional cyclic product code; 16x16 symbols; errors=%ld; seed=%ld\n",
           errors, seed);
    print_grid("Original codeword", &original);
    print_grid("Received word before correction", &received);
    print_grid("Error vector (received - original, mod 17)", &error);
    print_grid("Decoded word after correction", &corrected);
    print_grid("Residual error vector (decoded - original, mod 17)", &residual);

    if (!success) {
        puts("Decode status: failed (nonzero row/column syndrome remains).");
        return 1;
    }
    if (memcmp(&corrected, &original, sizeof(original)) != 0) {
        puts("Decode status: reached a codeword, but not the transmitted codeword.");
        return 1;
    }
    puts("Decode status: success; decoded word equals the original codeword.");
    return 0;
}
