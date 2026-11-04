#include <stdio.h>
#include <stdint.h>

#include "global.h"
#include "struct.h"

void cyclic_convolution(
    const vec a,
    const vec b,
    vec *c,
    int n
) {
    for (int k = 0; k < n; k++) {
        c->x[k] = 0;

        for (int i = 0; i < n; i++) {
            c->x[k] ^= (a.x[i] & b.x[(k - i + n) % n]);
        }
    }
}

int main(void) {
    int n = 8;

    vec a = {1,1,0,1,0,0,1,0};
    vec b = {1,0,1,0,1,0,0,0};
    vec c={0};

    cyclic_convolution(a, b, &c, n);

    printf("Result: ");
    for (int i = 0; i < n; i++)
        printf("%d ", c.x[i]);

    printf("\n");
    return 0;
}