#ifndef CHASH_P_H
#define CHASH_P_H

#include <stddef.h>

void seed(void);
unsigned long xor128(void);
void random_shuffle(unsigned short *array, size_t size);
int mlt(int x, int y);
int mltn(int n, int x);
int mltnc(int n, int x);
int mltn2(int n, int x);
unsigned int mltn0(unsigned int n, unsigned int u);

#endif
