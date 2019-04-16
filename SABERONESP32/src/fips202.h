#ifndef FIPS202_H
#define FIPS202_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define SHAKE128_RATE 168
#define SHAKE256_RATE 136
#define SHA3_256_RATE 136
#define SHA3_512_RATE 72

void keccak_squeezeblocks(unsigned char *h, unsigned long long int nblocks,
                                 uint64_t *s,
                                 unsigned int r);
void keccak_absorb(uint64_t *s,
                          unsigned int r,
                          const unsigned char *m, unsigned long long int mlen,
                          unsigned char p);

void shake128(unsigned char *output, unsigned long long outlen, const unsigned char *input,  unsigned long long inlen);

//void shake128(unsigned char *output, unsigned long long outlen, const unsigned char *input,  unsigned long long inlen);
void sha3_256(unsigned char *output, const unsigned char *input,  unsigned long long inlen);
void sha3_512(unsigned char *output, const unsigned char *input,  unsigned long long inlen);

#ifdef __cplusplus
}
#endif

#endif
