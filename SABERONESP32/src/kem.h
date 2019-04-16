#ifndef KEM_H
#define KEM_H

#include "stdint.h"

#ifdef __cplusplus
extern "C" {
#endif

void indcpa_keypair(uint8_t *pk, uint8_t *sk);

void indcpa_client(uint8_t *pk, uint8_t *b_prime, uint8_t *c, uint8_t *key);

void indcpa_server(uint8_t *pk, uint8_t *b_prime, uint8_t *c, uint8_t *key);

void indcpa_kem_keypair(uint8_t *pk, uint8_t *sk);
void indcpa_kem_enc(uint8_t *message, uint8_t *noiseseed,const uint8_t *pk,  uint8_t *ciphertext);
void indcpa_kem_dec(const uint8_t *sk,const uint8_t *ciphertext, uint8_t message_dec[]);

int crypto_kem_keypair(unsigned char *pk, unsigned char *sk);
int crypto_kem_enc(unsigned char *c, unsigned char *k, const unsigned char *pk);
int crypto_kem_dec(unsigned char *k, const unsigned char *c, const unsigned char *sk);


uint64_t clock1,clock2;
uint64_t clock_kp_mv,clock_cl_mv, clock_kp_sm, clock_cl_sm;

#ifdef __cplusplus
}
#endif

#endif

