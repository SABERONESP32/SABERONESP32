/*
 * kem_unrolled.h
 *
 *  Created on: 2019年3月8日
 *      Author: Mai
 */

#ifndef KEM_OPT_H_
#define KEM_OPT_H_

#include <string.h>
#include <stdint.h>
#include <stdio.h>
#include "kem.h"
#include "kem.h"
#include "verify.h"
#include "rng.h"
#include "poly.h"
#include "pack_unpack.h"
#include "poly_mul.h"
#include "recon.h"
#include "rng.h"
#include "fips202.h"
#include "SABER_params.h"
#include "polymul_hw.h"
#include "cbd.h"
#include "dualcore.h"

#define KEM_OPT_TYPE_SPEED     0
#define KEM_OPT_TYPE_SIGLECORE 1
#define KEM_OPT_TYPE_DUALCORE  2

#define KEM_OPT_TYPE 2

#ifdef __cplusplus
extern "C" {
#endif

int crypto_kem_keypair_opt(unsigned char *pk, unsigned char *sk);
int crypto_kem_enc_opt(unsigned char *c, unsigned char *k, const unsigned char *pk);
int crypto_kem_dec_opt(unsigned char *k, const unsigned char *c, const unsigned char *sk);

int dualcore_crypto_kem_keypair_opt(unsigned char *pk, unsigned char *sk);
int dualcore_crypto_kem_enc_opt(unsigned char *c, unsigned char *k, const unsigned char *pk);
int dualcore_crypto_kem_dec_opt(unsigned char *k, const unsigned char *c, const unsigned char *sk);

int crypto_kem_keypair_opt_speedtest(unsigned char *pk, unsigned char *sk);
int crypto_kem_enc_opt_speedtest(unsigned char *c, unsigned char *k,
		const unsigned char *pk) ;
int crypto_kem_dec_opt_speedtest(unsigned char *k, const unsigned char *c,
		const unsigned char *sk);
#ifdef __cplusplus
}
#endif

#endif /* KEM_OPT_H_ */
