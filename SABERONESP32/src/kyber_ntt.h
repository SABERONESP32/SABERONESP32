/*
 * kyber_ntt.h
 *
 *  Created on: 2019年1月7日
 */

#ifndef KYBER_NTT_H_
#define KYBER_NTT_H_

#include <stdint.h>

#define KYBER_N 256
#define KYBER_Q 7681


#ifdef __cplusplus
extern "C" {
#endif



void ntt(uint16_t* poly);
void invntt(uint16_t* poly);

void poly_pointwise(uint16_t *r, const uint16_t *a, const uint16_t *b);

#ifdef __cplusplus
}
#endif

#endif /* KYBER_NTT_H_ */
