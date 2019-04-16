/*
 * polyext_mul_speed.h
 *
 *  Created on: 2019年1月7日
 *      Author: Mai
 */

#ifndef POLYMUL_HW_SPEED_H_
#define POLYMUL_HW_SPEED_H_

#include "stdio.h"
#include "stdint.h"
#include "string.h"
#include "cpucycles.h"
#include "soc/hwcrypto_reg.h"
#include "esp_system.h"
#include "soc/dport_reg.h"

#ifdef __cplusplus
extern "C" {
#endif

void speedtest_reset_startcpucycles();
void speedtest_print_cpucycles(char* name);
void speedtest_print_title(char* name, int begin) ;
void speedtest_print_splitline();

//-------------


void ks1_speedtest(uint16_t* a, uint16_t* b);
void ks2_speedtest();

void ntt_speedtest(uint16_t* a, uint16_t* b, uint16_t* res, int16_t n);

void toom4_kyber_parallel_speedtest(uint16_t* a, uint16_t* b, uint16_t* result);

void karatsuba_parallel_speedtest(uint16_t* a, uint16_t* b, uint16_t* res);




#ifdef __cplusplus
}
#endif


#endif /* POLYMUL_HW_SPEED_H_ */
