/*
 * poly_mul_ext.h
 *
 *  Created on: 2019年1月7日
 *      Author: Mai
 */

#ifndef POLYMUL_HW_H_
#define POLYMUL_HW_H_

#include "stdio.h"
#include "stdint.h"
#include "string.h"
#include "cpucycles.h"
#include "polymul_hw_common.h"
#include "polymul_hw_speed.h"
#include "poly.h"

#ifdef __cplusplus
extern "C" {
#endif

void schoolbook_mul(uint16_t* a, uint16_t* b, uint16_t* res, uint16_t q);


void poly_mul_ks1_64_reduce(uint16_t *a, uint16_t *b, uint16_t* res);

//void karatsuba_root(uint16_t* a, uint16_t* b, int64_t* res, int16_t n);

void karatsuba(uint16_t* a, uint16_t* b, uint16_t* res);

void reduce_karatsuba_sub_w(uint32_t* ks1_res32, uint16_t* res);

//void karatsuba_sub_w0w1w2(uint16_t* w0, uint16_t* w1, uint16_t* w2, int32_t* res);

void karatsuba_parallel(uint16_t* a, uint16_t* b, uint16_t* res);

//串行版本
void toom4_kyber(uint16_t* a, uint16_t* b, uint16_t* res);
//并行版本
void toom4_kyber_parallel(uint16_t* a, uint16_t* b, uint16_t* res);

void saber_karatsuba_parallel(uint16_t* a, uint16_t* b, uint16_t* res);
void saber_karatsuba_parallel_vector(uint16_t* a0, uint16_t* a1, uint16_t* a2,
		uint16_t* b0, uint16_t* b1, uint16_t* b2, uint16_t* acc);
void saber_karatsuba_parallel_vector_in_matrix(uint16_t* a0, uint16_t* a1, uint16_t* a2,
		uint16_t* b0, uint16_t* b1, uint16_t* b2, uint16_t* acc, uint16_t* root_res_full,
		uint16_t* acc_previous_step);


void saber_karatsuba_parallel_justintime(uint16_t* a, uint16_t* b, uint16_t* root_res_full,
		uint16_t* acc_previous, uint16_t acc_previous_need_rounding);

void reduce_root_res_full_to_acc(uint16_t* root_res_full, uint16_t* acc) ;
void rounding_acc(uint16_t* acc);

void reduce_root_res_full_to_acc_and_rounding(uint16_t* root_res_full, uint16_t* acc);

void MatrixVectorMul_parallel(polyvec *a, uint16_t skpv[SABER_K][SABER_N], uint16_t res[SABER_K][SABER_N], uint16_t mod,
		int16_t transpose);

//void MatrixVectorMul_rounding_parallel(polyvec *a, uint16_t skpv[SABER_K][SABER_N], uint16_t res[SABER_K][SABER_N], uint16_t mod,
//		int16_t transpose);

void toom4_general_ks2_mpi_init();
void toom4_general_ks2_mpi_free();
//通用版本，中间过程不约化,中间过程会使用mpi，使用前必须调用esp_hw_ks2_init()
void toom4_general(uint16_t* a1, uint16_t* b1, uint16_t* result);

void MatrixVectorMul_rounding_parallel(polyvec *a, polyvec *skpv, polyvec* res, uint16_t mod,
		int16_t transpose);
void dualcore_MatrixVectorMul_rounding_parallel(polyvec *a, polyvec *skpv, polyvec* res, uint16_t mod,
		int16_t transpose);

void toom4_saber_10bits_parallel(uint16_t* a, uint16_t* b, uint16_t* result) ;
void toom4_saber_10bits(uint16_t* a, uint16_t* b, uint16_t* result);

void toom4_saber_10bits_parallel_vector(uint16_t* va0, uint16_t* va1, uint16_t* va2,
		uint16_t* vb0, uint16_t* vb1, uint16_t* vb2, uint16_t* acc);

void toom4_saber_10bits_parallel_vector_speedtest(uint16_t* va0, uint16_t* va1, uint16_t* va2,
		uint16_t* vb0, uint16_t* vb1, uint16_t* vb2, uint16_t* acc);
#ifdef __cplusplus
}
#endif

#endif /* POLYMUL_HW_H_ */
