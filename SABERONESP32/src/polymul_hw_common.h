/*
 * poly_common.h
 *
 *  Created on: 2019年1月7日
 *      Author: Mai
 */

#ifndef POLYMUL_HW_COMMON_H_
#define POLYMUL_HW_COMMON_H_

#include "stdio.h"
#include "stdint.h"
#include "string.h"
#include "cpucycles.h"
#include "soc/hwcrypto_reg.h"
#include "esp_system.h"
#include "soc/dport_reg.h"

#define POLY_PRINT_MAX 8


#define SABER_N 256
#define KYBER_N 256

#define REDUCE_Q 8192// 8192//7681//8192
#define SABER_Q  8192
#define SABER_MOD_MASK 8191 //0x1FFFF


#define REG_BASE_X ((uint32_t) RSA_MEM_X_BLOCK_BASE)
#define REG_BASE_Y ((uint32_t) RSA_MEM_Z_BLOCK_BASE)
#define REG_BASE_Z ((uint32_t) RSA_MEM_Z_BLOCK_BASE)

//write value to register
#define REG_WRITE_NOCHECK(_r, _v) ({                                                                                           \
            (*(volatile uint32_t *)(_r)) = (_v);                                                                       \
        })

//read value from register
#define REG_READ_NOCHECK(_r) ({                                                                                                \
            (*(volatile uint32_t *)(_r));                                                                              \
        })

#ifdef __cplusplus
extern "C" {
#endif

//-----------------------
void esp32_acquire_hardware( void );
void esp32_release_hardware( void );

//--------------------
void ks1_full_int32_int32(int32_t* a, int32_t* b, int32_t* res);
void ks1_full_uint16_int32(uint16_t *a, uint16_t *b, int32_t* res);
void ks1_full_uint16_uint16(uint16_t *a, uint16_t *b, uint16_t* res);

void ks1_write_uint16(uint16_t* a, uint16_t* b);
void ks1_write_int32(int32_t* a, int32_t* b);
void ks1_read(uint32_t* ks1_res32);
void ks1_saber_uint16_uint16(uint16_t *a, uint16_t *b, uint16_t* res, uint16_t mask);


//--------------------

void poly_print_uint16(uint16_t *p, char* name, int n, int full);
void poly_print_int32(int32_t *p, char* name, int n, int full);
void poly_print_int64(int64_t *p, char* name, int n, int full);
void check_equal_int16(int16_t* a, int16_t* b, uint32_t n, char* name, int w_index);
void check_equal_uint16(uint16_t* a, uint16_t* b, uint32_t n, char* name);
void check_equal_int64(int64_t* a, int64_t* b, uint32_t n, char* name);


#ifdef __cplusplus
}
#endif

#endif /* POLYMUL_HW_COMMON_H_ */
