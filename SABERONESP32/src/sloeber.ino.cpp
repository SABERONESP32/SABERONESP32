#ifdef __IN_ECLIPSE__
//This is a automatic generated file
//Please do not modify this file
//If you touch this file your change will be overwritten during the next build
//This file has been generated on 2019-04-16 09:31:19

#include "Arduino.h"
#include "Arduino.h"
#include <WiFi.h>
#include "stdio.h"
#include "time.h"
#include "cpucycles.h"
#include "stdio.h"
#include "stdint.h"
#include "a_test.h"
#include "polymul_hw.h"
#include "polymul_hw_common.h"
#include "polymul_hw_speed.h"
#include "api.h"
#include "poly.h"
#include "rng.h"
#include "SABER_indcpa.h"
#include "cpucycles.h"
#include "verify.h"
#include "kem.h"
#include "dualcore.h"
#include "kem_opt.h"
#include "kyber_ntt.h"

void setup() ;
void print_runtime(uint32_t start, uint32_t finish, char* name) ;
void loop() ;
void alert_if_not_equal_uint16(uint16_t* a, uint16_t* b, uint32_t n, char* name) ;
void test_poly_multiplication() ;
int test_regs_read_write() ;
void test_ks1_ks2() ;
void sample_random(uint16_t* a, uint16_t n, uint16_t mod) ;
void check_cca_equal(uint8_t* ss_a, uint8_t* ss_b) ;
int test_kem_cca() ;


#include "Sabar_on_esp32.ino"

#endif
