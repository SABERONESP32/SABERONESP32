/*
 * poly_mul_ext.c
 *
 *  Created on: 2019年1月7日
 *      Author: Mai
 */

#include "kyber_ntt.h"
#include "polymul_hw.h"
#include "polymul_hw_speed.h"
#include "polymul_hw_common.h"

static uint32_t speedtest_start, speedtest_end;

void speedtest_print_cpucycles(char* name) {
	speedtest_end = cpucycles();
	int32_t t = (speedtest_end - speedtest_start);
	printf("%30s  %10d\n", name, t);
	speedtest_start = cpucycles();
}

void speedtest_print_title(char* name, int begin) {
	if (begin) {
		printf("\n");
	}
	printf("============== ");
	printf("%s", name);
	if(begin){
		printf(".begin =============\n");
	}else{
		printf(".end ===============\n\n");
	}
}

void speedtest_print_splitline() {
	speedtest_print_cpucycles("-------------------");
}

//重置起始时间
void speedtest_reset_startcpucycles() {
	speedtest_start = cpucycles();
}

void ks1_read_wait_speedtest(uint32_t* ks1_res32) {
	speedtest_reset_startcpucycles();
	while (REG_READ_NOCHECK(RSA_INTERRUPT_REG) != 1) {
		//printf(".\n");
	}
	speedtest_print_cpucycles("ks1_wait");
	//由于 ks1_res32是uint32，所以可以用内存拷贝 ,//注意res的长度是2*n-1
	//void *memcpy(void *dest, const void *src, size_t n);
	memcpy(ks1_res32, ((uint32_t *) REG_BASE_Z), 127 * 4);	//pbase_z->ks1_res32
	REG_WRITE_NOCHECK(RSA_INTERRUPT_REG, 1);
}

void ks1_speedtest(uint16_t* a, uint16_t* b) {
	REG_WRITE(RSA_MULT_MODE_REG, 15);
	speedtest_print_title("ks1_speedtest", 1);
	uint16_t* aLL = &a[0];	//64个
	uint16_t* bLL = &b[0];	//64个
	uint32_t ks1_res32[127];	//存放从寄存器读取的ks1结果

	ks1_write_uint16(aLL, bLL);
	//test wait
	ks1_read_wait_speedtest(ks1_res32);

	speedtest_start = cpucycles();
	ks1_write_uint16(aLL, bLL);
	speedtest_print_cpucycles("ks1_write_uint16");
	ks1_read(ks1_res32);
	speedtest_print_cpucycles("ks1_read_uint32");

	speedtest_start = cpucycles();
	ks1_write_uint16(aLL, bLL);
	ks1_read(ks1_res32);
	speedtest_print_cpucycles("ks1_full");
	speedtest_print_title("ks1_speedtest", 0);
}


void ntt_speedtest(uint16_t* a, uint16_t* b, uint16_t* res, int16_t n) {
	uint16_t aa[KYBER_N];
	uint16_t bb[KYBER_N];
	int i;
	for (i = 0; i < KYBER_N; i++) {
		aa[i] = a[i];
		bb[i] = b[i];
	}
	//speedtest_print_title("ntt", 1);
	speedtest_reset_startcpucycles();
	ntt(aa);
	speedtest_print_cpucycles("ntt_a");
	ntt(bb);
	speedtest_print_cpucycles("ntt_b");
	poly_pointwise(res, aa, bb);
	speedtest_print_cpucycles("poly_pointwise");
	invntt(res);
	speedtest_print_cpucycles("invntt");

	//speedtest_print_splitline();

	for (i = 0; i < KYBER_N; i++) {
		aa[i] = a[i];
		bb[i] = b[i];
	}

	speedtest_reset_startcpucycles();
	ntt(aa);
	ntt(bb);
	poly_pointwise(res, aa, bb);
	invntt(res);
	speedtest_print_cpucycles("ntt_mul");

	for(i=0; i < KYBER_N; i++){
		res[i] = res[i] % KYBER_Q;
	}

	//speedtest_print_title("ntt", 0);

}

