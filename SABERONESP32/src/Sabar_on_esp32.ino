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

//Arduino setup function
void setup() {
	Serial.begin(115200);
	WiFi.mode(WIFI_OFF);
	Serial.print("boot over!\n");
	//ESP.getCycleCount();
	dualcoreInitSemaphore();
}

#define SPEED_CPUCYCLES 1

#if SPEED_CPUCYCLES
#define  clock_wrap cpucycles
uint32_t start, end;
#else
#define  clock_wrap clock
clock_t start, end;
#endif

int NTEST = 100;

#if SPEED_CPUCYCLES
void print_runtime(uint32_t start, uint32_t finish, char* name) {
	int32_t t = (finish - start) / NTEST;
	printf("%30s  %10d\n", name, t);
}
#else
void print_runtime(clock_t start, clock_t finish, char* name) {
	double t = (finish - start) / (double) NTEST;
	printf("%29s  %.3f ms\n", name, t);
}
#endif

//Arduino loop function.
void loop() {
	NTEST = 10000;
	test_regs_read_write();

	NTEST = 100;
	test_ks1_ks2();
	test_poly_multiplication();
	test_kem_cca();
	printf("============== loop over ==============\n");
	delay(1000);
}

void alert_if_not_equal_uint16(uint16_t* a, uint16_t* b, uint32_t n, char* name) {
	int i;
	for (i = 0; i < n; i++) {
		if (a[i] != b[i]) {
			printf("%17s-> ", name);
			printf("NOT equal!!!!!\n");
			break;
		}
	}
}
//256 degree poly multiplication
void test_poly_multiplication() {

	speedtest_print_title((char*) "256_degree_poly_mul", 1);
	int i;
	const int n = SABER_N;
	int q;
	uint16_t a[n];
	uint16_t b[n];
	esp32_acquire_hardware();

	uint16_t nttres[n];
	uint16_t kara[n];
	uint16_t kara_parallel[n];
	uint16_t toom4[n];
	uint16_t toom4_parallel[n];
	uint16_t schoolbook[n];

	//============================
	// ntt(7681) based on Kyber code
	//============================
	q = 7681;
	sample_random(a, n, q);
	sample_random(b, n, q);
	ntt_speedtest(a, b, nttres, n);

	uint16_t at[n], bt[n], zt[n];
	memcpy(at, a, sizeof(a[0]) * n);
	memcpy(bt, b, sizeof(a[0]) * n);
	start = clock_wrap();
	for (i = 0; i < NTEST; i++) {
		ntt(at);
		ntt(bt);
		poly_pointwise(zt, at, bt);
		invntt(zt);
	}
	end = clock_wrap();
	print_runtime(start, end, (char*) "ntt(7681)");

	//============================
	// toom4(7681)
	//============================

	start = clock_wrap();
	for (i = 0; i < NTEST; i++) {
		toom4_kyber(a, b, toom4);
	}
	end = clock_wrap();
	print_runtime(start, end, (char*) "toom4(7681)");

	start = clock_wrap();
	for (i = 0; i < NTEST; i++) {
		toom4_kyber_parallel(a, b, toom4_parallel);
	}
	end = clock_wrap();
	print_runtime(start, end, (char*) "toom4_parallel(7681)");

	//	poly_print_uint16(ntt, "ntt", n, 1);
	//	poly_print_uint16(toom4, "toom4", n, 1);
	//	poly_print_uint16(toom4_parallel, "toom4_parallel", n, 1);

	alert_if_not_equal_uint16(nttres, toom4, n, (char*) "ntt vs. toom4");
	alert_if_not_equal_uint16(nttres, toom4_parallel, n, (char*) "ntt vs. toom4_parallel");

	//============================
	// karatsuba(8192)
	//============================
	q = 8192;
	sample_random(a, n, q);
	sample_random(b, n, q);

	schoolbook_mul(a, b, schoolbook, q);

	start = clock_wrap();
	for (i = 0; i < NTEST; i++) {
		karatsuba(a, b, kara);
	}
	end = clock_wrap();
	print_runtime(start, end, (char*) "karatsuba(8192)");

	start = clock_wrap();
	for (i = 0; i < NTEST; i++) {
		saber_karatsuba_parallel(a, b, kara_parallel);		//karatsuba_parallel(a, b, kara_parallel);
	}
	end = clock_wrap();
	print_runtime(start, end, (char*) "karatsuba_parallel(8192)");

	alert_if_not_equal_uint16(kara, schoolbook, n, (char*) "kara vs. schoolbook");
	alert_if_not_equal_uint16(kara_parallel, schoolbook, n, (char*) "kara_parallel vs. schoolbook");

	//============================
	// toom4(1024)
	//============================
	q = 1024;
	sample_random(a, n, q);
	sample_random(b, n, q);
	schoolbook_mul(a, b, schoolbook, q);

	start = clock_wrap();
	for (i = 0; i < NTEST; i++) {
		toom4_saber_10bits(a, b, toom4);
	}
	end = clock_wrap();
	print_runtime(start, end, (char*) "toom4(1024)");

	start = clock_wrap();
	for (i = 0; i < NTEST; i++) {
		toom4_saber_10bits_parallel(a, b, toom4_parallel);
	}
	end = clock_wrap();
	print_runtime(start, end, (char*) "toom4_parallel(1024)");

	alert_if_not_equal_uint16(schoolbook, toom4, n, (char*) "toom4 vs. schoolbook");
	alert_if_not_equal_uint16(schoolbook, toom4_parallel, n, (char*) "toom4_parallel vs. schoolbook");

	esp32_release_hardware();
	speedtest_print_title((char*) "256_degree_poly_mul", 0);
}

//read and write regs
int test_regs_read_write() {

	speedtest_print_title((char*) "regs_read_write", 1);
	esp32_acquire_hardware();
	const int write_n = 64;
	const int read_n = 64 * 2 - 1;
	uint16_t w16[write_n];
	uint32_t w32[write_n];
	uint16_t r16[read_n];
	uint32_t r32[read_n];
	int i, j;

	//write 64 regs

	start = clock_wrap();
	for (i = 0; i < NTEST; i++) {
		for (j = 0; j < write_n; j++) {
			REG_WRITE_NOCHECK(((uint32_t *) REG_BASE_X) + i, (uint32_t ) w16[j]);
		}
	}
	end = clock_wrap();
	print_runtime(start, end, (char*) "reg_write_loop(64)");

	start = clock_wrap();
	for (i = 0; i < NTEST; i++) {
		memcpy(((uint32_t *) REG_BASE_Y), w32, write_n * 4);
	}
	end = clock_wrap();
	print_runtime(start, end, (char*) "reg_write_memcpy(64)");

	REG_WRITE_NOCHECK(RSA_MULT_START_REG, 1);	//start
	while (REG_READ_NOCHECK(RSA_INTERRUPT_REG) != 1) {	//wait...
	}

	//read 127 regs

	start = clock_wrap();
	for (i = 0; i < NTEST; i++) {
		for (j = 0; j < read_n; j++) {
			r16[j] = REG_READ_NOCHECK(((uint32_t *) REG_BASE_Z) + j);
		}
	}
	end = clock_wrap();

	print_runtime(start, end, (char*) "reg_read_loop(127)");

	start = clock_wrap();
	for (i = 0; i < NTEST; i++) {
		memcpy(r32, ((uint32_t *) REG_BASE_Z), read_n * 4);
	}
	end = clock_wrap();
	print_runtime(start, end, (char*) "reg_read_memcpy(127)");

	//use the variables once, for avoiding gcc to delete the code block
	for (i = 0; i < read_n; i++) {
		if (r32[i] != r16[i]) {
			j++;
		}
	}
	for (i = 0; i < write_n; i++) {
		if (w32[i] != w16[i]) {
			j++;
		}
	}
	printf("test holder: %d\n", j);

	esp32_release_hardware();

	//printf("read write results: %d %d %d %d \n", r16[0], r32[0], w16[0], w32[0]);
	speedtest_print_title((char*) "regs_read_write", 0);
	return r16[0] + r32[0] + w16[0] + w32[0];
}
//KS1MUL and KS2MUL
void test_ks1_ks2() {
	uint16_t* a = kyber_sample_a;
	uint16_t* b = kyber_sample_b;
	esp32_acquire_hardware();
	toom4_general_ks2_mpi_init();
	ks1_speedtest(a, b);
	ks1_speedtest(a, b);	//two times
	ks2_speedtest();
	toom4_general_ks2_mpi_free();
	esp32_release_hardware();
}

void sample_random(uint16_t* a, uint16_t n, uint16_t mod) {
	int i;
	for (i = 0; i < n; i++) {
		a[i] = esp_random() % mod;
	}
}

void check_cca_equal(uint8_t* ss_a, uint8_t* ss_b) {
	int i;
	int equal = 1;
	for (i = 0; i < SABER_KEYBYTES; i++) {
		if (ss_a[i] != ss_b[i]) {
			equal = 0;
			printf(" !!!!! ERR CCA KEM !!!!!\n");
			break;
		}
	}
	if (equal) {
		printf(" ----- SABER CCA KEM OK ------\n");
	}
}

//Saber KEM (CCA)
int test_kem_cca() {

	uint8_t pk[CRYPTO_PUBLICKEYBYTES];
	uint8_t sk[CRYPTO_SECRETKEYBYTES];
	uint8_t ct[CRYPTO_CIPHERTEXTBYTES];
	uint8_t ss_a[CRYPTO_BYTES], ss_b[CRYPTO_BYTES];

//	unsigned char entropy_input[48];
//	for (i = 0; i < 48; i++){
//		entropy_input[i] = i;
//	}
	//randombytes_init(entropy_input, NULL, 256);

	esp32_acquire_hardware();
	dualcoreInitSemaphore();

	REG_WRITE_NOCHECK(RSA_MULT_MODE_REG, 15);

	//===========================
	// Single Core
	//===========================

	//Generation of secret key sk and public key pk pair
	int i;
	start = clock_wrap();
	for (i = 0; i < NTEST; i++) {
		crypto_kem_keypair_opt(pk, sk);
	}
	end = clock_wrap();
	print_runtime(start, end, (char*) "Gen(SingleCore)");

	//Key-Encapsulation call; input: pk; output: ciphertext c, shared-secret ss_a;
	start = clock_wrap();
	for (i = 0; i < NTEST; i++) {
		crypto_kem_enc_opt(ct, ss_a, pk);
	}
	end = clock_wrap();
	print_runtime(start, end, (char*) "Enc(SingleCore)");

	//Key-Decapsulation call; input: sk, c; output: shared-secret ss_b;
	start = clock_wrap();
	for (i = 0; i < NTEST; i++) {
		crypto_kem_dec_opt(ss_b, ct, sk);
	}
	end = clock_wrap();
	print_runtime(start, end, (char*) "Dec(SingleCore)");

	// Functional verification: check if ss_a == ss_b?
	check_cca_equal(ss_a, ss_b);

	//===========================
	// Dual Core
	//===========================

	start = clock_wrap();
	for (i = 0; i < NTEST; i++) {
		dualcore_crypto_kem_keypair_opt(pk, sk);
	}
	end = clock_wrap();
	print_runtime(start, end, (char*) "Gen(DualCore)");

	start = clock_wrap();
	for (i = 0; i < NTEST; i++) {
		dualcore_crypto_kem_enc_opt(ct, ss_a, pk);
	}
	end = clock_wrap();
	print_runtime(start, end, (char*) "Enc(DualCore)");

	start = clock_wrap();
	for (i = 0; i < NTEST; i++) {
		dualcore_crypto_kem_dec_opt(ss_b, ct, sk);
	}
	end = clock_wrap();
	print_runtime(start, end, (char*) "Dec(DualCore)");

	check_cca_equal(ss_a, ss_b);

	//===========================
	// Step speed test
	//===========================
	//	crypto_kem_keypair_opt_speedtest(pk, sk);
	//	crypto_kem_enc_opt_speedtest(ct, ss_a, pk);
	//	crypto_kem_dec_opt_speedtest(ss_b, ct, sk);
	//	check_cca_equal(ss_a, ss_b);

	//speedtest_print_title((char*)"Saber KEM", 0);

	esp32_release_hardware();

	return 0;
}

