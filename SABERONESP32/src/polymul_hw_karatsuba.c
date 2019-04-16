/*
 * poly_mul_ext.c
 *
 *  Created on: 2019年1月7日
 *      Author: Mai
 */

#include "polymul_hw.h"
#include "bignum_lite.h"
#include "polymul_hw_speed.h"

//===========================================================

//约化硬件乘法器算出来的uint32，约化到int32即可，我们的w是以int32来计算的
static inline int32_t reduce_sub_w_uint32(uint32_t v, uint16_t q) {
#if REDUCE_Q == SABER_Q
	return v & SABER_MOD_MASK;
#else
	//return (uint16_t) (v % q);
//	uint32_t u;
//	u = v >> 13; //第一次必须m=1，不然v*m会超过uint32
//	u *= q;
//	v -= u;//能约化到max_v_reduce=124232553=2^26.9
//	return v;
	return v - (v >> 13) * q;
#endif
}

static inline uint16_t reduce_q_int32(int32_t v, uint16_t q) {
#if REDUCE_Q == SABER_Q
	return v & SABER_MOD_MASK;
#else
	int16_t t = v % REDUCE_Q;
	if ((t >> 15) == 0) {	//比判断t>0要快,与if((t & 0x8000) == 0)相同
		return t;
	}
	return t + REDUCE_Q;
#endif
}

static inline uint16_t reduce_q_add(uint16_t low, uint16_t high, uint16_t q) {
#if REDUCE_Q == SABER_Q
	return (low+high) & SABER_MOD_MASK;
#else
	uint16_t t = low + high;
	if (t > q) {
		return t - q;
	}
	return t;
#endif
}

//ks1结果约Q，长度为2 * 64 - 1=127
void karatsuba_reduce_sub_w(uint32_t* ks1_res32, int32_t* sub_w) {
	int i;
	for (i = 0; i < 127; i++) {
		sub_w[i] = reduce_sub_w_uint32(ks1_res32[i], REDUCE_Q);
	}
}

//(w0,w1,w2)->res, for kara-parallel. w_len=127, res_len=255, half_len=64
void karatsuba_sub_interpolation(int32_t* w0, uint32_t* w1_uint32, int32_t* w2, int32_t* res) {
	int w_len = 127;
	int half_len = 64;
	//memset(res, 0, 1020);	//res_len = 255, 4B * 255=1020B
	int i;
	memcpy(res, w0, w_len * 4);
	memcpy(res + w_len + 1, w2, w_len * 4);
	res[w_len] = 0;
	for (i = 0; i < w_len; i++) {
		res[i + half_len] += ((int32_t) reduce_sub_w_uint32(w1_uint32[i], REDUCE_Q)) - w0[i] - w2[i];
	}
}

//a，b长度为64（每个占32位，32*64=2048），res长度为64*2-1=127
//void ks1_full_uint16_uint16(uint16_t* a, uint16_t* b, uint16_t* res) {
//	//printf("-----poly_mul_ks1_64_hw2048_reduce-----\n");
//
//	REG_WRITE(RSA_MULT_MODE_REG, 15);
//	int i = 0;
//
//	//X低到高->64个X，64个0，
//	uint32_t *pbase_x = (uint32_t *) ((uint32_t) RSA_MEM_X_BLOCK_BASE);
//	uint32_t *pbase_y = (uint32_t *) ((uint32_t) RSA_MEM_Z_BLOCK_BASE);	//文档上是Z_BASE？？？
//	uint32_t *pbase_z = (uint32_t *) ((uint32_t) RSA_MEM_Z_BLOCK_BASE);
//
//	for (i = 0; i < 64; i++) {
//		REG_WRITE(pbase_x + i, (uint32_t )a[i]);
//	}
//	bzero(pbase_x + 64, 64 * 4);
//
//	//Y低到高->64个0，64个Y
//	bzero(pbase_y, 64 * 4);	// 64个0
//	for (i = 0; i < 64; i++) {
//		REG_WRITE(pbase_y + i + 64, (uint32_t )b[i]);
//	}
//	REG_WRITE(RSA_MULT_START_REG, 1);	//start
//	while (REG_READ(RSA_INTERRUPT_REG) != 1) {
//	}
//	for (i = 0; i < 127; i++) {	//注意res的长度是2*n-1
//		res[i] = reduce_q_uint32(REG_READ(pbase_z + i), REDUCE_Q);
//	}
//	REG_WRITE(RSA_INTERRUPT_REG, 1);	//clear interrupt
//}
//simple school book
void poly_mul_schoolbook_64(uint16_t *a, uint16_t *b, int32_t* res) {
	int n = 64;
	int i, j;
	int res_len = 2 * n - 1;
	uint32_t c[res_len];
//	for (i = 0; i < res_len; i++) {
//		c[i] = 0;
//	}
	memset(c, 0, sizeof(c[0]) * res_len);

	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			c[i + j] += (a[i] * b[j]);
		}
	}
	for (i = 0; i < res_len; i++) {
		res[i] = reduce_sub_w_uint32(c[i], REDUCE_Q);
	}
}

//第二次的karatsuba,n=128, half_len=64
void karatsuba_sub(uint16_t *a, uint16_t *b, int32_t* res, int16_t n) {
	int16_t half_len = n / 2;
	int16_t w_len = half_len * 2 - 1;	//127
	int16_t res_len = n * 2 - 1;		// 2* 128 - 1= 255
	uint16_t* aL = &a[0];
	uint16_t* bL = &b[0];
	uint16_t aLH[half_len];
	uint16_t bLH[half_len];
	uint16_t* aH = &a[half_len];
	uint16_t* bH = &b[half_len];
	int32_t w2[w_len], w1[w_len], w0[w_len];

	int i;
	for (i = 0; i < half_len; i++) {
		aLH[i] = reduce_q_add(aL[i], aH[i], REDUCE_Q);
		bLH[i] = reduce_q_add(bL[i], bH[i], REDUCE_Q);
	}
	//SW
	//poly_mul_ks1_64_reduce(aL, bL, w0);
	//poly_mul_ks1_64_reduce(aLH, bLH, w1);
	//poly_mul_ks1_64_reduce(aH, bH, w2);

	//HW
	ks1_full_uint16_int32(aL, bL, w0);
	ks1_full_uint16_int32(aLH, bLH, w1);
	ks1_full_uint16_int32(aH, bH, w2);

//	poly_print_uint16(w0, "sub_w0", w_len, 1);
//	poly_print_uint16(w1, "sub_w1", w_len, 1);
//	poly_print_uint16(w2, "sub_w2", w_len, 1);

	memset(res, 0, sizeof(res[0]) * res_len);

	memcpy(res, w0, w_len * 4);
	memcpy(res + w_len + 1, w2, w_len * 4);
	res[w_len] = 0; //余下一个255需要置零

	for (i = 0; i < w_len; i++) {
		res[half_len + i] += (int32_t) w1[i] - (int32_t) w0[i] - (int32_t) w2[i];
	}
}

//n=256
void karatsuba(uint16_t* a, uint16_t* b, uint16_t* res) {
	//printf("----karatsuba_root_reduce----\n");

	REG_WRITE_NOCHECK(RSA_MULT_MODE_REG, 15);
	int n = 256;

	int16_t half_len = 128;	//n / 2
	int16_t w_len = 255;	//half_len * 2 - 12^8
	int16_t res_full_len = 511;	//n * 2 - 1;	//res_full_len=511，2^9
	uint16_t* aL = &a[0];
	uint16_t* bL = &b[0];
	uint16_t aLH[half_len];
	uint16_t bLH[half_len];
	uint16_t* aH = &a[half_len];
	uint16_t* bH = &b[half_len];

	int i;
	for (i = 0; i < half_len; i++) {
		aLH[i] = reduce_q_add(aL[i], aH[i], REDUCE_Q);
		bLH[i] = reduce_q_add(bL[i], bH[i], REDUCE_Q);
	}

	int32_t w2[w_len], w1[w_len], w0[w_len];

	karatsuba_sub(aL, bL, w0, half_len);
	karatsuba_sub(aLH, bLH, w1, half_len);
	karatsuba_sub(aH, bH, w2, half_len);

	//poly_print_int32(w0, "w0", w_len, 1);
	//poly_print_int32(w1, "w1", w_len, 1);
	//poly_print_int32(w2, "w2", w_len, 1);

	int32_t res_full[res_full_len]; //res_full_len=511,最多21+9=30位

	memcpy(res_full, w0, w_len * 4);
	memcpy(res_full + w_len + 1, w2, w_len * 4);
	res_full[w_len] = 0; //余下一个255需要置零

	for (i = 0; i < w_len; i++) {
		res_full[half_len + i] += (w1[i] - w0[i] - w2[i]);
	}

//	poly_print_int32(w0, "w0", w_len, 1);
//	poly_print_int32(w1, "w1", w_len, 1);
//	poly_print_int32(w2, "w2", w_len, 1);
	//----------reduction-------
	for (i = 0; i < 255; i++) { //res[0 ~ (2n-2)-n=n-2=256-2=254]，剩下res[n-1]
		res[i] = reduce_q_int32(res_full[i] - res_full[i + n], REDUCE_Q);
	}
	res[n - 1] = reduce_q_int32(res_full[n - 1], REDUCE_Q);
}

//karatsuba_parallel， n=256
void karatsuba_parallel(uint16_t* a, uint16_t* b, uint16_t* res) {
	//printf("---karatsuba_parallel---\n");
	REG_WRITE_NOCHECK(RSA_MULT_MODE_REG, 15);
	int16_t n = 256;
	int16_t root_ab_len = 128;	// n / 2; //128
	int16_t root_w_len = 255;	// root_ab_len * 2 - 1; //root_w_len = 255
	int16_t root_res_full_len = 511;	// n * 2 - 1;	//res_full_len = 511，2^9

	int16_t sub_ab_len = 64;	// root_ab_len / 2;
	int16_t sub_w_len = 127;	//sub_ab_len * 2 - 1;

	int i;

	//-----------aL, bL---------------------------------
	uint16_t* aLL = &a[0];
	uint16_t* aLH = &a[64];	//sub_ab_len
	uint16_t* bLL = &b[0];
	uint16_t* bLH = &b[64];	//sub_ab_len

	//-----------aH, bH------------------------------
	uint16_t* aHL = &a[128];	//sub_ab_len*2
	uint16_t* aHH = &a[192];	//sub_ab_len*3
	uint16_t* bHL = &b[128];
	uint16_t* bHH = &b[192];

	//=== aL,bL->root_w0 ====

	ks1_write_uint16(aLL, bLL);
	int32_t root_w0[root_w_len], root_w1[root_w_len], root_w2[root_w_len];	//最多13+9=22位
	int32_t root_res_full[root_res_full_len]; //res_full_len=511,无需这里置零，使用w0和w2再加一个对其赋值即可
	uint16_t aS[root_ab_len];	//len=128, a的前一半加后一半
	uint16_t bS[root_ab_len];

	uint16_t* aSL = &aS[0];
	uint16_t* aSH = &aS[64];
	uint16_t* bSL = &bS[0];
	uint16_t* bSH = &bS[64];

	int32_t aTT[sub_ab_len];	//暂存下一次的write用的数据
	int32_t bTT[sub_ab_len];
	for (i = 0; i < root_ab_len; i++) {	//aS,bS
		aS[i] = reduce_q_add(a[i], a[i + root_ab_len], REDUCE_Q);
		bS[i] = reduce_q_add(b[i], b[i + root_ab_len], REDUCE_Q);
	}
	for (i = 0; i < sub_ab_len; i++) {
		aTT[i] = aLH[i];	//aLH = &a[64];
		bTT[i] = bLH[i];	//bLH = &b[64];
	}
	uint32_t ks1_res32[127];	//存放从寄存器读取的ks1结果
	ks1_read(ks1_res32);	// -> w0

	//-------

	ks1_write_int32(aTT, bTT);	//aLH, bLH
	for (i = 0; i < sub_ab_len; i++) {	//aLS, bLS
		aTT[i] = reduce_q_add(aLL[i], aLH[i], REDUCE_Q);	//aLS
		bTT[i] = reduce_q_add(bLL[i], bLH[i], REDUCE_Q);	//bLS
	}
	int32_t sub_w0[sub_w_len], sub_w2[sub_w_len];	//
	karatsuba_reduce_sub_w(ks1_res32, sub_w0);
	ks1_read(ks1_res32);	// -> w2

	//----------

	ks1_write_int32(aTT, bTT);	//aLS, bLS
	karatsuba_reduce_sub_w(ks1_res32, sub_w2);
	for (i = 0; i < sub_ab_len; i++) {
		aTT[i] = aHL[i];	//aHL
		bTT[i] = bHL[i];	//bHL
	}
	ks1_read(ks1_res32); // -> w1   (w0, ks1_res32, w2)->root_w0

	//==== aH,bH->root_w2 ====

	ks1_write_int32(aTT, bTT);	//aHL, bHL
	karatsuba_sub_interpolation(sub_w0, ks1_res32, sub_w2, root_w0);	// ->root_w0

	memcpy(root_res_full, root_w0, root_w_len * 4);	// copy root_w0->root_res_full
	root_res_full[root_w_len] = 0;

	for (i = 0; i < sub_ab_len; i++) {
		aTT[i] = aHH[i]; //aHH
		bTT[i] = bHH[i]; //bHH
	}
	ks1_read(ks1_res32);	// -> w0

	//-----

	ks1_write_int32(aTT, bTT);	//aHH, bHH
	karatsuba_reduce_sub_w(ks1_res32, sub_w0);
	for (i = 0; i < sub_ab_len; i++) {	//aHS, bHS
		aTT[i] = reduce_q_add(aHL[i], aHH[i], REDUCE_Q);	//aHS
		bTT[i] = reduce_q_add(bHL[i], bHH[i], REDUCE_Q);	//bHS
	}
	ks1_read(ks1_res32); // -> w2

	//-----

	ks1_write_int32(aTT, bTT); //aHS, bHS
	karatsuba_reduce_sub_w(ks1_res32, sub_w2);

	for (i = 0; i < sub_ab_len; i++) {
		aTT[i] = aSL[i];	//aSL
		bTT[i] = bSL[i];	//bSL
	}

	ks1_read(ks1_res32);	// -> w1   (w0, ks1_res32, w2)->root_w2

	//==== aS,bS->root_w1 ====

	ks1_write_int32(aTT, bTT); //aSL, bSL
	karatsuba_sub_interpolation(sub_w0, ks1_res32, sub_w2, root_w2); // -> root_w2

	memcpy(root_res_full + root_w_len + 1, root_w2, root_w_len * 4); // copy root_w2->root_res_full
	for (i = 0; i < sub_ab_len; i++) {
		aTT[i] = aSH[i];
		bTT[i] = bSH[i];
	}
	ks1_read(ks1_res32);	// -> w0

	//----------

	ks1_write_int32(aTT, bTT);	//aSH, bSH
	karatsuba_reduce_sub_w(ks1_res32, sub_w0);
	for (i = 0; i < root_w_len; i++) {	//(比较耗时，约3580)
		//先把root_w0和root_w2减掉，最后算出来root_w1之后再加上
		root_res_full[root_ab_len + i] += (-root_w0[i] - root_w2[i]);
	}
	for (i = 0; i < sub_ab_len; i++) {	//aSS, bSS
		aTT[i] = reduce_q_add(aSL[i], aSH[i], REDUCE_Q);	//aSS
		bTT[i] = reduce_q_add(bSL[i], bSH[i], REDUCE_Q);	//bSS
	}
	ks1_read(ks1_res32);	// -> w2

	//---------

	ks1_write_int32(aTT, bTT);	//aSS, bSS
	karatsuba_reduce_sub_w(ks1_res32, sub_w2);
	memcpy(root_w1, sub_w0, sub_w_len * 4);
	root_w1[sub_w_len] = 0;
	memcpy(root_w1 + sub_w_len + 1, sub_w2, sub_w_len * 4);
	for (i = 0; i < sub_w_len; i++) {
		//先把sub_w0和sub_w2减掉，最后再加sub_w1即可
		root_w1[i + sub_ab_len] += (-sub_w0[i] - sub_w2[i]);
	}
	ks1_read(ks1_res32);	// -> w1  (w0, ks1_res32, w2)->root_w1

	//=================以下不能异步了=====================
#ifdef DEBUG_T
	return ;
#endif
	for (i = 0; i < sub_w_len; i++) {	//并行中已减掉sub_w0和sub_w2
		root_w1[i + sub_ab_len] += reduce_sub_w_uint32(ks1_res32[i], REDUCE_Q);
	}

	for (i = 0; i < root_w_len; i++) {	//并行中已减掉root_w0和root_w2
		root_res_full[root_ab_len + i] += root_w1[i];
	}
	//----------reduction-------

	for (i = 0; i < 255; i++) { //res[0 ~ (2n-2)-n=n-2=256-2=254]，剩下res[n-1]
		res[i] = reduce_q_int32(root_res_full[i] - root_res_full[i + n], REDUCE_Q);
	}
	//n=256
	res[n - 1] = reduce_q_int32(root_res_full[n - 1], REDUCE_Q);
}

void karatsuba_parallel_speedtest(uint16_t* a, uint16_t* b, uint16_t* res) {
	speedtest_print_title("karatsuba_parallel_speedtest", 1);
	printf("REDUCE_Q=%d\n", REDUCE_Q);
	speedtest_reset_startcpucycles();
	speedtest_print_splitline();
	REG_WRITE_NOCHECK(RSA_MULT_MODE_REG, 15);
	speedtest_print_cpucycles("RSA_MULT_MODE_REG");
	int16_t n = 256;
	int16_t root_ab_len = 128;	// n / 2; //128
	int16_t root_w_len = 255;	// root_ab_len * 2 - 1; //root_w_len = 255
	int16_t root_res_full_len = 511;	// n * 2 - 1;	//res_full_len = 511，2^9

	int16_t sub_ab_len = 64;	// root_ab_len / 2;
	int16_t sub_w_len = 127;	//sub_ab_len * 2 - 1;

	int i;

	//-----------aL, bL---------------------------------
	uint16_t* aLL = &a[0];
	uint16_t* aLH = &a[64];	//sub_ab_len
	uint16_t* bLL = &b[0];
	uint16_t* bLH = &b[64];	//sub_ab_len

	//-----------aH, bH------------------------------
	uint16_t* aHL = &a[128];	//sub_ab_len*2
	uint16_t* aHH = &a[192];	//sub_ab_len*3
	uint16_t* bHL = &b[128];
	uint16_t* bHH = &b[192];

	speedtest_print_cpucycles("init_all_pointers");
	speedtest_print_splitline();

	//=== aL,bL->root_w0 ====

	ks1_write_uint16(aLL, bLL);
	speedtest_print_cpucycles("aLL_bLL_write");
	int32_t root_w0[root_w_len], root_w1[root_w_len], root_w2[root_w_len];	//最多13+9=22位
	int32_t root_res_full[root_res_full_len]; //res_full_len=511,无需这里置零，使用w0和w2再加一个对其赋值即可
	uint16_t aS[root_ab_len];	//len=128, a的前一半加后一半
	uint16_t bS[root_ab_len];

	uint16_t* aSL = &aS[0];
	uint16_t* aSH = &aS[64];
	uint16_t* bSL = &bS[0];
	uint16_t* bSH = &bS[64];

	int32_t aTT[sub_ab_len];	//暂存下一次的write用的数据
	int32_t bTT[sub_ab_len];
	for (i = 0; i < root_ab_len; i++) {	//aS,bS
		aS[i] = reduce_q_add(a[i], a[i + root_ab_len], REDUCE_Q);
		bS[i] = reduce_q_add(b[i], b[i + root_ab_len], REDUCE_Q);
	}
	for (i = 0; i < sub_ab_len; i++) {
		aTT[i] = aLH[i];	//aLH = &a[64];
		bTT[i] = bLH[i];	//bLH = &b[64];
	}
	uint32_t ks1_res32[127];	//存放从寄存器读取的ks1结果
	speedtest_print_cpucycles("aLL_bLL_wait_cpu_used");
	ks1_read(ks1_res32);	// -> w0
	speedtest_print_cpucycles("aLL_bLL_read");
	speedtest_print_splitline();

	//-------

	ks1_write_int32(aTT, bTT);	//aLH, bLH
	speedtest_print_cpucycles("aLH_bLH_write");
	for (i = 0; i < sub_ab_len; i++) {	//aLS, bLS
		aTT[i] = reduce_q_add(aLL[i], aLH[i], REDUCE_Q);	//aLS
		bTT[i] = reduce_q_add(bLL[i], bLH[i], REDUCE_Q);	//bLS
	}
	int32_t sub_w0[sub_w_len], sub_w2[sub_w_len];	//
	karatsuba_reduce_sub_w(ks1_res32, sub_w0);
	speedtest_print_cpucycles("aLH_bLH_wait_cpu_used");
	ks1_read(ks1_res32);	// -> w2
	speedtest_print_cpucycles("aLH_bLH_read");
	speedtest_print_splitline();

	//----------

	ks1_write_int32(aTT, bTT);	//aLS, bLS
	speedtest_print_cpucycles("aLS_bLS_write");
	karatsuba_reduce_sub_w(ks1_res32, sub_w2);
	for (i = 0; i < sub_ab_len; i++) {
		aTT[i] = aHL[i];	//aHL
		bTT[i] = bHL[i];	//bHL
	}
	speedtest_print_cpucycles("aLS_bLS_wait_cpu_used");
	ks1_read(ks1_res32); // -> w1   (w0, ks1_res32, w2)->root_w0
	speedtest_print_cpucycles("aLS_bLS_read");
	speedtest_print_splitline();
	speedtest_print_splitline();

	//==== aH,bH->root_w2 ====

	ks1_write_int32(aTT, bTT);	//aHL, bHL
	speedtest_print_cpucycles("aHL_bHL_write");
	karatsuba_sub_interpolation(sub_w0, ks1_res32, sub_w2, root_w0);	// ->root_w0

	memcpy(root_res_full, root_w0, root_w_len * 4);	// copy root_w0->root_res_full
	root_res_full[root_w_len] = 0;

	for (i = 0; i < sub_ab_len; i++) {
		aTT[i] = aHH[i]; //aHH
		bTT[i] = bHH[i]; //bHH
	}
	speedtest_print_cpucycles("aHL_bHL_wait_cpu_used");
	ks1_read(ks1_res32);	// -> w0
	speedtest_print_cpucycles("aHL_bHL_read");
	speedtest_print_splitline();

	//-----

	ks1_write_int32(aTT, bTT);	//aHH, bHH
	speedtest_print_cpucycles("aHH_bHH_write");
	karatsuba_reduce_sub_w(ks1_res32, sub_w0);
	for (i = 0; i < sub_ab_len; i++) {	//aHS, bHS
		aTT[i] = reduce_q_add(aHL[i], aHH[i], REDUCE_Q);	//aHS
		bTT[i] = reduce_q_add(bHL[i], bHH[i], REDUCE_Q);	//bHS
	}
	speedtest_print_cpucycles("aHH_bHH_wait_cpu_used");
	ks1_read(ks1_res32); // -> w2
	speedtest_print_cpucycles("aHH_bHH_read");
	speedtest_print_splitline();

	//-----

	ks1_write_int32(aTT, bTT); //aHS, bHS
	speedtest_print_cpucycles("aHS_bHS_write");
	karatsuba_reduce_sub_w(ks1_res32, sub_w2);

	for (i = 0; i < sub_ab_len; i++) {
		aTT[i] = aSL[i];	//aSL
		bTT[i] = bSL[i];	//bSL
	}
	speedtest_print_cpucycles("aHS_bHS_wait_cpu_used");
	ks1_read(ks1_res32);	// -> w1   (w0, ks1_res32, w2)->root_w2
	speedtest_print_cpucycles("aHS_bHS_read");
	speedtest_print_splitline();
	speedtest_print_splitline();

	//==== aS,bS->root_w1 ====

	ks1_write_int32(aTT, bTT); //aSL, bSL
	speedtest_print_cpucycles("aSL_bSL_write");
	karatsuba_sub_interpolation(sub_w0, ks1_res32, sub_w2, root_w2); // -> root_w2
	memcpy(root_res_full + root_w_len + 1, root_w2, root_w_len * 4); // copy root_w2->root_res_full
	for (i = 0; i < sub_ab_len; i++) {
		aTT[i] = aSH[i];
		bTT[i] = bSH[i];
	}
	speedtest_print_cpucycles("aSL_bSL_wait_cpu_used");
	ks1_read(ks1_res32);	// -> w0
	speedtest_print_cpucycles("aSL_bSL_read");
	speedtest_print_splitline();

	//----------

	ks1_write_int32(aTT, bTT);	//aSH, bSH
	speedtest_print_cpucycles("aSH_bSH_write");
	karatsuba_reduce_sub_w(ks1_res32, sub_w0);
	for (i = 0; i < root_w_len; i++) {	//(比较耗时，约3580)
		//先把root_w0和root_w2减掉，最后算出来root_w1之后再加上
		root_res_full[root_ab_len + i] += (-root_w0[i] - root_w2[i]);
	}
	for (i = 0; i < sub_ab_len; i++) {	//aSS, bSS
		aTT[i] = reduce_q_add(aSL[i], aSH[i], REDUCE_Q);	//aSS
		bTT[i] = reduce_q_add(bSL[i], bSH[i], REDUCE_Q);	//bSS
	}
	speedtest_print_cpucycles("aSH_bSH_wait_cpu_used");
	ks1_read(ks1_res32);	// -> w2
	speedtest_print_cpucycles("aSH_bSH_read");
	speedtest_print_splitline();

	//---------

	ks1_write_int32(aTT, bTT);	//aSS, bSS
	speedtest_print_cpucycles("aSS_bSS_write");
	karatsuba_reduce_sub_w(ks1_res32, sub_w2);
	memcpy(root_w1, sub_w0, sub_w_len * 4);
	root_w1[sub_w_len] = 0;
	memcpy(root_w1 + sub_w_len + 1, sub_w2, sub_w_len * 4);
	for (i = 0; i < sub_w_len; i++) {
		//先把sub_w0和sub_w2减掉，最后再加sub_w1即可
		root_w1[i + sub_ab_len] += (-sub_w0[i] - sub_w2[i]);
	}
	speedtest_print_cpucycles("aSS_bSS_wait_cpu_used");
	ks1_read(ks1_res32);	// -> w1  (w0, ks1_res32, w2)->root_w1
	speedtest_print_cpucycles("aSS_bSS_read");
	speedtest_print_splitline();

	//=================以下不能异步了=====================

//	poly_mul_schoolbook_64(aTT, bTT, sub_w2);
//	speedtest_print_cpucycles("poly_mul_schoolbook_64");

	for (i = 0; i < sub_w_len; i++) {	//并行中已减掉sub_w0和sub_w2
		root_w1[i + sub_ab_len] += ((int32_t) reduce_sub_w_uint32(ks1_res32[i], REDUCE_Q));
	}
	speedtest_print_cpucycles("root_w1+=sub_w1");

	for (i = 0; i < root_w_len; i++) {	//并行中已减掉root_w0和root_w2
		root_res_full[root_ab_len + i] += root_w1[i];
	}
	speedtest_print_cpucycles("root_res_full+=root_w1");
	//----------reduction-------

	for (i = 0; i < 255; i++) { //res[0 ~ (2n-2)-n=n-2=256-2=254]，剩下res[n-1]
		res[i] = reduce_q_int32(root_res_full[i] - root_res_full[i + n], REDUCE_Q);
	}
	//n=256
	res[n - 1] = reduce_q_int32(root_res_full[n - 1], REDUCE_Q);
	speedtest_print_cpucycles("reduction(X^N+1)");
	speedtest_print_title("karatsuba_parallel_speedtest", 0);
}

