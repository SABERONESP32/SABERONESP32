#include "kem_opt.h"
#include "polymul_hw.h"
#include "SABER_params.h"
#include "poly.h"

static inline uint16_t saber_reduce_uint16(uint16_t v) {
	return v & SABER_Q_MASK;
}

//ks1结果约Q，长度为2 * 64 - 1=127
void saber_karatsuba_reduce_sub_w(uint32_t* ks1_res32, uint16_t* sub_w) {
	int i;
	for (i = 0; i < 127; i++) {
		sub_w[i] = ks1_res32[i] & SABER_Q_MASK; //saber_reduce_uint16(ks1_res32[i]);
	}
}

//(w0,w1,w2)->res, for kara-parallel. w_len=127, res_len=255, half_len=64
void saber_karatsuba_sub_interpolation(uint16_t* w0, uint32_t* w1_uint32, uint16_t* w2, uint16_t* res) {
	int w_len = 127;
	int half_len = 64;
	//memset(res, 0, 1020);	//res_len = 255, 4B * 255=1020B
	int i;
	memcpy(res, w0, w_len * sizeof(w0[0])); //4);
	memcpy(res + w_len + 1, w2, w_len * sizeof(w2[0])); //4);
	res[w_len] = 0;
	for (i = 0; i < w_len; i++) {
		//res[i + half_len] += ((int32_t) saber_karatsuba_reduce_sub_w(w1_uint32[i], REDUCE_Q)) - w0[i] - w2[i];
		res[i + half_len] = ((res[i + half_len] + (w1_uint32[i] - w0[i] - w2[i])) & SABER_Q_MASK);
	}
}

//karatsuba_parallel， n=256
void saber_karatsuba_parallel(uint16_t* a, uint16_t* b, uint16_t* res) {
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
	uint16_t root_w0[root_w_len], root_w1[root_w_len], root_w2[root_w_len];	//最多13+9=22位
	uint16_t root_res_full[root_res_full_len]; //res_full_len=511,无需这里置零，使用w0和w2再加一个对其赋值即可
	uint16_t aS[root_ab_len];	//len=128, a的前一半加后一半
	uint16_t bS[root_ab_len];

	uint16_t* aSL = &aS[0];
	uint16_t* aSH = &aS[64];
	uint16_t* bSL = &bS[0];
	uint16_t* bSH = &bS[64];

	int32_t aTT[sub_ab_len];	//暂存下一次的write用的数据
	int32_t bTT[sub_ab_len];
	for (i = 0; i < root_ab_len; i++) {	//aS,bS
		aS[i] = saber_reduce_uint16(a[i] + a[i + root_ab_len]);
		bS[i] = saber_reduce_uint16(b[i] + b[i + root_ab_len]);
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
		aTT[i] = saber_reduce_uint16(aLL[i] + aLH[i]);	//aLS
		bTT[i] = saber_reduce_uint16(bLL[i] + bLH[i]);	//bLS
	}
	uint16_t sub_w0[sub_w_len], sub_w2[sub_w_len];	//
	saber_karatsuba_reduce_sub_w(ks1_res32, sub_w0);
	ks1_read(ks1_res32);	// -> w2

	//----------

	ks1_write_int32(aTT, bTT);	//aLS, bLS
	saber_karatsuba_reduce_sub_w(ks1_res32, sub_w2);
	for (i = 0; i < sub_ab_len; i++) {
		aTT[i] = aHL[i];	//aHL
		bTT[i] = bHL[i];	//bHL
	}
	ks1_read(ks1_res32); // -> w1   (w0, ks1_res32, w2)->root_w0

	//==== aH,bH->root_w2 ====

	ks1_write_int32(aTT, bTT);	//aHL, bHL
	saber_karatsuba_sub_interpolation(sub_w0, ks1_res32, sub_w2, root_w0);	// ->root_w0

	memcpy(root_res_full, root_w0, root_w_len * 4);	// copy root_w0->root_res_full
	root_res_full[root_w_len] = 0;

	for (i = 0; i < sub_ab_len; i++) {
		aTT[i] = aHH[i]; //aHH
		bTT[i] = bHH[i]; //bHH
	}
	ks1_read(ks1_res32);	// -> w0

	//-----

	ks1_write_int32(aTT, bTT);	//aHH, bHH
	saber_karatsuba_reduce_sub_w(ks1_res32, sub_w0);
	for (i = 0; i < sub_ab_len; i++) {	//aHS, bHS
		aTT[i] = saber_reduce_uint16(aHL[i] + aHH[i]);	//aHS
		bTT[i] = saber_reduce_uint16(bHL[i] + bHH[i]);	//bHS
	}
	ks1_read(ks1_res32); // -> w2

	//-----

	ks1_write_int32(aTT, bTT); //aHS, bHS
	saber_karatsuba_reduce_sub_w(ks1_res32, sub_w2);

	for (i = 0; i < sub_ab_len; i++) {
		aTT[i] = aSL[i];	//aSL
		bTT[i] = bSL[i];	//bSL
	}

	ks1_read(ks1_res32);	// -> w1   (w0, ks1_res32, w2)->root_w2

	//==== aS,bS->root_w1 ====

	ks1_write_int32(aTT, bTT); //aSL, bSL
	saber_karatsuba_sub_interpolation(sub_w0, ks1_res32, sub_w2, root_w2); // -> root_w2

	memcpy(root_res_full + root_w_len + 1, root_w2, root_w_len * sizeof(root_w2[0])); // copy root_w2->root_res_full
	for (i = 0; i < sub_ab_len; i++) {
		aTT[i] = aSH[i];
		bTT[i] = bSH[i];
	}
	ks1_read(ks1_res32);	// -> w0

	//----------

	ks1_write_int32(aTT, bTT);	//aSH, bSH
	saber_karatsuba_reduce_sub_w(ks1_res32, sub_w0);
	for (i = 0; i < root_w_len; i++) {	//(比较耗时，约3580)
		//先把root_w0和root_w2减掉，最后算出来root_w1之后再加上
		root_res_full[root_ab_len + i] += (-root_w0[i] - root_w2[i]);
	}
	for (i = 0; i < sub_ab_len; i++) {	//aSS, bSS
		aTT[i] = saber_reduce_uint16(aSL[i] + aSH[i]);	//aSS
		bTT[i] = saber_reduce_uint16(bSL[i] + bSH[i]);	//bSS
	}
	ks1_read(ks1_res32);	// -> w2

	//---------

	ks1_write_int32(aTT, bTT);	//aSS, bSS
	saber_karatsuba_reduce_sub_w(ks1_res32, sub_w2);
	memcpy(root_w1, sub_w0, sub_w_len * sizeof(sub_w0[0]));
	root_w1[sub_w_len] = 0;
	memcpy(root_w1 + sub_w_len + 1, sub_w2, sub_w_len * sizeof(sub_w2[0]));
	for (i = 0; i < sub_w_len; i++) {
		//先把sub_w0和sub_w2减掉，最后再加sub_w1即可
		root_w1[i + sub_ab_len] += (-sub_w0[i] - sub_w2[i]);
	}
	ks1_read(ks1_res32);	// -> w1  (w0, ks1_res32, w2)->root_w1

	//=================以下不能异步了=====================
	//要传给下一次的数据包括root_w1, ks1_res32, root_res_full, res (和整个的acc)

	for (i = 0; i < sub_w_len; i++) {	//并行中已减掉sub_w0和sub_w2
		root_w1[i + sub_ab_len] += (ks1_res32[i]);	// & SABER_Q_MASK);//reduce_sub_w_uint32(ks1_res32[i], REDUCE_Q);
	}

	for (i = 0; i < root_w_len; i++) {	//并行中已减掉root_w0和root_w2
		root_res_full[root_ab_len + i] += root_w1[i];
	}
	//----------reduction-------

	for (i = 0; i < 255; i++) { //res[0 ~ (2n-2)-n=n-2=256-2=254]，剩下res[n-1]
		//这里可能有问题
		res[i] = saber_reduce_uint16(root_res_full[i] - root_res_full[i + n]);
	}
	//n=256
	res[n - 1] = saber_reduce_uint16(root_res_full[n - 1]);
}

//====================================================================================

void reduce_root_res_full_to_acc_and_rounding(uint16_t* root_res_full, uint16_t* acc) {
	int i;
	const int n = SABER_N;
//	for (i = 0; i < 255; i++) { //res[0 ~ (2n-2)-n=n-2=256-2=254]，剩下res[n-1]
//		//res[i] = saber_reduce_uint16(root_res_full[i] - root_res_full[i + n]);
//		root_res_full[i] -= root_res_full[i + n];
//	}
//	for (i = 0; i < n; i++) {
//		acc[i] = ((acc[i] + root_res_full[i]) & SABER_Q_MASK);
//	}
	//res[i][j] = (res[i][j] + 4) >> 3;//rounding
	//两个循环合并一下
	for (i = 0; i < 255; i++) { //res[0 ~ (2n-2)-n=n-2=256-2=254]，剩下res[n-1]
		acc[i] = ((acc[i] + root_res_full[i] - root_res_full[i + n]) & SABER_Q_MASK);
		acc[i] = acc[i] + 4;
		acc[i] = acc[i] >> 3;
	}
	acc[n - 1] = ((acc[n - 1] + root_res_full[n - 1]) & SABER_Q_MASK);
	acc[n - 1] = acc[n - 1] + 4;
	acc[n - 1] = acc[n - 1] >> 3;
}

void reduce_root_res_full_to_acc_and_rounding_step0(uint16_t* root_res_full, uint16_t* acc) {
	int i;
	const int n = SABER_N;
	for (i = 0; i < 255; i++) { //res[0 ~ (2n-2)-n=n-2=256-2=254]，剩下res[n-1]
		acc[i] = ((acc[i] + root_res_full[i] - root_res_full[i + n]) & SABER_Q_MASK);
	}
	acc[n - 1] = ((acc[n - 1] + root_res_full[n - 1]) & SABER_Q_MASK);
}
void reduce_root_res_full_to_acc_and_rounding_step1(uint16_t* acc) {
	int i;
	const int n = SABER_N;
	for (i = 0; i < n; i++) { //res[0 ~ (2n-2)-n=n-2=256-2=254]，剩下res[n-1]
		acc[i] = acc[i] + 4;
		acc[i] = acc[i] >> 3;
	}
}

//acc_previous_step是上一次的acc指针，如果是NULL就是第一次的乘法，没有上一步的acc
void saber_karatsuba_parallel_vector_in_matrix_and_rounding(uint16_t* a0, uint16_t* a1, uint16_t* a2,
		uint16_t* b0, uint16_t* b1, uint16_t* b2, uint16_t* acc,
		uint16_t* root_res_full, uint16_t* acc_previous_step) {
//printf("---karatsuba_parallel---\n");
	REG_WRITE_NOCHECK(RSA_MULT_MODE_REG, 15);
	const int16_t n = 256;
	int16_t root_ab_len = 128;	// n / 2; //128
	int16_t root_w_len = 255;	// root_ab_len * 2 - 1; //root_w_len = 255
//int16_t root_res_full_len = 511;	// n * 2 - 1;	//res_full_len = 511，2^9

	int16_t sub_ab_len = 64;	// root_ab_len / 2;
	int16_t sub_w_len = 127;	//sub_ab_len * 2 - 1;

	int i;

//============ a0, b0==================

//-----------aL, bL---------------------------------
	uint16_t* aL = &a0[0];
	uint16_t* aH = &a0[128];  //root_ab_len
	uint16_t* bL = &b0[0];
	uint16_t* bH = &b0[128];

	uint16_t* aLL = &a0[0];
	uint16_t* aLH = &a0[64];	//sub_ab_len
	uint16_t* bLL = &b0[0];
	uint16_t* bLH = &b0[64];	//sub_ab_len

//-----------aH, bH------------------------------
	uint16_t* aHL = &a0[128];	//sub_ab_len*2
	uint16_t* aHH = &a0[192];	//sub_ab_len*3
	uint16_t* bHL = &b0[128];
	uint16_t* bHH = &b0[192];

//=== aL,bL->root_w0 ====

	ks1_write_uint16(aLL, bLL);
	uint16_t root_w0[root_w_len], root_w1[root_w_len], root_w2[root_w_len];	//最多13+9=22位
//uint16_t root_res_full[root_res_full_len]; //res_full_len=511,无需这里置零，使用w0和w2再加一个对其赋值即可
	uint16_t aS[root_ab_len];	//len=128, a的前一半加后一半
	uint16_t bS[root_ab_len];

	uint16_t* aSL = &aS[0];
	uint16_t* aSH = &aS[64];
	uint16_t* bSL = &bS[0];
	uint16_t* bSH = &bS[64];

	int32_t aTT[sub_ab_len];	//暂存下一次的write用的数据
	int32_t bTT[sub_ab_len];
	for (i = 0; i < root_ab_len; i++) {	//aS,bS
		aS[i] = saber_reduce_uint16(aL[i] + aH[i]);
		bS[i] = saber_reduce_uint16(bL[i] + bH[i]);
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
		aTT[i] = saber_reduce_uint16(aLL[i] + aLH[i]);	//aLS
		bTT[i] = saber_reduce_uint16(bLL[i] + bLH[i]);	//bLS
	}
	uint16_t sub_w0[sub_w_len], sub_w2[sub_w_len];	//
	saber_karatsuba_reduce_sub_w(ks1_res32, sub_w0);
	ks1_read(ks1_res32);	// -> w2

//----------

	ks1_write_int32(aTT, bTT);	//aLS, bLS
	saber_karatsuba_reduce_sub_w(ks1_res32, sub_w2);
	for (i = 0; i < sub_ab_len; i++) {
		aTT[i] = aHL[i];	//aHL
		bTT[i] = bHL[i];	//bHL
	}
	ks1_read(ks1_res32); // -> w1   (w0, ks1_res32, w2)->root_w0

//==== aH,bH->root_w2 ====

	ks1_write_int32(aTT, bTT);	//aHL, bHL
	//speedtest_reset_startcpucycles();

	saber_karatsuba_sub_interpolation(sub_w0, ks1_res32, sub_w2, root_w0);	// ->root_w0
	if (acc_previous_step != NULL) {
		reduce_root_res_full_to_acc_and_rounding_step0(root_res_full, acc_previous_step);
	}
	for (i = 0; i < sub_ab_len; i++) {
		aTT[i] = aHH[i]; //aHH
		bTT[i] = bHH[i]; //bHH
	}
	//speedtest_print_cpucycles("aHL, bHL cpu uesd");
	//7953,8922,10204消耗过头了，把copy root_w0->root_res_full放到下一步
	//7273,9476消耗还是多，把reduce_root_res_full_to_acc_and_rounding拆分为acc和rounding两步
	//已OK，这里是3415, 6754, 8472, 8206
	ks1_read(ks1_res32);	// -> w0

//-----

	ks1_write_int32(aTT, bTT);	//aHH, bHH
	//speedtest_reset_startcpucycles();
	saber_karatsuba_reduce_sub_w(ks1_res32, sub_w0);
	if (acc_previous_step != NULL) {
		reduce_root_res_full_to_acc_and_rounding_step1(acc_previous_step);
	}

	memcpy(root_res_full, root_w0, root_w_len * 4);	// copy root_w0->root_res_full
	root_res_full[root_w_len] = 0;
	for (i = 0; i < sub_ab_len; i++) {	//aHS, bHS
		aTT[i] = saber_reduce_uint16(aHL[i] + aHH[i]);	//aHS
		bTT[i] = saber_reduce_uint16(bHL[i] + bHH[i]);	//bHS
	}
	//speedtest_print_cpucycles("aHH, bHH cpu uesd");//2846, 5131
	ks1_read(ks1_res32); // -> w2

//-----

	ks1_write_int32(aTT, bTT); //aHS, bHS
	saber_karatsuba_reduce_sub_w(ks1_res32, sub_w2);

	for (i = 0; i < sub_ab_len; i++) {
		aTT[i] = aSL[i];	//aSL
		bTT[i] = bSL[i];	//bSL
	}

	ks1_read(ks1_res32);	// -> w1   (w0, ks1_res32, w2)->root_w2

//==== aS,bS->root_w1 ====

	ks1_write_int32(aTT, bTT); //aSL, bSL
	saber_karatsuba_sub_interpolation(sub_w0, ks1_res32, sub_w2, root_w2); // -> root_w2

	memcpy(root_res_full + root_w_len + 1, root_w2, root_w_len * sizeof(root_w2[0])); // copy root_w2->root_res_full
	for (i = 0; i < sub_ab_len; i++) {
		aTT[i] = aSH[i];
		bTT[i] = bSH[i];
	}
	ks1_read(ks1_res32);	// -> w0

//----------

	ks1_write_int32(aTT, bTT);	//aSH, bSH
	saber_karatsuba_reduce_sub_w(ks1_res32, sub_w0);
	for (i = 0; i < root_w_len; i++) {	//(比较耗时，约3580)
		//先把root_w0和root_w2减掉，最后算出来root_w1之后再加上
		root_res_full[root_ab_len + i] += (-root_w0[i] - root_w2[i]);
	}
	for (i = 0; i < sub_ab_len; i++) {	//aSS, bSS
		aTT[i] = saber_reduce_uint16(aSL[i] + aSH[i]);	//aSS
		bTT[i] = saber_reduce_uint16(bSL[i] + bSH[i]);	//bSS
	}
	ks1_read(ks1_res32);	// -> w2

//---------

	ks1_write_int32(aTT, bTT);	//aSS, bSS
	saber_karatsuba_reduce_sub_w(ks1_res32, sub_w2);
	memcpy(root_w1, sub_w0, sub_w_len * sizeof(sub_w0[0]));
	root_w1[sub_w_len] = 0;
	memcpy(root_w1 + sub_w_len + 1, sub_w2, sub_w_len * sizeof(sub_w2[0]));
	for (i = 0; i < sub_w_len; i++) {
		//先把sub_w0和sub_w2减掉，最后再加sub_w1即可
		root_w1[i + sub_ab_len] += (-sub_w0[i] - sub_w2[i]);
	}

// a1, b1 pointer
	aL = &a1[0];
	aH = &a1[128];  //root_ab_len
	bL = &b1[0];
	bH = &b1[128];

//-----------aL, bL---------------------------------

	aLL = &a1[0];
	aLH = &a1[64];	//sub_ab_len
	bLL = &b1[0];
	bLH = &b1[64];	//sub_ab_len

//-----------aH, bH------------------------------
	aHL = &a1[128];	//sub_ab_len*2
	aHH = &a1[192];	//sub_ab_len*3
	bHL = &b1[128];
	bHH = &b1[192];

	for (i = 0; i < sub_ab_len; i++) {
		aTT[i] = aLL[i];
		bTT[i] = bLL[i];
	}

	ks1_read(ks1_res32);	// -> w1  (w0, ks1_res32, w2)->root_w1

//#########################################################

//============ a1, b1 ==============

//=== aL,bL->root_w0 ====

//ks1_write_uint16(aLL, bLL);
	ks1_write_int32(aTT, bTT);		//aLL,bLL已经处理过了

	for (i = 0; i < root_ab_len; i++) {	//aS,bS
		aS[i] = saber_reduce_uint16(aL[i] + aH[i]);
		bS[i] = saber_reduce_uint16(bL[i] + bH[i]);
	}
	for (i = 0; i < sub_ab_len; i++) {
		aTT[i] = aLH[i];	//aLH = &a[64];
		bTT[i] = bLH[i];	//bLH = &b[64];
	}

//==== a0, b0 rest step.0
	for (i = 0; i < sub_w_len; i++) {	//并行中已减掉sub_w0和sub_w2
		root_w1[i + sub_ab_len] += (ks1_res32[i]);	// & SABER_Q_MASK);//reduce_sub_w_uint32(ks1_res32[i], REDUCE_Q);
	}
//======================

	ks1_read(ks1_res32);	// -> w0

//-------

	ks1_write_int32(aTT, bTT);	//aLH, bLH
	for (i = 0; i < sub_ab_len; i++) {	//aLS, bLS
		aTT[i] = saber_reduce_uint16(aLL[i] + aLH[i]);	//aLS
		bTT[i] = saber_reduce_uint16(bLL[i] + bLH[i]);	//bLS
	}
	saber_karatsuba_reduce_sub_w(ks1_res32, sub_w0);

//==== a0, b0 rest step.1
	for (i = 0; i < root_w_len; i++) {	//并行中已减掉root_w0和root_w2
		root_res_full[root_ab_len + i] += root_w1[i];
	}
//=======================

	ks1_read(ks1_res32);	// -> w2

//----------

	ks1_write_int32(aTT, bTT);	//aLS, bLS
	saber_karatsuba_reduce_sub_w(ks1_res32, sub_w2);
	for (i = 0; i < sub_ab_len; i++) {
		aTT[i] = aHL[i];	//aHL
		bTT[i] = bHL[i];	//bHL
	}

//==== a0, b0 rest step.2
	for (i = 0; i < 255; i++) { //reduction,res[0 ~ (2n-2)-n=n-2=256-2=254]，剩下res[n-1]
		//res[i] = saber_reduce_uint16(root_res_full[i] - root_res_full[i + n]);
		root_res_full[i] -= root_res_full[i + n];
	}
//n=256
//res[n - 1] = saber_reduce_uint16(root_res_full[n - 1]);

	memcpy(acc, root_res_full, sizeof(root_res_full[0]) * n);

//======================

	ks1_read(ks1_res32); // -> w1   (w0, ks1_res32, w2)->root_w0

//==== aH,bH->root_w2 ====

	ks1_write_int32(aTT, bTT);	//aHL, bHL
	saber_karatsuba_sub_interpolation(sub_w0, ks1_res32, sub_w2, root_w0);	// ->root_w0

	memcpy(root_res_full, root_w0, root_w_len * 4);	// copy root_w0->root_res_full
	root_res_full[root_w_len] = 0;

	for (i = 0; i < sub_ab_len; i++) {
		aTT[i] = aHH[i]; //aHH
		bTT[i] = bHH[i]; //bHH
	}
	ks1_read(ks1_res32);	// -> w0

//-----

	ks1_write_int32(aTT, bTT);	//aHH, bHH
	saber_karatsuba_reduce_sub_w(ks1_res32, sub_w0);
	for (i = 0; i < sub_ab_len; i++) {	//aHS, bHS
		aTT[i] = saber_reduce_uint16(aHL[i] + aHH[i]);	//aHS
		bTT[i] = saber_reduce_uint16(bHL[i] + bHH[i]);	//bHS
	}
	ks1_read(ks1_res32); // -> w2

//-----

	ks1_write_int32(aTT, bTT); //aHS, bHS
	saber_karatsuba_reduce_sub_w(ks1_res32, sub_w2);

	for (i = 0; i < sub_ab_len; i++) {
		aTT[i] = aSL[i];	//aSL
		bTT[i] = bSL[i];	//bSL
	}

	ks1_read(ks1_res32);	// -> w1   (w0, ks1_res32, w2)->root_w2

//==== aS,bS->root_w1 ====

	ks1_write_int32(aTT, bTT); //aSL, bSL
	saber_karatsuba_sub_interpolation(sub_w0, ks1_res32, sub_w2, root_w2); // -> root_w2

	memcpy(root_res_full + root_w_len + 1, root_w2, root_w_len * sizeof(root_w2[0])); // copy root_w2->root_res_full
	for (i = 0; i < sub_ab_len; i++) {
		aTT[i] = aSH[i];
		bTT[i] = bSH[i];
	}
	ks1_read(ks1_res32);	// -> w0

//----------

	ks1_write_int32(aTT, bTT);	//aSH, bSH
	saber_karatsuba_reduce_sub_w(ks1_res32, sub_w0);
	for (i = 0; i < root_w_len; i++) {	//(比较耗时，约3580)
		//先把root_w0和root_w2减掉，最后算出来root_w1之后再加上
		root_res_full[root_ab_len + i] += (-root_w0[i] - root_w2[i]);
	}
	for (i = 0; i < sub_ab_len; i++) {	//aSS, bSS
		aTT[i] = saber_reduce_uint16(aSL[i] + aSH[i]);	//aSS
		bTT[i] = saber_reduce_uint16(bSL[i] + bSH[i]);	//bSS
	}
	ks1_read(ks1_res32);	// -> w2

//---------

	ks1_write_int32(aTT, bTT);	//aSS, bSS
	saber_karatsuba_reduce_sub_w(ks1_res32, sub_w2);
	memcpy(root_w1, sub_w0, sub_w_len * sizeof(sub_w0[0]));
	root_w1[sub_w_len] = 0;
	memcpy(root_w1 + sub_w_len + 1, sub_w2, sub_w_len * sizeof(sub_w2[0]));
	for (i = 0; i < sub_w_len; i++) {
		//先把sub_w0和sub_w2减掉，最后再加sub_w1即可
		root_w1[i + sub_ab_len] += (-sub_w0[i] - sub_w2[i]);
	}

// a2, b2 pointer
	aL = &a2[0];
	aH = &a2[128];  //root_ab_len
	bL = &b2[0];
	bH = &b2[128];

//-----------aL, bL---------------------------------

	aLL = &a2[0];
	aLH = &a2[64];	//sub_ab_len
	bLL = &b2[0];
	bLH = &b2[64];	//sub_ab_len

//-----------aH, bH------------------------------
	aHL = &a2[128];	//sub_ab_len*2
	aHH = &a2[192];	//sub_ab_len*3
	bHL = &b2[128];
	bHH = &b2[192];

	for (i = 0; i < sub_ab_len; i++) {
		aTT[i] = aLL[i];
		bTT[i] = bLL[i];
	}

	ks1_read(ks1_res32);	// -> w1  (w0, ks1_res32, w2)->root_w1

//#########################################################

//============ a2, b2 ==============

//=== aL,bL->root_w0 ====

//ks1_write_uint16(aLL, bLL);
	ks1_write_int32(aTT, bTT);		//aLL,bLL已经处理过了

	for (i = 0; i < root_ab_len; i++) {	//aS,bS
		aS[i] = saber_reduce_uint16(aL[i] + aH[i]);
		bS[i] = saber_reduce_uint16(bL[i] + bH[i]);
	}
	for (i = 0; i < sub_ab_len; i++) {
		aTT[i] = aLH[i];	//aLH = &a[64];
		bTT[i] = bLH[i];	//bLH = &b[64];
	}

//==== a1, b1 rest step.0
	for (i = 0; i < sub_w_len; i++) {	//并行中已减掉sub_w0和sub_w2
		root_w1[i + sub_ab_len] += (ks1_res32[i]);	// & SABER_Q_MASK);//reduce_sub_w_uint32(ks1_res32[i], REDUCE_Q);
	}
//======================

	ks1_read(ks1_res32);	// -> w0

//-------

	ks1_write_int32(aTT, bTT);	//aLH, bLH
	for (i = 0; i < sub_ab_len; i++) {	//aLS, bLS
		aTT[i] = saber_reduce_uint16(aLL[i] + aLH[i]);	//aLS
		bTT[i] = saber_reduce_uint16(bLL[i] + bLH[i]);	//bLS
	}
	saber_karatsuba_reduce_sub_w(ks1_res32, sub_w0);

//==== a1, b1 rest step.1
	for (i = 0; i < root_w_len; i++) {	//并行中已减掉root_w0和root_w2
		root_res_full[root_ab_len + i] += root_w1[i];
	}
//=======================

	ks1_read(ks1_res32);	// -> w2

//----------

	ks1_write_int32(aTT, bTT);	//aLS, bLS
	saber_karatsuba_reduce_sub_w(ks1_res32, sub_w2);
	for (i = 0; i < sub_ab_len; i++) {
		aTT[i] = aHL[i];	//aHL
		bTT[i] = bHL[i];	//bHL
	}

//==== a1, b1 rest step.2
	for (i = 0; i < 255; i++) { //reduction,res[0 ~ (2n-2)-n=n-2=256-2=254]，剩下res[n-1]
		//这里可能有问题
		//res[i] = saber_reduce_uint16(root_res_full[i] - root_res_full[i + n]);
		//root_res_full[i] = saber_reduce_uint16(root_res_full[i] - root_res_full[i + n]);
		root_res_full[i] -= root_res_full[i + n];
	}
//n=256
//res[n - 1] = saber_reduce_uint16(root_res_full[n - 1]);		// a1, b1, over

	for (i = 0; i < n; i++) {
		//acc[i] += res[i];
		acc[i] += root_res_full[i];
	}

//======================

	ks1_read(ks1_res32); // -> w1   (w0, ks1_res32, w2)->root_w0

//==== aH,bH->root_w2 ====

	ks1_write_int32(aTT, bTT);	//aHL, bHL
	saber_karatsuba_sub_interpolation(sub_w0, ks1_res32, sub_w2, root_w0);	// ->root_w0

	memcpy(root_res_full, root_w0, root_w_len * 4);	// copy root_w0->root_res_full
	root_res_full[root_w_len] = 0;

	for (i = 0; i < sub_ab_len; i++) {
		aTT[i] = aHH[i]; //aHH
		bTT[i] = bHH[i]; //bHH
	}
	ks1_read(ks1_res32);	// -> w0

//-----

	ks1_write_int32(aTT, bTT);	//aHH, bHH
	saber_karatsuba_reduce_sub_w(ks1_res32, sub_w0);
	for (i = 0; i < sub_ab_len; i++) {	//aHS, bHS
		aTT[i] = saber_reduce_uint16(aHL[i] + aHH[i]);	//aHS
		bTT[i] = saber_reduce_uint16(bHL[i] + bHH[i]);	//bHS
	}
	ks1_read(ks1_res32); // -> w2

//-----

	ks1_write_int32(aTT, bTT); //aHS, bHS
	saber_karatsuba_reduce_sub_w(ks1_res32, sub_w2);

	for (i = 0; i < sub_ab_len; i++) {
		aTT[i] = aSL[i];	//aSL
		bTT[i] = bSL[i];	//bSL
	}

	ks1_read(ks1_res32);	// -> w1   (w0, ks1_res32, w2)->root_w2

//==== aS,bS->root_w1 ====

	ks1_write_int32(aTT, bTT); //aSL, bSL
	saber_karatsuba_sub_interpolation(sub_w0, ks1_res32, sub_w2, root_w2); // -> root_w2

	memcpy(root_res_full + root_w_len + 1, root_w2, root_w_len * sizeof(root_w2[0])); // copy root_w2->root_res_full
	for (i = 0; i < sub_ab_len; i++) {
		aTT[i] = aSH[i];
		bTT[i] = bSH[i];
	}
	ks1_read(ks1_res32);	// -> w0

//----------

	ks1_write_int32(aTT, bTT);	//aSH, bSH
	saber_karatsuba_reduce_sub_w(ks1_res32, sub_w0);
	for (i = 0; i < root_w_len; i++) {	//(比较耗时，约3580)
		//先把root_w0和root_w2减掉，最后算出来root_w1之后再加上
		root_res_full[root_ab_len + i] += (-root_w0[i] - root_w2[i]);
	}
	for (i = 0; i < sub_ab_len; i++) {	//aSS, bSS
		aTT[i] = saber_reduce_uint16(aSL[i] + aSH[i]);	//aSS
		bTT[i] = saber_reduce_uint16(bSL[i] + bSH[i]);	//bSS
	}
	ks1_read(ks1_res32);	// -> w2

//---------

	ks1_write_int32(aTT, bTT);	//aSS, bSS
	saber_karatsuba_reduce_sub_w(ks1_res32, sub_w2);
	memcpy(root_w1, sub_w0, sub_w_len * sizeof(sub_w0[0]));
	root_w1[sub_w_len] = 0;
	memcpy(root_w1 + sub_w_len + 1, sub_w2, sub_w_len * sizeof(sub_w2[0]));
	for (i = 0; i < sub_w_len; i++) {
		//先把sub_w0和sub_w2减掉，最后再加sub_w1即可
		root_w1[i + sub_ab_len] += (-sub_w0[i] - sub_w2[i]);
	}
	ks1_read(ks1_res32);	// -> w1  (w0, ks1_res32, w2)->root_w1

//=================以下不能异步了=====================
//要传给下一次的数据包括root_w1, ks1_res32, root_res_full, res (和整个的acc)

	for (i = 0; i < sub_w_len; i++) {	//并行中已减掉sub_w0和sub_w2
		root_w1[i + sub_ab_len] += (ks1_res32[i]);	// & SABER_Q_MASK);//reduce_sub_w_uint32(ks1_res32[i], REDUCE_Q);
	}

	for (i = 0; i < root_w_len; i++) {	//并行中已减掉root_w0和root_w2
		root_res_full[root_ab_len + i] += root_w1[i];
	}
//----------reduction-------
//2019.03.06 把把reduction也异步了，MatrixVectorMul_parallel

}

//void MatrixVectorMul_rounding_parallel(polyvec *a, uint16_t skpv[SABER_K][SABER_N], uint16_t res[SABER_K][SABER_N],
//		uint16_t mod,
//		int16_t transpose) {
//	uint16_t i;
//	//printf("MatrixVectorMul_parallel, transpose=%d\n", transpose);
//	const int root_res_full_len = 511;		// 511;
//	uint16_t root_res_full[root_res_full_len]; //res_full_len=511
//	if (transpose == 1) {
//		//for (i = 0; i < SABER_K; i++) {
//		i = 0;
//		saber_karatsuba_parallel_vector_in_matrix_and_rounding((uint16_t *) &a[0].vec[i], (uint16_t *) &a[1].vec[i],
//				(uint16_t *) &a[2].vec[i], skpv[0], skpv[1], skpv[2], res[i], root_res_full, NULL);
//		i = 1;
//		saber_karatsuba_parallel_vector_in_matrix_and_rounding((uint16_t *) &a[0].vec[i], (uint16_t *) &a[1].vec[i],
//				(uint16_t *) &a[2].vec[i], skpv[0], skpv[1], skpv[2], res[i], root_res_full, res[i - 1]);
//		i = 2;
//		saber_karatsuba_parallel_vector_in_matrix_and_rounding((uint16_t *) &a[0].vec[i], (uint16_t *) &a[1].vec[i],
//				(uint16_t *) &a[2].vec[i], skpv[0], skpv[1], skpv[2], res[i], root_res_full, res[i - 1]);
//		reduce_root_res_full_to_acc_and_rounding(root_res_full, res[i]);
//
//		//}
//	} else {
//		//transpose=0
//		//for (i = 0; i < SABER_K; i++) {
//		i = 0;
//		saber_karatsuba_parallel_vector_in_matrix_and_rounding((uint16_t *) &a[i].vec[0], (uint16_t *) &a[i].vec[1],
//				(uint16_t *) &a[i].vec[2], skpv[0], skpv[1], skpv[2], res[i], root_res_full, NULL);
//		i = 1;
//		saber_karatsuba_parallel_vector_in_matrix_and_rounding((uint16_t *) &a[i].vec[0], (uint16_t *) &a[i].vec[1],
//				(uint16_t *) &a[i].vec[2], skpv[0], skpv[1], skpv[2], res[i], root_res_full, res[i - 1]);
//		i = 2;
//		saber_karatsuba_parallel_vector_in_matrix_and_rounding((uint16_t *) &a[i].vec[0], (uint16_t *) &a[i].vec[1],
//				(uint16_t *) &a[i].vec[2], skpv[0], skpv[1], skpv[2], res[i], root_res_full, res[i - 1]);
//		reduce_root_res_full_to_acc_and_rounding(root_res_full, res[i]);
//		//}
//	}
//
//}

//dualcore run on CORE_MAIN
void dualcore_MatrixVectorMul_rounding_parallel(polyvec *a, polyvec *skpv, polyvec* res, uint16_t mod,
		int16_t transpose) {
	//printf("MatrixVectorMul_rounding_parallel_unrolled, transpose=%d\n", transpose);
	//skpv[SABER_K][SABER_N]
	//uint16_t res[SABER_K][SABER_N]
	uint16_t i;
	//printf("MatrixVectorMul_parallel, transpose=%d\n", transpose);
	const int root_res_full_len = 511;		// 511;
	uint16_t root_res_full[root_res_full_len]; //res_full_len=511
	if (transpose == 1) {
		dualcoreTakeSemaphore(semaphore_ext);
		i = 0;
		//printf("begin mul A[0]\n");
		saber_karatsuba_parallel_vector_in_matrix_and_rounding(a[0].vec[i].coeffs, a[1].vec[i].coeffs,
				a[2].vec[i].coeffs, skpv->vec[0].coeffs, skpv->vec[1].coeffs, skpv->vec[2].coeffs, res->vec[i].coeffs,
				root_res_full, NULL);

		dualcoreTakeSemaphore(semaphore_ext);
		i = 1;
		//printf("begin mul A[1]\n");
		saber_karatsuba_parallel_vector_in_matrix_and_rounding(a[0].vec[i].coeffs, a[1].vec[i].coeffs,
				a[2].vec[i].coeffs, skpv->vec[0].coeffs, skpv->vec[1].coeffs, skpv->vec[2].coeffs,
				res->vec[i].coeffs, root_res_full, res->vec[i - 1].coeffs);

		dualcoreTakeSemaphore(semaphore_ext);
		i = 2;
		//printf("begin mul A[2]\n");
		saber_karatsuba_parallel_vector_in_matrix_and_rounding(a[0].vec[i].coeffs, a[1].vec[i].coeffs,
				a[2].vec[i].coeffs, skpv->vec[0].coeffs, skpv->vec[1].coeffs, skpv->vec[2].coeffs,
				res->vec[i].coeffs, root_res_full, res->vec[i - 1].coeffs);
		reduce_root_res_full_to_acc_and_rounding(root_res_full, res->vec[i].coeffs);
	} else { //transpose=0
		dualcoreTakeSemaphore(semaphore_ext);
		i = 0;
		saber_karatsuba_parallel_vector_in_matrix_and_rounding(a[i].vec[0].coeffs, a[i].vec[1].coeffs,
				a[i].vec[2].coeffs, skpv->vec[0].coeffs, skpv->vec[1].coeffs, skpv->vec[2].coeffs,
				res->vec[i].coeffs, root_res_full, NULL);
		dualcoreTakeSemaphore(semaphore_ext);
		i = 1;
		saber_karatsuba_parallel_vector_in_matrix_and_rounding(a[i].vec[0].coeffs, a[i].vec[1].coeffs,
				a[i].vec[2].coeffs, skpv->vec[0].coeffs, skpv->vec[1].coeffs, skpv->vec[2].coeffs,
				res->vec[i].coeffs, root_res_full, res->vec[i - 1].coeffs);
		dualcoreTakeSemaphore(semaphore_ext);
		i = 2;
		saber_karatsuba_parallel_vector_in_matrix_and_rounding(a[i].vec[0].coeffs, a[i].vec[1].coeffs,
				a[i].vec[2].coeffs, skpv->vec[0].coeffs, skpv->vec[1].coeffs, skpv->vec[2].coeffs,
				res->vec[i].coeffs, root_res_full, res->vec[i - 1].coeffs);
		reduce_root_res_full_to_acc_and_rounding(root_res_full, res->vec[i].coeffs);
	}
}
void MatrixVectorMul_rounding_parallel(polyvec *a, polyvec *skpv, polyvec* res, uint16_t mod,
		int16_t transpose) {
	uint16_t i;
	const int root_res_full_len = 511;	// 511;
	uint16_t root_res_full[root_res_full_len];	//res_full_len=511
	if (transpose == 1) {
		i = 0;
		saber_karatsuba_parallel_vector_in_matrix_and_rounding(a[0].vec[i].coeffs, a[1].vec[i].coeffs,
				a[2].vec[i].coeffs, skpv->vec[0].coeffs, skpv->vec[1].coeffs, skpv->vec[2].coeffs, res->vec[i].coeffs,
				root_res_full, NULL);
		i = 1;
		saber_karatsuba_parallel_vector_in_matrix_and_rounding(a[0].vec[i].coeffs, a[1].vec[i].coeffs,
				a[2].vec[i].coeffs, skpv->vec[0].coeffs, skpv->vec[1].coeffs, skpv->vec[2].coeffs,
				res->vec[i].coeffs, root_res_full, res->vec[i - 1].coeffs);
		i = 2;
		saber_karatsuba_parallel_vector_in_matrix_and_rounding(a[0].vec[i].coeffs, a[1].vec[i].coeffs,
				a[2].vec[i].coeffs, skpv->vec[0].coeffs, skpv->vec[1].coeffs, skpv->vec[2].coeffs,
				res->vec[i].coeffs, root_res_full, res->vec[i - 1].coeffs);
		reduce_root_res_full_to_acc_and_rounding(root_res_full, res->vec[i].coeffs);
	} else { //transpose=0
		i = 0;
		saber_karatsuba_parallel_vector_in_matrix_and_rounding(a[i].vec[0].coeffs, a[i].vec[1].coeffs,
				a[i].vec[2].coeffs, skpv->vec[0].coeffs, skpv->vec[1].coeffs, skpv->vec[2].coeffs,
				res->vec[i].coeffs, root_res_full, NULL);
		i = 1;
		saber_karatsuba_parallel_vector_in_matrix_and_rounding(a[i].vec[0].coeffs, a[i].vec[1].coeffs,
				a[i].vec[2].coeffs, skpv->vec[0].coeffs, skpv->vec[1].coeffs, skpv->vec[2].coeffs,
				res->vec[i].coeffs, root_res_full, res->vec[i - 1].coeffs);
		i = 2;
		saber_karatsuba_parallel_vector_in_matrix_and_rounding(a[i].vec[0].coeffs, a[i].vec[1].coeffs,
				a[i].vec[2].coeffs, skpv->vec[0].coeffs, skpv->vec[1].coeffs, skpv->vec[2].coeffs,
				res->vec[i].coeffs, root_res_full, res->vec[i - 1].coeffs);
		reduce_root_res_full_to_acc_and_rounding(root_res_full, res->vec[i].coeffs);
	}
}


//====================================================================================

void reduce_root_res_full_to_acc(uint16_t* root_res_full, uint16_t* acc) {
	int i;
	const int n = SABER_N;
//	for (i = 0; i < 255; i++) { //res[0 ~ (2n-2)-n=n-2=256-2=254]，剩下res[n-1]
//		//res[i] = saber_reduce_uint16(root_res_full[i] - root_res_full[i + n]);
//		root_res_full[i] -= root_res_full[i + n];
//	}
//	for (i = 0; i < n; i++) {
//		acc[i] = ((acc[i] + root_res_full[i]) & SABER_Q_MASK);
//	}
	//两个循环合并一下
	for (i = 0; i < 255; i++) { //res[0 ~ (2n-2)-n=n-2=256-2=254]，剩下res[n-1]
		acc[i] = ((acc[i] + root_res_full[i] - root_res_full[i + n]) & SABER_Q_MASK);
	}
	acc[n - 1] = ((acc[n - 1] + root_res_full[n - 1]) & SABER_Q_MASK);
}

void MatrixVectorMul_parallel(polyvec *a, uint16_t skpv[SABER_K][SABER_N], uint16_t res[SABER_K][SABER_N], uint16_t mod,
		int16_t transpose) {
	uint16_t i;
	printf("MatrixVectorMul_parallel, transpose=%d\n", transpose);
	const int root_res_full_len = 511;		// 511;
	uint16_t root_res_full[root_res_full_len]; //res_full_len=511
	if (transpose == 1) {
		//for (i = 0; i < SABER_K; i++) {
		i = 0;
		saber_karatsuba_parallel_vector_in_matrix((uint16_t *) &a[0].vec[i], (uint16_t *) &a[1].vec[i],
				(uint16_t *) &a[2].vec[i], skpv[0], skpv[1], skpv[2], res[i], root_res_full, NULL);
		i = 1;
		saber_karatsuba_parallel_vector_in_matrix((uint16_t *) &a[0].vec[i], (uint16_t *) &a[1].vec[i],
				(uint16_t *) &a[2].vec[i], skpv[0], skpv[1], skpv[2], res[i], root_res_full, res[i - 1]);
		i = 2;
		saber_karatsuba_parallel_vector_in_matrix((uint16_t *) &a[0].vec[i], (uint16_t *) &a[1].vec[i],
				(uint16_t *) &a[2].vec[i], skpv[0], skpv[1], skpv[2], res[i], root_res_full, res[i - 1]);
		reduce_root_res_full_to_acc(root_res_full, res[i]);

		//}
	} else {
		//transpose=0
		//for (i = 0; i < SABER_K; i++) {
		i = 0;
		saber_karatsuba_parallel_vector_in_matrix((uint16_t *) &a[i].vec[0], (uint16_t *) &a[i].vec[1],
				(uint16_t *) &a[i].vec[2], skpv[0], skpv[1], skpv[2], res[i], root_res_full, NULL);
		i = 1;
		saber_karatsuba_parallel_vector_in_matrix((uint16_t *) &a[i].vec[0], (uint16_t *) &a[i].vec[1],
				(uint16_t *) &a[i].vec[2], skpv[0], skpv[1], skpv[2], res[i], root_res_full, res[i - 1]);
		i = 2;
		saber_karatsuba_parallel_vector_in_matrix((uint16_t *) &a[i].vec[0], (uint16_t *) &a[i].vec[1],
				(uint16_t *) &a[i].vec[2], skpv[0], skpv[1], skpv[2], res[i], root_res_full, res[i - 1]);
		reduce_root_res_full_to_acc(root_res_full, res[i]);
		//}
	}

}

void rounding_acc(uint16_t* acc) {
	int i;
	const int n = SABER_N;
	for (i = 0; i < n; i++) { //res[0 ~ (2n-2)-n=n-2=256-2=254]，剩下res[n-1]
		acc[i] = acc[i] + 4;
		acc[i] = acc[i] >> 3;
	}
}

//n=256, 2019.03.05写的单个相乘的版本，用于justintime, root_res_full_len = 511;
//计算结果保存在root_res_full中，需要下一步才会把root_res_full在redcue和rounding后放在acc_previous_step中
void saber_karatsuba_parallel_justintime(uint16_t* a, uint16_t* b, uint16_t* root_res_full,
		uint16_t* acc_previous, uint16_t acc_previous_need_rounding) {
	REG_WRITE_NOCHECK(RSA_MULT_MODE_REG, 15);
	int16_t n = 256;
	int16_t root_ab_len = 128;	// n / 2; //128
	int16_t root_w_len = 255;	// root_ab_len * 2 - 1; //root_w_len = 255
	int16_t root_res_full_len = 511;	// n * 2 - 1;	//res_full_len = 511，2^9

	int16_t sub_ab_len = 64;	// root_ab_len / 2;
	int16_t sub_w_len = 127;	//sub_ab_len * 2 - 1;

	int i;

	//-----------aL, bL---------------------------------
	uint16_t* aL = &a[0];
	uint16_t* aH = &a[128];  //root_ab_len
	uint16_t* bL = &b[0];
	uint16_t* bH = &b[128];

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
	uint16_t root_w0[root_w_len], root_w1[root_w_len], root_w2[root_w_len];	//最多13+9=22位
	//uint16_t root_res_full[root_res_full_len]; //res_full_len=511,无需这里置零，使用w0和w2再加一个对其赋值即可
	uint16_t aS[root_ab_len];	//len=128, a的前一半加后一半
	uint16_t bS[root_ab_len];

	uint16_t* aSL = &aS[0];
	uint16_t* aSH = &aS[64];
	uint16_t* bSL = &bS[0];
	uint16_t* bSH = &bS[64];

	int32_t aTT[sub_ab_len];	//暂存下一次的write用的数据
	int32_t bTT[sub_ab_len];
	for (i = 0; i < root_ab_len; i++) {	//aS,bS
		aS[i] = saber_reduce_uint16(aL[i] + aH[i]);
		bS[i] = saber_reduce_uint16(bL[i] + bH[i]);
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
		aTT[i] = saber_reduce_uint16(aLL[i] + aLH[i]);	//aLS
		bTT[i] = saber_reduce_uint16(bLL[i] + bLH[i]);	//bLS
	}
	uint16_t sub_w0[sub_w_len], sub_w2[sub_w_len];	//
	saber_karatsuba_reduce_sub_w(ks1_res32, sub_w0);
	ks1_read(ks1_res32);	// -> w2

	//----------

	ks1_write_int32(aTT, bTT);	//aLS, bLS
	saber_karatsuba_reduce_sub_w(ks1_res32, sub_w2);
	for (i = 0; i < sub_ab_len; i++) {
		aTT[i] = aHL[i];	//aHL
		bTT[i] = bHL[i];	//bHL
	}
	ks1_read(ks1_res32); // -> w1   (w0, ks1_res32, w2)->root_w0

	//==== aH,bH->root_w2 ====

	ks1_write_int32(aTT, bTT);	//aHL, bHL
	saber_karatsuba_sub_interpolation(sub_w0, ks1_res32, sub_w2, root_w0);	// ->root_w0
	if (acc_previous != NULL) {
		reduce_root_res_full_to_acc(root_res_full, acc_previous);
	}

	for (i = 0; i < sub_ab_len; i++) {
		aTT[i] = aHH[i]; //aHH
		bTT[i] = bHH[i]; //bHH
	}
	ks1_read(ks1_res32);	// -> w0

	//-----

	ks1_write_int32(aTT, bTT);	//aHH, bHH
	saber_karatsuba_reduce_sub_w(ks1_res32, sub_w0);

	if (acc_previous_need_rounding) {
		rounding_acc(acc_previous);
	}

	memcpy(root_res_full, root_w0, root_w_len * 4);	// copy root_w0->root_res_full
	root_res_full[root_w_len] = 0;

	for (i = 0; i < sub_ab_len; i++) {	//aHS, bHS
		aTT[i] = saber_reduce_uint16(aHL[i] + aHH[i]);	//aHS
		bTT[i] = saber_reduce_uint16(bHL[i] + bHH[i]);	//bHS
	}
	ks1_read(ks1_res32); // -> w2

	//-----

	ks1_write_int32(aTT, bTT); //aHS, bHS
	saber_karatsuba_reduce_sub_w(ks1_res32, sub_w2);

	for (i = 0; i < sub_ab_len; i++) {
		aTT[i] = aSL[i];	//aSL
		bTT[i] = bSL[i];	//bSL
	}

	ks1_read(ks1_res32);	// -> w1   (w0, ks1_res32, w2)->root_w2

	//==== aS,bS->root_w1 ====

	ks1_write_int32(aTT, bTT); //aSL, bSL
	saber_karatsuba_sub_interpolation(sub_w0, ks1_res32, sub_w2, root_w2); // -> root_w2

	memcpy(root_res_full + root_w_len + 1, root_w2, root_w_len * sizeof(root_w2[0])); // copy root_w2->root_res_full
	for (i = 0; i < sub_ab_len; i++) {
		aTT[i] = aSH[i];
		bTT[i] = bSH[i];
	}
	ks1_read(ks1_res32);	// -> w0

	//----------

	ks1_write_int32(aTT, bTT);	//aSH, bSH
	saber_karatsuba_reduce_sub_w(ks1_res32, sub_w0);
	for (i = 0; i < root_w_len; i++) {	//(比较耗时，约3580)
		//先把root_w0和root_w2减掉，最后算出来root_w1之后再加上
		root_res_full[root_ab_len + i] += (-root_w0[i] - root_w2[i]);
	}
	for (i = 0; i < sub_ab_len; i++) {	//aSS, bSS
		aTT[i] = saber_reduce_uint16(aSL[i] + aSH[i]);	//aSS
		bTT[i] = saber_reduce_uint16(bSL[i] + bSH[i]);	//bSS
	}
	ks1_read(ks1_res32);	// -> w2

	//---------

	ks1_write_int32(aTT, bTT);	//aSS, bSS
	saber_karatsuba_reduce_sub_w(ks1_res32, sub_w2);
	memcpy(root_w1, sub_w0, sub_w_len * sizeof(sub_w0[0]));
	root_w1[sub_w_len] = 0;
	memcpy(root_w1 + sub_w_len + 1, sub_w2, sub_w_len * sizeof(sub_w2[0]));
	for (i = 0; i < sub_w_len; i++) {
		//先把sub_w0和sub_w2减掉，最后再加sub_w1即可
		root_w1[i + sub_ab_len] += (-sub_w0[i] - sub_w2[i]);
	}
	ks1_read(ks1_res32);	// -> w1  (w0, ks1_res32, w2)->root_w1

	//=================以下不能异步了=====================

	///rest step.0
	for (i = 0; i < sub_w_len; i++) {	//并行中已减掉sub_w0和sub_w2
		root_w1[i + sub_ab_len] += (ks1_res32[i]);	// & SABER_Q_MASK);//reduce_sub_w_uint32(ks1_res32[i], REDUCE_Q);
	}
	// rest step.1
	for (i = 0; i < root_w_len; i++) {	//并行中已减掉root_w0和root_w2
		root_res_full[root_ab_len + i] += root_w1[i];
	}
	//----------reduction-------
//	//rest step.2
//	for (i = 0; i < 255; i++) { //res[0 ~ (2n-2)-n=n-2=256-2=254]，剩下res[n-1]
//		acc[i] += (root_res_full[i] - root_res_full[i + n], REDUCE_Q);
//	}
//	//n=256
//	acc[n - 1] += reduce_q_int32(root_res_full[n - 1], REDUCE_Q);
	//---------- rounding----------
}

//acc_previous_step是上一次的acc指针，如果是NULL就是第一次的乘法，没有上一步的acc
void saber_karatsuba_parallel_vector_in_matrix(uint16_t* a0, uint16_t* a1, uint16_t* a2,
		uint16_t* b0, uint16_t* b1, uint16_t* b2, uint16_t* acc, uint16_t* root_res_full,
		uint16_t* acc_previous_step) {
//printf("---karatsuba_parallel---\n");
	REG_WRITE_NOCHECK(RSA_MULT_MODE_REG, 15);
	const int16_t n = 256;
	int16_t root_ab_len = 128;	// n / 2; //128
	int16_t root_w_len = 255;	// root_ab_len * 2 - 1; //root_w_len = 255
//int16_t root_res_full_len = 511;	// n * 2 - 1;	//res_full_len = 511，2^9

	int16_t sub_ab_len = 64;	// root_ab_len / 2;
	int16_t sub_w_len = 127;	//sub_ab_len * 2 - 1;

	int i;

//============ a0, b0==================

//-----------aL, bL---------------------------------
	uint16_t* aL = &a0[0];
	uint16_t* aH = &a0[128];  //root_ab_len
	uint16_t* bL = &b0[0];
	uint16_t* bH = &b0[128];

	uint16_t* aLL = &a0[0];
	uint16_t* aLH = &a0[64];	//sub_ab_len
	uint16_t* bLL = &b0[0];
	uint16_t* bLH = &b0[64];	//sub_ab_len

//-----------aH, bH------------------------------
	uint16_t* aHL = &a0[128];	//sub_ab_len*2
	uint16_t* aHH = &a0[192];	//sub_ab_len*3
	uint16_t* bHL = &b0[128];
	uint16_t* bHH = &b0[192];

//=== aL,bL->root_w0 ====

	ks1_write_uint16(aLL, bLL);
	uint16_t root_w0[root_w_len], root_w1[root_w_len], root_w2[root_w_len];	//最多13+9=22位
//uint16_t root_res_full[root_res_full_len]; //res_full_len=511,无需这里置零，使用w0和w2再加一个对其赋值即可
	uint16_t aS[root_ab_len];	//len=128, a的前一半加后一半
	uint16_t bS[root_ab_len];

	uint16_t* aSL = &aS[0];
	uint16_t* aSH = &aS[64];
	uint16_t* bSL = &bS[0];
	uint16_t* bSH = &bS[64];

	int32_t aTT[sub_ab_len];	//暂存下一次的write用的数据
	int32_t bTT[sub_ab_len];
	for (i = 0; i < root_ab_len; i++) {	//aS,bS
		aS[i] = saber_reduce_uint16(aL[i] + aH[i]);
		bS[i] = saber_reduce_uint16(bL[i] + bH[i]);
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
		aTT[i] = saber_reduce_uint16(aLL[i] + aLH[i]);	//aLS
		bTT[i] = saber_reduce_uint16(bLL[i] + bLH[i]);	//bLS
	}
	uint16_t sub_w0[sub_w_len], sub_w2[sub_w_len];	//
	saber_karatsuba_reduce_sub_w(ks1_res32, sub_w0);
	ks1_read(ks1_res32);	// -> w2

//----------

	ks1_write_int32(aTT, bTT);	//aLS, bLS
	saber_karatsuba_reduce_sub_w(ks1_res32, sub_w2);
	for (i = 0; i < sub_ab_len; i++) {
		aTT[i] = aHL[i];	//aHL
		bTT[i] = bHL[i];	//bHL
	}
	ks1_read(ks1_res32); // -> w1   (w0, ks1_res32, w2)->root_w0

//==== aH,bH->root_w2 ====

	ks1_write_int32(aTT, bTT);	//aHL, bHL
	saber_karatsuba_sub_interpolation(sub_w0, ks1_res32, sub_w2, root_w0);	// ->root_w0
	if (acc_previous_step != NULL) {
		reduce_root_res_full_to_acc(root_res_full, acc_previous_step);
	}

	memcpy(root_res_full, root_w0, root_w_len * 4);	// copy root_w0->root_res_full
	root_res_full[root_w_len] = 0;

	for (i = 0; i < sub_ab_len; i++) {
		aTT[i] = aHH[i]; //aHH
		bTT[i] = bHH[i]; //bHH
	}
	ks1_read(ks1_res32);	// -> w0

//-----

	ks1_write_int32(aTT, bTT);	//aHH, bHH
	saber_karatsuba_reduce_sub_w(ks1_res32, sub_w0);
	for (i = 0; i < sub_ab_len; i++) {	//aHS, bHS
		aTT[i] = saber_reduce_uint16(aHL[i] + aHH[i]);	//aHS
		bTT[i] = saber_reduce_uint16(bHL[i] + bHH[i]);	//bHS
	}
	ks1_read(ks1_res32); // -> w2

//-----

	ks1_write_int32(aTT, bTT); //aHS, bHS
	saber_karatsuba_reduce_sub_w(ks1_res32, sub_w2);

	for (i = 0; i < sub_ab_len; i++) {
		aTT[i] = aSL[i];	//aSL
		bTT[i] = bSL[i];	//bSL
	}

	ks1_read(ks1_res32);	// -> w1   (w0, ks1_res32, w2)->root_w2

//==== aS,bS->root_w1 ====

	ks1_write_int32(aTT, bTT); //aSL, bSL
	saber_karatsuba_sub_interpolation(sub_w0, ks1_res32, sub_w2, root_w2); // -> root_w2

	memcpy(root_res_full + root_w_len + 1, root_w2, root_w_len * sizeof(root_w2[0])); // copy root_w2->root_res_full
	for (i = 0; i < sub_ab_len; i++) {
		aTT[i] = aSH[i];
		bTT[i] = bSH[i];
	}
	ks1_read(ks1_res32);	// -> w0

//----------

	ks1_write_int32(aTT, bTT);	//aSH, bSH
	saber_karatsuba_reduce_sub_w(ks1_res32, sub_w0);
	for (i = 0; i < root_w_len; i++) {	//(比较耗时，约3580)
		//先把root_w0和root_w2减掉，最后算出来root_w1之后再加上
		root_res_full[root_ab_len + i] += (-root_w0[i] - root_w2[i]);
	}
	for (i = 0; i < sub_ab_len; i++) {	//aSS, bSS
		aTT[i] = saber_reduce_uint16(aSL[i] + aSH[i]);	//aSS
		bTT[i] = saber_reduce_uint16(bSL[i] + bSH[i]);	//bSS
	}
	ks1_read(ks1_res32);	// -> w2

//---------

	ks1_write_int32(aTT, bTT);	//aSS, bSS
	saber_karatsuba_reduce_sub_w(ks1_res32, sub_w2);
	memcpy(root_w1, sub_w0, sub_w_len * sizeof(sub_w0[0]));
	root_w1[sub_w_len] = 0;
	memcpy(root_w1 + sub_w_len + 1, sub_w2, sub_w_len * sizeof(sub_w2[0]));
	for (i = 0; i < sub_w_len; i++) {
		//先把sub_w0和sub_w2减掉，最后再加sub_w1即可
		root_w1[i + sub_ab_len] += (-sub_w0[i] - sub_w2[i]);
	}

// a1, b1 pointer
	aL = &a1[0];
	aH = &a1[128];  //root_ab_len
	bL = &b1[0];
	bH = &b1[128];

//-----------aL, bL---------------------------------

	aLL = &a1[0];
	aLH = &a1[64];	//sub_ab_len
	bLL = &b1[0];
	bLH = &b1[64];	//sub_ab_len

//-----------aH, bH------------------------------
	aHL = &a1[128];	//sub_ab_len*2
	aHH = &a1[192];	//sub_ab_len*3
	bHL = &b1[128];
	bHH = &b1[192];

	for (i = 0; i < sub_ab_len; i++) {
		aTT[i] = aLL[i];
		bTT[i] = bLL[i];
	}

	ks1_read(ks1_res32);	// -> w1  (w0, ks1_res32, w2)->root_w1

//#########################################################

//============ a1, b1 ==============

//=== aL,bL->root_w0 ====

//ks1_write_uint16(aLL, bLL);
	ks1_write_int32(aTT, bTT);		//aLL,bLL已经处理过了

	for (i = 0; i < root_ab_len; i++) {	//aS,bS
		aS[i] = saber_reduce_uint16(aL[i] + aH[i]);
		bS[i] = saber_reduce_uint16(bL[i] + bH[i]);
	}
	for (i = 0; i < sub_ab_len; i++) {
		aTT[i] = aLH[i];	//aLH = &a[64];
		bTT[i] = bLH[i];	//bLH = &b[64];
	}

//==== a0, b0 rest step.0
	for (i = 0; i < sub_w_len; i++) {	//并行中已减掉sub_w0和sub_w2
		root_w1[i + sub_ab_len] += (ks1_res32[i]);	// & SABER_Q_MASK);//reduce_sub_w_uint32(ks1_res32[i], REDUCE_Q);
	}
//======================

	ks1_read(ks1_res32);	// -> w0

//-------

	ks1_write_int32(aTT, bTT);	//aLH, bLH
	for (i = 0; i < sub_ab_len; i++) {	//aLS, bLS
		aTT[i] = saber_reduce_uint16(aLL[i] + aLH[i]);	//aLS
		bTT[i] = saber_reduce_uint16(bLL[i] + bLH[i]);	//bLS
	}
	saber_karatsuba_reduce_sub_w(ks1_res32, sub_w0);

//==== a0, b0 rest step.1
	for (i = 0; i < root_w_len; i++) {	//并行中已减掉root_w0和root_w2
		root_res_full[root_ab_len + i] += root_w1[i];
	}
//=======================

	ks1_read(ks1_res32);	// -> w2

//----------

	ks1_write_int32(aTT, bTT);	//aLS, bLS
	saber_karatsuba_reduce_sub_w(ks1_res32, sub_w2);
	for (i = 0; i < sub_ab_len; i++) {
		aTT[i] = aHL[i];	//aHL
		bTT[i] = bHL[i];	//bHL
	}

//==== a0, b0 rest step.2
	for (i = 0; i < 255; i++) { //reduction,res[0 ~ (2n-2)-n=n-2=256-2=254]，剩下res[n-1]
		//res[i] = saber_reduce_uint16(root_res_full[i] - root_res_full[i + n]);
		root_res_full[i] -= root_res_full[i + n];
	}
//n=256
//res[n - 1] = saber_reduce_uint16(root_res_full[n - 1]);

	memcpy(acc, root_res_full, sizeof(root_res_full[0]) * n);

//======================

	ks1_read(ks1_res32); // -> w1   (w0, ks1_res32, w2)->root_w0

//==== aH,bH->root_w2 ====

	ks1_write_int32(aTT, bTT);	//aHL, bHL
	saber_karatsuba_sub_interpolation(sub_w0, ks1_res32, sub_w2, root_w0);	// ->root_w0

	memcpy(root_res_full, root_w0, root_w_len * 4);	// copy root_w0->root_res_full
	root_res_full[root_w_len] = 0;

	for (i = 0; i < sub_ab_len; i++) {
		aTT[i] = aHH[i]; //aHH
		bTT[i] = bHH[i]; //bHH
	}
	ks1_read(ks1_res32);	// -> w0

//-----

	ks1_write_int32(aTT, bTT);	//aHH, bHH
	saber_karatsuba_reduce_sub_w(ks1_res32, sub_w0);
	for (i = 0; i < sub_ab_len; i++) {	//aHS, bHS
		aTT[i] = saber_reduce_uint16(aHL[i] + aHH[i]);	//aHS
		bTT[i] = saber_reduce_uint16(bHL[i] + bHH[i]);	//bHS
	}
	ks1_read(ks1_res32); // -> w2

//-----

	ks1_write_int32(aTT, bTT); //aHS, bHS
	saber_karatsuba_reduce_sub_w(ks1_res32, sub_w2);

	for (i = 0; i < sub_ab_len; i++) {
		aTT[i] = aSL[i];	//aSL
		bTT[i] = bSL[i];	//bSL
	}

	ks1_read(ks1_res32);	// -> w1   (w0, ks1_res32, w2)->root_w2

//==== aS,bS->root_w1 ====

	ks1_write_int32(aTT, bTT); //aSL, bSL
	saber_karatsuba_sub_interpolation(sub_w0, ks1_res32, sub_w2, root_w2); // -> root_w2

	memcpy(root_res_full + root_w_len + 1, root_w2, root_w_len * sizeof(root_w2[0])); // copy root_w2->root_res_full
	for (i = 0; i < sub_ab_len; i++) {
		aTT[i] = aSH[i];
		bTT[i] = bSH[i];
	}
	ks1_read(ks1_res32);	// -> w0

//----------

	ks1_write_int32(aTT, bTT);	//aSH, bSH
	saber_karatsuba_reduce_sub_w(ks1_res32, sub_w0);
	for (i = 0; i < root_w_len; i++) {	//(比较耗时，约3580)
		//先把root_w0和root_w2减掉，最后算出来root_w1之后再加上
		root_res_full[root_ab_len + i] += (-root_w0[i] - root_w2[i]);
	}
	for (i = 0; i < sub_ab_len; i++) {	//aSS, bSS
		aTT[i] = saber_reduce_uint16(aSL[i] + aSH[i]);	//aSS
		bTT[i] = saber_reduce_uint16(bSL[i] + bSH[i]);	//bSS
	}
	ks1_read(ks1_res32);	// -> w2

//---------

	ks1_write_int32(aTT, bTT);	//aSS, bSS
	saber_karatsuba_reduce_sub_w(ks1_res32, sub_w2);
	memcpy(root_w1, sub_w0, sub_w_len * sizeof(sub_w0[0]));
	root_w1[sub_w_len] = 0;
	memcpy(root_w1 + sub_w_len + 1, sub_w2, sub_w_len * sizeof(sub_w2[0]));
	for (i = 0; i < sub_w_len; i++) {
		//先把sub_w0和sub_w2减掉，最后再加sub_w1即可
		root_w1[i + sub_ab_len] += (-sub_w0[i] - sub_w2[i]);
	}

// a2, b2 pointer
	aL = &a2[0];
	aH = &a2[128];  //root_ab_len
	bL = &b2[0];
	bH = &b2[128];

//-----------aL, bL---------------------------------

	aLL = &a2[0];
	aLH = &a2[64];	//sub_ab_len
	bLL = &b2[0];
	bLH = &b2[64];	//sub_ab_len

//-----------aH, bH------------------------------
	aHL = &a2[128];	//sub_ab_len*2
	aHH = &a2[192];	//sub_ab_len*3
	bHL = &b2[128];
	bHH = &b2[192];

	for (i = 0; i < sub_ab_len; i++) {
		aTT[i] = aLL[i];
		bTT[i] = bLL[i];
	}

	ks1_read(ks1_res32);	// -> w1  (w0, ks1_res32, w2)->root_w1

//#########################################################

//============ a2, b2 ==============

//=== aL,bL->root_w0 ====

//ks1_write_uint16(aLL, bLL);
	ks1_write_int32(aTT, bTT);		//aLL,bLL已经处理过了

	for (i = 0; i < root_ab_len; i++) {	//aS,bS
		aS[i] = saber_reduce_uint16(aL[i] + aH[i]);
		bS[i] = saber_reduce_uint16(bL[i] + bH[i]);
	}
	for (i = 0; i < sub_ab_len; i++) {
		aTT[i] = aLH[i];	//aLH = &a[64];
		bTT[i] = bLH[i];	//bLH = &b[64];
	}

//==== a1, b1 rest step.0
	for (i = 0; i < sub_w_len; i++) {	//并行中已减掉sub_w0和sub_w2
		root_w1[i + sub_ab_len] += (ks1_res32[i]);	// & SABER_Q_MASK);//reduce_sub_w_uint32(ks1_res32[i], REDUCE_Q);
	}
//======================

	ks1_read(ks1_res32);	// -> w0

//-------

	ks1_write_int32(aTT, bTT);	//aLH, bLH
	for (i = 0; i < sub_ab_len; i++) {	//aLS, bLS
		aTT[i] = saber_reduce_uint16(aLL[i] + aLH[i]);	//aLS
		bTT[i] = saber_reduce_uint16(bLL[i] + bLH[i]);	//bLS
	}
	saber_karatsuba_reduce_sub_w(ks1_res32, sub_w0);

//==== a1, b1 rest step.1
	for (i = 0; i < root_w_len; i++) {	//并行中已减掉root_w0和root_w2
		root_res_full[root_ab_len + i] += root_w1[i];
	}
//=======================

	ks1_read(ks1_res32);	// -> w2

//----------

	ks1_write_int32(aTT, bTT);	//aLS, bLS
	saber_karatsuba_reduce_sub_w(ks1_res32, sub_w2);
	for (i = 0; i < sub_ab_len; i++) {
		aTT[i] = aHL[i];	//aHL
		bTT[i] = bHL[i];	//bHL
	}

//==== a1, b1 rest step.2
	for (i = 0; i < 255; i++) { //reduction,res[0 ~ (2n-2)-n=n-2=256-2=254]，剩下res[n-1]
		//这里可能有问题
		//res[i] = saber_reduce_uint16(root_res_full[i] - root_res_full[i + n]);
		//root_res_full[i] = saber_reduce_uint16(root_res_full[i] - root_res_full[i + n]);
		root_res_full[i] -= root_res_full[i + n];
	}
//n=256
//res[n - 1] = saber_reduce_uint16(root_res_full[n - 1]);		// a1, b1, over

	for (i = 0; i < n; i++) {
		//acc[i] += res[i];
		acc[i] += root_res_full[i];
	}

//======================

	ks1_read(ks1_res32); // -> w1   (w0, ks1_res32, w2)->root_w0

//==== aH,bH->root_w2 ====

	ks1_write_int32(aTT, bTT);	//aHL, bHL
	saber_karatsuba_sub_interpolation(sub_w0, ks1_res32, sub_w2, root_w0);	// ->root_w0

	memcpy(root_res_full, root_w0, root_w_len * 4);	// copy root_w0->root_res_full
	root_res_full[root_w_len] = 0;

	for (i = 0; i < sub_ab_len; i++) {
		aTT[i] = aHH[i]; //aHH
		bTT[i] = bHH[i]; //bHH
	}
	ks1_read(ks1_res32);	// -> w0

//-----

	ks1_write_int32(aTT, bTT);	//aHH, bHH
	saber_karatsuba_reduce_sub_w(ks1_res32, sub_w0);
	for (i = 0; i < sub_ab_len; i++) {	//aHS, bHS
		aTT[i] = saber_reduce_uint16(aHL[i] + aHH[i]);	//aHS
		bTT[i] = saber_reduce_uint16(bHL[i] + bHH[i]);	//bHS
	}
	ks1_read(ks1_res32); // -> w2

//-----

	ks1_write_int32(aTT, bTT); //aHS, bHS
	saber_karatsuba_reduce_sub_w(ks1_res32, sub_w2);

	for (i = 0; i < sub_ab_len; i++) {
		aTT[i] = aSL[i];	//aSL
		bTT[i] = bSL[i];	//bSL
	}

	ks1_read(ks1_res32);	// -> w1   (w0, ks1_res32, w2)->root_w2

//==== aS,bS->root_w1 ====

	ks1_write_int32(aTT, bTT); //aSL, bSL
	saber_karatsuba_sub_interpolation(sub_w0, ks1_res32, sub_w2, root_w2); // -> root_w2

	memcpy(root_res_full + root_w_len + 1, root_w2, root_w_len * sizeof(root_w2[0])); // copy root_w2->root_res_full
	for (i = 0; i < sub_ab_len; i++) {
		aTT[i] = aSH[i];
		bTT[i] = bSH[i];
	}
	ks1_read(ks1_res32);	// -> w0

//----------

	ks1_write_int32(aTT, bTT);	//aSH, bSH
	saber_karatsuba_reduce_sub_w(ks1_res32, sub_w0);
	for (i = 0; i < root_w_len; i++) {	//(比较耗时，约3580)
		//先把root_w0和root_w2减掉，最后算出来root_w1之后再加上
		root_res_full[root_ab_len + i] += (-root_w0[i] - root_w2[i]);
	}
	for (i = 0; i < sub_ab_len; i++) {	//aSS, bSS
		aTT[i] = saber_reduce_uint16(aSL[i] + aSH[i]);	//aSS
		bTT[i] = saber_reduce_uint16(bSL[i] + bSH[i]);	//bSS
	}
	ks1_read(ks1_res32);	// -> w2

//---------

	ks1_write_int32(aTT, bTT);	//aSS, bSS
	saber_karatsuba_reduce_sub_w(ks1_res32, sub_w2);
	memcpy(root_w1, sub_w0, sub_w_len * sizeof(sub_w0[0]));
	root_w1[sub_w_len] = 0;
	memcpy(root_w1 + sub_w_len + 1, sub_w2, sub_w_len * sizeof(sub_w2[0]));
	for (i = 0; i < sub_w_len; i++) {
		//先把sub_w0和sub_w2减掉，最后再加sub_w1即可
		root_w1[i + sub_ab_len] += (-sub_w0[i] - sub_w2[i]);
	}
	ks1_read(ks1_res32);	// -> w1  (w0, ks1_res32, w2)->root_w1

//=================以下不能异步了=====================
//要传给下一次的数据包括root_w1, ks1_res32, root_res_full, res (和整个的acc)

	for (i = 0; i < sub_w_len; i++) {	//并行中已减掉sub_w0和sub_w2
		root_w1[i + sub_ab_len] += (ks1_res32[i]);	// & SABER_Q_MASK);//reduce_sub_w_uint32(ks1_res32[i], REDUCE_Q);
	}

	for (i = 0; i < root_w_len; i++) {	//并行中已减掉root_w0和root_w2
		root_res_full[root_ab_len + i] += root_w1[i];
	}
//----------reduction-------
//2019.03.06 把把reduction也异步了，MatrixVectorMul_parallel

//	for (i = 0; i < 255; i++) { //res[0 ~ (2n-2)-n=n-2=256-2=254]，剩下res[n-1]
//		//res[i] = saber_reduce_uint16(root_res_full[i] - root_res_full[i + n]);
//		root_res_full[i] -= root_res_full[i + n];
//	}
//	//n=256
//	//res[n - 1] = saber_reduce_uint16(root_res_full[n - 1]);
//
//	for (i = 0; i < n; i++) {
//		acc[i] = ((acc[i] + root_res_full[i]) & SABER_Q_MASK);
//	}

}

void saber_karatsuba_parallel_vector(uint16_t* a0, uint16_t* a1, uint16_t* a2,
		uint16_t* b0, uint16_t* b1, uint16_t* b2, uint16_t* acc) {
//printf("---karatsuba_parallel---\n");
	REG_WRITE_NOCHECK(RSA_MULT_MODE_REG, 15);
	const int16_t n = 256;
	int16_t root_ab_len = 128;	// n / 2; //128
	int16_t root_w_len = 255;	// root_ab_len * 2 - 1; //root_w_len = 255
	int16_t root_res_full_len = 511;	// n * 2 - 1;	//res_full_len = 511，2^9

	int16_t sub_ab_len = 64;	// root_ab_len / 2;
	int16_t sub_w_len = 127;	//sub_ab_len * 2 - 1;

	int i;

//============ a0, b0==================

//-----------aL, bL---------------------------------
	uint16_t* aL = &a0[0];
	uint16_t* aH = &a0[128];  //root_ab_len
	uint16_t* bL = &b0[0];
	uint16_t* bH = &b0[128];

	uint16_t* aLL = &a0[0];
	uint16_t* aLH = &a0[64];	//sub_ab_len
	uint16_t* bLL = &b0[0];
	uint16_t* bLH = &b0[64];	//sub_ab_len

//-----------aH, bH------------------------------
	uint16_t* aHL = &a0[128];	//sub_ab_len*2
	uint16_t* aHH = &a0[192];	//sub_ab_len*3
	uint16_t* bHL = &b0[128];
	uint16_t* bHH = &b0[192];

//=== aL,bL->root_w0 ====

	ks1_write_uint16(aLL, bLL);
	uint16_t root_w0[root_w_len], root_w1[root_w_len], root_w2[root_w_len];	//最多13+9=22位
	uint16_t root_res_full[root_res_full_len]; //res_full_len=511,无需这里置零，使用w0和w2再加一个对其赋值即可
	uint16_t aS[root_ab_len];	//len=128, a的前一半加后一半
	uint16_t bS[root_ab_len];

	uint16_t* aSL = &aS[0];
	uint16_t* aSH = &aS[64];
	uint16_t* bSL = &bS[0];
	uint16_t* bSH = &bS[64];

	int32_t aTT[sub_ab_len];	//暂存下一次的write用的数据
	int32_t bTT[sub_ab_len];
	for (i = 0; i < root_ab_len; i++) {	//aS,bS
		aS[i] = saber_reduce_uint16(aL[i] + aH[i]);
		bS[i] = saber_reduce_uint16(bL[i] + bH[i]);
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
		aTT[i] = saber_reduce_uint16(aLL[i] + aLH[i]);	//aLS
		bTT[i] = saber_reduce_uint16(bLL[i] + bLH[i]);	//bLS
	}
	uint16_t sub_w0[sub_w_len], sub_w2[sub_w_len];	//
	saber_karatsuba_reduce_sub_w(ks1_res32, sub_w0);
	ks1_read(ks1_res32);	// -> w2

//----------

	ks1_write_int32(aTT, bTT);	//aLS, bLS
	saber_karatsuba_reduce_sub_w(ks1_res32, sub_w2);
	for (i = 0; i < sub_ab_len; i++) {
		aTT[i] = aHL[i];	//aHL
		bTT[i] = bHL[i];	//bHL
	}
	ks1_read(ks1_res32); // -> w1   (w0, ks1_res32, w2)->root_w0

//==== aH,bH->root_w2 ====

	ks1_write_int32(aTT, bTT);	//aHL, bHL
	saber_karatsuba_sub_interpolation(sub_w0, ks1_res32, sub_w2, root_w0);	// ->root_w0

	memcpy(root_res_full, root_w0, root_w_len * 4);	// copy root_w0->root_res_full
	root_res_full[root_w_len] = 0;

	for (i = 0; i < sub_ab_len; i++) {
		aTT[i] = aHH[i]; //aHH
		bTT[i] = bHH[i]; //bHH
	}
	ks1_read(ks1_res32);	// -> w0

//-----

	ks1_write_int32(aTT, bTT);	//aHH, bHH
	saber_karatsuba_reduce_sub_w(ks1_res32, sub_w0);
	for (i = 0; i < sub_ab_len; i++) {	//aHS, bHS
		aTT[i] = saber_reduce_uint16(aHL[i] + aHH[i]);	//aHS
		bTT[i] = saber_reduce_uint16(bHL[i] + bHH[i]);	//bHS
	}
	ks1_read(ks1_res32); // -> w2

//-----

	ks1_write_int32(aTT, bTT); //aHS, bHS
	saber_karatsuba_reduce_sub_w(ks1_res32, sub_w2);

	for (i = 0; i < sub_ab_len; i++) {
		aTT[i] = aSL[i];	//aSL
		bTT[i] = bSL[i];	//bSL
	}

	ks1_read(ks1_res32);	// -> w1   (w0, ks1_res32, w2)->root_w2

//==== aS,bS->root_w1 ====

	ks1_write_int32(aTT, bTT); //aSL, bSL
	saber_karatsuba_sub_interpolation(sub_w0, ks1_res32, sub_w2, root_w2); // -> root_w2

	memcpy(root_res_full + root_w_len + 1, root_w2, root_w_len * sizeof(root_w2[0])); // copy root_w2->root_res_full
	for (i = 0; i < sub_ab_len; i++) {
		aTT[i] = aSH[i];
		bTT[i] = bSH[i];
	}
	ks1_read(ks1_res32);	// -> w0

//----------

	ks1_write_int32(aTT, bTT);	//aSH, bSH
	saber_karatsuba_reduce_sub_w(ks1_res32, sub_w0);
	for (i = 0; i < root_w_len; i++) {	//(比较耗时，约3580)
		//先把root_w0和root_w2减掉，最后算出来root_w1之后再加上
		root_res_full[root_ab_len + i] += (-root_w0[i] - root_w2[i]);
	}
	for (i = 0; i < sub_ab_len; i++) {	//aSS, bSS
		aTT[i] = saber_reduce_uint16(aSL[i] + aSH[i]);	//aSS
		bTT[i] = saber_reduce_uint16(bSL[i] + bSH[i]);	//bSS
	}
	ks1_read(ks1_res32);	// -> w2

//---------

	ks1_write_int32(aTT, bTT);	//aSS, bSS
	saber_karatsuba_reduce_sub_w(ks1_res32, sub_w2);
	memcpy(root_w1, sub_w0, sub_w_len * sizeof(sub_w0[0]));
	root_w1[sub_w_len] = 0;
	memcpy(root_w1 + sub_w_len + 1, sub_w2, sub_w_len * sizeof(sub_w2[0]));
	for (i = 0; i < sub_w_len; i++) {
		//先把sub_w0和sub_w2减掉，最后再加sub_w1即可
		root_w1[i + sub_ab_len] += (-sub_w0[i] - sub_w2[i]);
	}

// a1, b1 pointer
	aL = &a1[0];
	aH = &a1[128];  //root_ab_len
	bL = &b1[0];
	bH = &b1[128];

//-----------aL, bL---------------------------------

	aLL = &a1[0];
	aLH = &a1[64];	//sub_ab_len
	bLL = &b1[0];
	bLH = &b1[64];	//sub_ab_len

//-----------aH, bH------------------------------
	aHL = &a1[128];	//sub_ab_len*2
	aHH = &a1[192];	//sub_ab_len*3
	bHL = &b1[128];
	bHH = &b1[192];

	for (i = 0; i < sub_ab_len; i++) {
		aTT[i] = aLL[i];
		bTT[i] = bLL[i];
	}

	ks1_read(ks1_res32);	// -> w1  (w0, ks1_res32, w2)->root_w1

//#########################################################

//============ a1, b1 ==============

//=== aL,bL->root_w0 ====

//ks1_write_uint16(aLL, bLL);
	ks1_write_int32(aTT, bTT);		//aLL,bLL已经处理过了

	for (i = 0; i < root_ab_len; i++) {	//aS,bS
		aS[i] = saber_reduce_uint16(aL[i] + aH[i]);
		bS[i] = saber_reduce_uint16(bL[i] + bH[i]);
	}
	for (i = 0; i < sub_ab_len; i++) {
		aTT[i] = aLH[i];	//aLH = &a[64];
		bTT[i] = bLH[i];	//bLH = &b[64];
	}

//==== a0, b0 rest step.0
	for (i = 0; i < sub_w_len; i++) {	//并行中已减掉sub_w0和sub_w2
		root_w1[i + sub_ab_len] += (ks1_res32[i]);	// & SABER_Q_MASK);//reduce_sub_w_uint32(ks1_res32[i], REDUCE_Q);
	}
//======================

	ks1_read(ks1_res32);	// -> w0

//-------

	ks1_write_int32(aTT, bTT);	//aLH, bLH
	for (i = 0; i < sub_ab_len; i++) {	//aLS, bLS
		aTT[i] = saber_reduce_uint16(aLL[i] + aLH[i]);	//aLS
		bTT[i] = saber_reduce_uint16(bLL[i] + bLH[i]);	//bLS
	}
	saber_karatsuba_reduce_sub_w(ks1_res32, sub_w0);

//==== a0, b0 rest step.1
	for (i = 0; i < root_w_len; i++) {	//并行中已减掉root_w0和root_w2
		root_res_full[root_ab_len + i] += root_w1[i];
	}
//=======================

	ks1_read(ks1_res32);	// -> w2

//----------

	ks1_write_int32(aTT, bTT);	//aLS, bLS
	saber_karatsuba_reduce_sub_w(ks1_res32, sub_w2);
	for (i = 0; i < sub_ab_len; i++) {
		aTT[i] = aHL[i];	//aHL
		bTT[i] = bHL[i];	//bHL
	}

//==== a0, b0 rest step.2
	for (i = 0; i < 255; i++) { //reduction,res[0 ~ (2n-2)-n=n-2=256-2=254]，剩下res[n-1]
		//res[i] = saber_reduce_uint16(root_res_full[i] - root_res_full[i + n]);
		root_res_full[i] -= root_res_full[i + n];
	}
//n=256
//res[n - 1] = saber_reduce_uint16(root_res_full[n - 1]);

	memcpy(acc, root_res_full, sizeof(root_res_full[0]) * n);

//======================

	ks1_read(ks1_res32); // -> w1   (w0, ks1_res32, w2)->root_w0

//==== aH,bH->root_w2 ====

	ks1_write_int32(aTT, bTT);	//aHL, bHL
	saber_karatsuba_sub_interpolation(sub_w0, ks1_res32, sub_w2, root_w0);	// ->root_w0

	memcpy(root_res_full, root_w0, root_w_len * 4);	// copy root_w0->root_res_full
	root_res_full[root_w_len] = 0;

	for (i = 0; i < sub_ab_len; i++) {
		aTT[i] = aHH[i]; //aHH
		bTT[i] = bHH[i]; //bHH
	}
	ks1_read(ks1_res32);	// -> w0

//-----

	ks1_write_int32(aTT, bTT);	//aHH, bHH
	saber_karatsuba_reduce_sub_w(ks1_res32, sub_w0);
	for (i = 0; i < sub_ab_len; i++) {	//aHS, bHS
		aTT[i] = saber_reduce_uint16(aHL[i] + aHH[i]);	//aHS
		bTT[i] = saber_reduce_uint16(bHL[i] + bHH[i]);	//bHS
	}
	ks1_read(ks1_res32); // -> w2

//-----

	ks1_write_int32(aTT, bTT); //aHS, bHS
	saber_karatsuba_reduce_sub_w(ks1_res32, sub_w2);

	for (i = 0; i < sub_ab_len; i++) {
		aTT[i] = aSL[i];	//aSL
		bTT[i] = bSL[i];	//bSL
	}

	ks1_read(ks1_res32);	// -> w1   (w0, ks1_res32, w2)->root_w2

//==== aS,bS->root_w1 ====

	ks1_write_int32(aTT, bTT); //aSL, bSL
	saber_karatsuba_sub_interpolation(sub_w0, ks1_res32, sub_w2, root_w2); // -> root_w2

	memcpy(root_res_full + root_w_len + 1, root_w2, root_w_len * sizeof(root_w2[0])); // copy root_w2->root_res_full
	for (i = 0; i < sub_ab_len; i++) {
		aTT[i] = aSH[i];
		bTT[i] = bSH[i];
	}
	ks1_read(ks1_res32);	// -> w0

//----------

	ks1_write_int32(aTT, bTT);	//aSH, bSH
	saber_karatsuba_reduce_sub_w(ks1_res32, sub_w0);
	for (i = 0; i < root_w_len; i++) {	//(比较耗时，约3580)
		//先把root_w0和root_w2减掉，最后算出来root_w1之后再加上
		root_res_full[root_ab_len + i] += (-root_w0[i] - root_w2[i]);
	}
	for (i = 0; i < sub_ab_len; i++) {	//aSS, bSS
		aTT[i] = saber_reduce_uint16(aSL[i] + aSH[i]);	//aSS
		bTT[i] = saber_reduce_uint16(bSL[i] + bSH[i]);	//bSS
	}
	ks1_read(ks1_res32);	// -> w2

//---------

	ks1_write_int32(aTT, bTT);	//aSS, bSS
	saber_karatsuba_reduce_sub_w(ks1_res32, sub_w2);
	memcpy(root_w1, sub_w0, sub_w_len * sizeof(sub_w0[0]));
	root_w1[sub_w_len] = 0;
	memcpy(root_w1 + sub_w_len + 1, sub_w2, sub_w_len * sizeof(sub_w2[0]));
	for (i = 0; i < sub_w_len; i++) {
		//先把sub_w0和sub_w2减掉，最后再加sub_w1即可
		root_w1[i + sub_ab_len] += (-sub_w0[i] - sub_w2[i]);
	}

// a2, b2 pointer
	aL = &a2[0];
	aH = &a2[128];  //root_ab_len
	bL = &b2[0];
	bH = &b2[128];

//-----------aL, bL---------------------------------

	aLL = &a2[0];
	aLH = &a2[64];	//sub_ab_len
	bLL = &b2[0];
	bLH = &b2[64];	//sub_ab_len

//-----------aH, bH------------------------------
	aHL = &a2[128];	//sub_ab_len*2
	aHH = &a2[192];	//sub_ab_len*3
	bHL = &b2[128];
	bHH = &b2[192];

	for (i = 0; i < sub_ab_len; i++) {
		aTT[i] = aLL[i];
		bTT[i] = bLL[i];
	}

	ks1_read(ks1_res32);	// -> w1  (w0, ks1_res32, w2)->root_w1

//#########################################################

//============ a2, b2 ==============

//=== aL,bL->root_w0 ====

//ks1_write_uint16(aLL, bLL);
	ks1_write_int32(aTT, bTT);		//aLL,bLL已经处理过了

	for (i = 0; i < root_ab_len; i++) {	//aS,bS
		aS[i] = saber_reduce_uint16(aL[i] + aH[i]);
		bS[i] = saber_reduce_uint16(bL[i] + bH[i]);
	}
	for (i = 0; i < sub_ab_len; i++) {
		aTT[i] = aLH[i];	//aLH = &a[64];
		bTT[i] = bLH[i];	//bLH = &b[64];
	}

//==== a1, b1 rest step.0
	for (i = 0; i < sub_w_len; i++) {	//并行中已减掉sub_w0和sub_w2
		root_w1[i + sub_ab_len] += (ks1_res32[i]);	// & SABER_Q_MASK);//reduce_sub_w_uint32(ks1_res32[i], REDUCE_Q);
	}
//======================

	ks1_read(ks1_res32);	// -> w0

//-------

	ks1_write_int32(aTT, bTT);	//aLH, bLH
	for (i = 0; i < sub_ab_len; i++) {	//aLS, bLS
		aTT[i] = saber_reduce_uint16(aLL[i] + aLH[i]);	//aLS
		bTT[i] = saber_reduce_uint16(bLL[i] + bLH[i]);	//bLS
	}
	saber_karatsuba_reduce_sub_w(ks1_res32, sub_w0);

//==== a1, b1 rest step.1
	for (i = 0; i < root_w_len; i++) {	//并行中已减掉root_w0和root_w2
		root_res_full[root_ab_len + i] += root_w1[i];
	}
//=======================

	ks1_read(ks1_res32);	// -> w2

//----------

	ks1_write_int32(aTT, bTT);	//aLS, bLS
	saber_karatsuba_reduce_sub_w(ks1_res32, sub_w2);
	for (i = 0; i < sub_ab_len; i++) {
		aTT[i] = aHL[i];	//aHL
		bTT[i] = bHL[i];	//bHL
	}

//==== a1, b1 rest step.2
	for (i = 0; i < 255; i++) { //reduction,res[0 ~ (2n-2)-n=n-2=256-2=254]，剩下res[n-1]
		//这里可能有问题
		//res[i] = saber_reduce_uint16(root_res_full[i] - root_res_full[i + n]);
		//root_res_full[i] = saber_reduce_uint16(root_res_full[i] - root_res_full[i + n]);
		root_res_full[i] -= root_res_full[i + n];
	}
//n=256
//res[n - 1] = saber_reduce_uint16(root_res_full[n - 1]);		// a1, b1, over

	for (i = 0; i < n; i++) {
		//acc[i] += res[i];
		acc[i] += root_res_full[i];
	}

//======================

	ks1_read(ks1_res32); // -> w1   (w0, ks1_res32, w2)->root_w0

//==== aH,bH->root_w2 ====

	ks1_write_int32(aTT, bTT);	//aHL, bHL
	saber_karatsuba_sub_interpolation(sub_w0, ks1_res32, sub_w2, root_w0);	// ->root_w0

	memcpy(root_res_full, root_w0, root_w_len * 4);	// copy root_w0->root_res_full
	root_res_full[root_w_len] = 0;

	for (i = 0; i < sub_ab_len; i++) {
		aTT[i] = aHH[i]; //aHH
		bTT[i] = bHH[i]; //bHH
	}
	ks1_read(ks1_res32);	// -> w0

//-----

	ks1_write_int32(aTT, bTT);	//aHH, bHH
	saber_karatsuba_reduce_sub_w(ks1_res32, sub_w0);
	for (i = 0; i < sub_ab_len; i++) {	//aHS, bHS
		aTT[i] = saber_reduce_uint16(aHL[i] + aHH[i]);	//aHS
		bTT[i] = saber_reduce_uint16(bHL[i] + bHH[i]);	//bHS
	}
	ks1_read(ks1_res32); // -> w2

//-----

	ks1_write_int32(aTT, bTT); //aHS, bHS
	saber_karatsuba_reduce_sub_w(ks1_res32, sub_w2);

	for (i = 0; i < sub_ab_len; i++) {
		aTT[i] = aSL[i];	//aSL
		bTT[i] = bSL[i];	//bSL
	}

	ks1_read(ks1_res32);	// -> w1   (w0, ks1_res32, w2)->root_w2

//==== aS,bS->root_w1 ====

	ks1_write_int32(aTT, bTT); //aSL, bSL
	saber_karatsuba_sub_interpolation(sub_w0, ks1_res32, sub_w2, root_w2); // -> root_w2

	memcpy(root_res_full + root_w_len + 1, root_w2, root_w_len * sizeof(root_w2[0])); // copy root_w2->root_res_full
	for (i = 0; i < sub_ab_len; i++) {
		aTT[i] = aSH[i];
		bTT[i] = bSH[i];
	}
	ks1_read(ks1_res32);	// -> w0

//----------

	ks1_write_int32(aTT, bTT);	//aSH, bSH
	saber_karatsuba_reduce_sub_w(ks1_res32, sub_w0);
	for (i = 0; i < root_w_len; i++) {	//(比较耗时，约3580)
		//先把root_w0和root_w2减掉，最后算出来root_w1之后再加上
		root_res_full[root_ab_len + i] += (-root_w0[i] - root_w2[i]);
	}
	for (i = 0; i < sub_ab_len; i++) {	//aSS, bSS
		aTT[i] = saber_reduce_uint16(aSL[i] + aSH[i]);	//aSS
		bTT[i] = saber_reduce_uint16(bSL[i] + bSH[i]);	//bSS
	}
	ks1_read(ks1_res32);	// -> w2

//---------

	ks1_write_int32(aTT, bTT);	//aSS, bSS
	saber_karatsuba_reduce_sub_w(ks1_res32, sub_w2);
	memcpy(root_w1, sub_w0, sub_w_len * sizeof(sub_w0[0]));
	root_w1[sub_w_len] = 0;
	memcpy(root_w1 + sub_w_len + 1, sub_w2, sub_w_len * sizeof(sub_w2[0]));
	for (i = 0; i < sub_w_len; i++) {
		//先把sub_w0和sub_w2减掉，最后再加sub_w1即可
		root_w1[i + sub_ab_len] += (-sub_w0[i] - sub_w2[i]);
	}
	ks1_read(ks1_res32);	// -> w1  (w0, ks1_res32, w2)->root_w1

//=================以下不能异步了=====================
//要传给下一次的数据包括root_w1, ks1_res32, root_res_full, res (和整个的acc)

	for (i = 0; i < sub_w_len; i++) {	//并行中已减掉sub_w0和sub_w2
		root_w1[i + sub_ab_len] += (ks1_res32[i]);	// & SABER_Q_MASK);//reduce_sub_w_uint32(ks1_res32[i], REDUCE_Q);
	}

	for (i = 0; i < root_w_len; i++) {	//并行中已减掉root_w0和root_w2
		root_res_full[root_ab_len + i] += root_w1[i];
	}
//----------reduction-------

	for (i = 0; i < 255; i++) { //res[0 ~ (2n-2)-n=n-2=256-2=254]，剩下res[n-1]
		//res[i] = saber_reduce_uint16(root_res_full[i] - root_res_full[i + n]);
		root_res_full[i] -= root_res_full[i + n];
	}
//n=256
//res[n - 1] = saber_reduce_uint16(root_res_full[n - 1]);

	for (i = 0; i < n; i++) {
		acc[i] = ((acc[i] + root_res_full[i]) & SABER_Q_MASK);
	}

}
