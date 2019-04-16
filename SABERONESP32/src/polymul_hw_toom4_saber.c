#include "polymul_hw_common.h"
#include "polymul_hw.h"
#include "polymul_hw_speed.h"

void toom4_saber_10bits(uint16_t* a, uint16_t* b, uint16_t* result) {
	//printf("---toom4_saber_10bits---\n");
	REG_WRITE_NOCHECK(RSA_MULT_MODE_REG, 15);

	int16_t n = 256;
	int16_t i;
	uint16_t sub_len = 64;	//64,n / 4
	uint16_t sub_len2 = 128;
	uint16_t sub_len3 = 192;
	uint16_t sub_len4 = 256;
	uint16_t sub_len5 = 320;
	uint16_t sub_len6 = 384;
	uint16_t w_len = 127;	//2 * sub_len - 1;

	uint16_t w1[w_len], w2[w_len], w3[w_len], w4[w_len], w5[w_len], w6[w_len], w7[w_len];			//w_len=128
	//uint32_t ks1_res32[w_len];

	uint16_t at1, at2;  //temp for w
	uint16_t bt1, bt2;

	uint16_t* ainf = &a[sub_len3];	//A(inf),B(inf)->w1, inf:infinite
	uint16_t* binf = &b[sub_len3];
	uint16_t a2[sub_len], b2[sub_len];	//A(2)*B(2)->w2
	uint16_t a1[sub_len], b1[sub_len];	//A(1)*B(1)->w3
	uint16_t a_1[sub_len], b_1[sub_len];	//A(-1)*B(-1)->w4
	uint16_t ahalf[sub_len], bhalf[sub_len];	//A(1/2)*B(1/2)->w5
	uint16_t a_half[sub_len], b_half[sub_len];	//A(-1/2)*B(-1/2)->w6
	uint16_t* a0 = &a[0];	//A0,B0->w7
	uint16_t* b0 = &b[0];

	int16_t inv3 = 43691;
	int16_t inv9 = 36409;
	int16_t inv15 = 61167;

	int16_t int45 = 45;
	int16_t int30 = 30;

	const uint16_t toom4_precision_mask = SABER_P * 8 - 1;

	for (i = 0; i < sub_len; i++) {

		//---w5,w6=a±half*b±half=A(±1/2)*B(±1/2)=(8A0+2A2)±(4A1+A3)---
		at1 = a0[i] << 2;	//4*A0
		bt1 = b0[i] << 2;
		at1 = at1 + a[sub_len2 + i];	//4A0+A2
		bt1 = bt1 + b[sub_len2 + i];
		at1 = at1 << 1;	//8A0 + 2A2
		bt1 = bt1 << 1;
		at2 = a[sub_len + i];	//A1
		bt2 = b[sub_len + i];
		at2 = at2 << 2;	//4A1
		bt2 = bt2 << 2;
		at2 = at2 + a[sub_len3 + i];	//4A1 + A3
		bt2 = bt2 + b[sub_len3 + i];
		ahalf[i] = at1 + at2;	//(8A0 + 2A2) + (4A1 + A3) = 8A0 + 4A1 + 2A2 + A3
		bhalf[i] = bt1 + bt2;
		a_half[i] = at1 - at2;	//(8A0 + 2A2) - (4A1 + A3) = 8A0 - 4A1 - 2A2 + A3
		b_half[i] = bt1 - bt2;

		//--- w3,w4=a(±1)*b(±1)=(A2+A0)±(A3+A1) ---
		at1 = a[sub_len2 + i] + a[i];	//A2 + A0
		bt1 = b[sub_len2 + i] + b[i];
		at2 = a[sub_len3 + i] + a[sub_len + i];	//A3 + A1
		bt2 = b[sub_len3 + i] + b[sub_len + i];
		a1[i] = at1 + at2; //(A2 + A0) + (A3 + A1) = A0 + A1 + A2 + A3
		b1[i] = bt1 + bt2;
		a_1[i] = at1 - at2;	//(A2 + A0) - (A3 + A1) = A0 - A1 + A2 - A3
		b_1[i] = bt1 - bt2;

		//--- w2=a2*b2=8A3+4A2+2A1+A0 ---
		a2[i] = a[sub_len3 + i] << 1;	//2A3
		b2[i] = b[sub_len3 + i] << 1;
		a2[i] += a[sub_len2 + i];	//2A3+A2
		b2[i] += b[sub_len2 + i];
		a2[i] <<= 1;	//4A3+2A2
		b2[i] <<= 1;
		a2[i] += a[sub_len + i];	//4A3+2A2+A1
		b2[i] += b[sub_len + i];
		a2[i] <<= 1;	//8A3+2A2+2A1
		b2[i] <<= 1;
		a2[i] += a[i];	//8A3+2A2+2A1+A0
		b2[i] += b[i];

		a2[i] = a2[i] & toom4_precision_mask;
		b2[i] = b2[i] & toom4_precision_mask;

		a1[i] = a1[i] & toom4_precision_mask;
		b1[i] = b1[i] & toom4_precision_mask;

		a_1[i] = a_1[i] & toom4_precision_mask;
		b_1[i] = b_1[i] & toom4_precision_mask;

		ahalf[i] = ahalf[i] & toom4_precision_mask;
		bhalf[i] = bhalf[i] & toom4_precision_mask;

		a_half[i] = a_half[i] & toom4_precision_mask;
		b_half[i] = b_half[i] & toom4_precision_mask;

	}
	const uint16_t mask = toom4_precision_mask;
	ks1_saber_uint16_uint16(ainf, binf, w1, mask);
	ks1_saber_uint16_uint16(a2, b2, w2, mask);
	ks1_saber_uint16_uint16(a1, b1, w3, mask);
	ks1_saber_uint16_uint16(a_1, b_1, w4, mask);
	ks1_saber_uint16_uint16(ahalf, bhalf, w5, mask);
	ks1_saber_uint16_uint16(a_half, b_half, w6, mask);
	ks1_saber_uint16_uint16(a0, b0, w7, mask);

	for (i = 0; i < w_len; i++) {
		//step10: w2=w2+w5
		w2[i] = w2[i] + w5[i];	//w2 <- w2+w5
		//step11: w6=w6-w5
		w6[i] = w6[i] - w5[i];	// w6 <- w6-w5
		//step12: w4=(w4-w3)/2
		w4[i] = w4[i] - w3[i];
		w4[i] = w4[i] >> 1;	//toom4_w_mul_inv(w4[i], inv2);
		//step13: w5=w5-w1-64*w7
		w5[i] = w5[i] - w1[i] - (w7[i] << 6); // w5 <- w5-64*w7
		//step14: w3=w3+w4
		w3[i] = w3[i] + w4[i]; //w3 <- w3+w4
		//step15: w5=2*w5+w6
		w5[i] = w6[i] + (w5[i] << 1); //w5 <- 2*w5+w6
		//step16: w2=w2-65*w3
		w2[i] = w2[i] - (w3[i] << 6) - w3[i];
		//step17: w3=w3-w7-w1
		w3[i] = w3[i] - w7[i] - w1[i];
		//step18: w2=w2+45*w3
		w2[i] = w2[i] + w3[i] * int45; //w2 <- w2+45*w3
		//step19: w5=(w5-8*w3)/24
		w5[i] = w5[i] - (w3[i] << 3); //w5 <- w5-8*w3
		w5[i] = w5[i] * inv3;
		w5[i] = w5[i] >> 3;
		//step20: w6=w6+w2
		w6[i] = w2[i] + w6[i]; //w6 <- w6+w2
		//step21: w2=(w2+16*w4)/18
		w2[i] = w2[i] + (w4[i] << 4); //w2 <- w2+16*w4
		w2[i] = w2[i] * inv9;
		w2[i] = w2[i] >> 1;
		//step22. w3 <- w3-w5
		w3[i] = w3[i] - w5[i];
		//step23: w4=-(w4+w2/2)
		w4[i] = -(w4[i] + w2[i]);
		//step24:  w6=(30*w2-w6)/60
		w6[i] = w2[i] * int30 - w6[i]; //w6 <- 30*w2-w6
		w6[i] = w6[i] * inv15; //*inv15 / 4
		w6[i] = w6[i] >> 2;
		//step25: w2=w2-w6
		w2[i] = w2[i] - w6[i];
	}

	uint16_t res_full_len = 511; //2 * n - 1;
	uint16_t res_full[res_full_len];

	//sub_len=64,w_len=127,res_full_len=511
//	memset(res_full, 0, sizeof(res_full[0]) * res_full_len);
//	for (i = 0; i < w_len; i++) {
//		res_full[i] += w7[i];
//		res_full[sub_len + i] += w6[i];
//		res_full[sub_len2 + i] += w5[i];
//		res_full[sub_len3 + i] += w4[i];
//		res_full[sub_len4 + i] += w3[i];
//		res_full[sub_len5 + i] += w2[i];
//		res_full[sub_len6 + i] += w1[i];
//	}

//	poly_print_uint16(w1, "w1", SABER_N, 0);
//	poly_print_uint16(w2, "w2", SABER_N, 0);
//	poly_print_uint16(w3, "w3", SABER_N, 0);
//	poly_print_uint16(w4, "w4", SABER_N, 0);
//	poly_print_uint16(w5, "w5", SABER_N, 0);
//	poly_print_uint16(w6, "w6", SABER_N, 0);
//	poly_print_uint16(w7, "w7", SABER_N, 0);

	memcpy(res_full, w7, w_len * 2);
	memcpy(res_full + sub_len2, w5, w_len * 2);
	memcpy(res_full + sub_len4, w3, w_len * 2);
	memcpy(res_full + sub_len6, w1, w_len * 2);
	res_full[sub_len2 - 1] = 0;
	res_full[sub_len4 - 1] = 0;
	res_full[sub_len6 - 1] = 0;
	for (i = 0; i < w_len; i++) {
		res_full[sub_len + i] += w6[i];
		res_full[sub_len3 + i] += w4[i];
		res_full[sub_len5 + i] += w2[i];
	}

	//---------------reduction-------
	const uint16_t result_mask = SABER_P - 1;
	//i从n到2 * n - 1;即从256到511
	for (i = 0; i < 255; i++) { //res[0 ~ (2n-2)-n=n-2=256-2=254]，剩下res[n-1]
		//result[i] = toom4_reduce_int32_full(res_full[i] - res_full[i + n]);
		result[i] = (res_full[i] - res_full[i + n]) & result_mask;
	}
	//result[n - 1] = toom4_reduce_int32_full(res_full[n - 1]);
	result[n - 1] = res_full[n - 1] & result_mask;
}

void toom4_saber_10bits_reduce_w(uint16_t* w, uint32_t* ks1_res32, uint16_t mask) {
	int i;
	//const uint16_t mask = SABER_Q - 1;
	for (i = 0; i < 127; i++) {	//2 * 64 - 1
		w[i] = ks1_res32[i] & mask;
	}
}

void toom4_saber_10bits_parallel(uint16_t* a, uint16_t* b, uint16_t* result) {
	//printf("---toom4_saber_10bits_parallel---\n");
	//REG_WRITE_NOCHECK(RSA_MULT_MODE_REG, 15);//这个很神奇啊，不注释掉就不对了？？？
	const uint16_t toom4_precision_mask = SABER_P * 8 - 1;

	int16_t n = 256;
	int16_t i;
	uint16_t sub_len = 64;	//64,n / 4
	uint16_t sub_len2 = 128;
	uint16_t sub_len3 = 192;
	uint16_t sub_len4 = 256;
	uint16_t sub_len5 = 320;
	uint16_t sub_len6 = 384;
	uint16_t w_len = 127;	//2 * sub_len - 1;

	uint16_t at1, at2;  //temp for w
	uint16_t bt1, bt2;

	int16_t inv3 = 43691;
	int16_t inv9 = 36409;
	int16_t inv15 = 61167;

	int16_t int45 = 45;
	int16_t int30 = 30;

	uint16_t* ainf = &a[sub_len3];	//A(inf),B(inf)->w1, inf:infinite
	uint16_t* binf = &b[sub_len3];

	//---w1---
	ks1_write_uint16(ainf, binf);

	uint16_t a2[sub_len], b2[sub_len];	//A(2)*B(2)->w2
	uint16_t a1[sub_len], b1[sub_len];	//A(1)*B(1)->w3
	uint16_t a_1[sub_len], b_1[sub_len];	//A(-1)*B(-1)->w4
	uint16_t ahalf[sub_len], bhalf[sub_len];	//A(1/2)*B(1/2)->w5
	uint16_t a_half[sub_len], b_half[sub_len];	//A(-1/2)*B(-1/2)->w6
	uint16_t* a0 = &a[0];	//A0,B0->w7
	uint16_t* b0 = &b[0];

	uint16_t w1[w_len], w2[w_len], w3[w_len], w4[w_len], w5[w_len], w6[w_len], w7[w_len];			//w_len=128

	uint16_t w5_copy[w_len];	//备份w5给step10,11
	uint16_t w2_copy[w_len];	//在step21之前备份w5给step20

	uint32_t ks1_res32[127];	//存放从寄存器读取的ks1结果
	int32_t aT[sub_len]; //32-bit alignment
	int32_t bT[sub_len];

	for (i = 0; i < sub_len; i++) {
		//--- w3,w4=a(±1)*b(±1)=(A2+A0)±(A3+A1) ---
		at1 = a[sub_len2 + i] + a[i];	//A2 + A0
		bt1 = b[sub_len2 + i] + b[i];
		at2 = a[sub_len3 + i] + a[sub_len + i];	//A3 + A1
		bt2 = b[sub_len3 + i] + b[sub_len + i];
		a1[i] = at1 + at2; //(A2 + A0) + (A3 + A1) = A0 + A1 + A2 + A3
		b1[i] = bt1 + bt2;
		a_1[i] = at1 - at2;	//(A2 + A0) - (A3 + A1) = A0 - A1 + A2 - A3
		b_1[i] = bt1 - bt2;

		a1[i] = a1[i] & toom4_precision_mask;
		b1[i] = b1[i] & toom4_precision_mask;
		a_1[i] = a_1[i] & toom4_precision_mask;
		b_1[i] = b_1[i] & toom4_precision_mask;

		aT[i] = a1[i];
		bT[i] = b1[i];
	}
	ks1_read(ks1_res32);

	//---w3---
	ks1_write_int32(aT, bT);	//a1, b1;
	toom4_saber_10bits_reduce_w(w1, ks1_res32, toom4_precision_mask);	//->w1
	//poly_print_uint16(w1, "w1", 256, 0);
	for (i = 0; i < sub_len; i++) {
		//---w5,w6=a±half*b±half=A(±1/2)*B(±1/2)=(8A0+2A2)±(4A1+A3)---
		at1 = a0[i] << 2;	//4*A0
		bt1 = b0[i] << 2;
		at1 = at1 + a[sub_len2 + i];	//4A0+A2
		bt1 = bt1 + b[sub_len2 + i];
		at1 = at1 << 1;	//8A0 + 2A2
		bt1 = bt1 << 1;
		at2 = a[sub_len + i];	//A1
		bt2 = b[sub_len + i];
		at2 = at2 << 2;	//4A1
		bt2 = bt2 << 2;
		at2 = at2 + a[sub_len3 + i];	//4A1 + A3
		bt2 = bt2 + b[sub_len3 + i];
		ahalf[i] = at1 + at2;	//(8A0 + 2A2) + (4A1 + A3) = 8A0 + 4A1 + 2A2 + A3
		bhalf[i] = bt1 + bt2;
		a_half[i] = at1 - at2;	//(8A0 + 2A2) - (4A1 + A3) = 8A0 - 4A1 - 2A2 + A3
		b_half[i] = bt1 - bt2;

		ahalf[i] = ahalf[i] & toom4_precision_mask;
		bhalf[i] = bhalf[i] & toom4_precision_mask;
		a_half[i] = a_half[i] & toom4_precision_mask;
		b_half[i] = b_half[i] & toom4_precision_mask;

		aT[i] = a_1[i];
		bT[i] = b_1[i];
	}
	ks1_read(ks1_res32);

	//---w4---
	ks1_write_int32(aT, bT);	//(a_1, b_1);
	toom4_saber_10bits_reduce_w(w3, ks1_res32, toom4_precision_mask);	//->w3
	for (i = 0; i < sub_len; i++) {
		//--- w2=a2*b2=8A3+4A2+2A1+A0 ---
		a2[i] = a[sub_len3 + i] << 1;	//2A3
		b2[i] = b[sub_len3 + i] << 1;
		a2[i] += a[sub_len2 + i];	//2A3+A2
		b2[i] += b[sub_len2 + i];
		a2[i] <<= 1;	//4A3+2A2
		b2[i] <<= 1;
		a2[i] += a[sub_len + i];	//4A3+2A2+A1
		b2[i] += b[sub_len + i];
		a2[i] <<= 1;	//8A3+2A2+2A1
		b2[i] <<= 1;
		a2[i] += a[i];	//8A3+2A2+2A1+A0
		b2[i] += b[i];

		a2[i] = a2[i] & toom4_precision_mask;
		b2[i] = b2[i] & toom4_precision_mask;

		aT[i] = ahalf[i];
		bT[i] = bhalf[i];
	}
	const int16_t res_full_len = 2 * n - 1;
	uint16_t res_full[res_full_len];
	memset(res_full, 0, sizeof(res_full[0]) * res_full_len);
	for (i = 0; i < w_len; i++) {
		res_full[sub_len6 + i] += w1[i];
	}
	ks1_read(ks1_res32);

	//---w5---
	ks1_write_int32(aT, bT);	//(ahalf, bhalf);
	toom4_saber_10bits_reduce_w(w4, ks1_res32, toom4_precision_mask);	//->w4
	//[step12,14]
	for (i = 0; i < w_len; i++) {
		//step12: w4=(w4-w3)/2
		w4[i] = w4[i] - w3[i];
		w4[i] = w4[i] >> 1;
		//step14: w3=w3+w4
		w3[i] = w3[i] + w4[i]; //w3 <- w3+w4
	}
	for (i = 0; i < sub_len; i++) {
		aT[i] = a2[i];
		bT[i] = b2[i];
	}
	ks1_read(ks1_res32);

	//---w2---
	ks1_write_int32(aT, bT);	//(a2, b2);
	toom4_saber_10bits_reduce_w(w5, ks1_res32, toom4_precision_mask); //->w5
	//[step13.1,copy_w5]
	for (i = 0; i < w_len; i++) {
		w5_copy[i] = w5[i];	//copy_w5
		//step13.1:w5=w5-w1;[w5=w5-w1-64*w7]
		w5[i] = w5[i] - w1[i];
	}
	for (i = 0; i < sub_len; i++) {
		aT[i] = a0[i];
		bT[i] = b0[i];
	}
	ks1_read(ks1_res32);

	//--w7---
	ks1_write_int32(aT, bT);	//(a0, b0);
	toom4_saber_10bits_reduce_w(w2, ks1_res32, toom4_precision_mask);		//->w2
	//[step10,16,17.1]
	for (i = 0; i < w_len; i++) {
		//step10: w2=w2+w5
		w2[i] = w2[i] + w5_copy[i];	//w2 <- w2+w5
		//step16: w2=w2-65*w3
		w2[i] = w2[i] - (w3[i] << 6) - w3[i];
		//step17.1:w3=w3-w1; [w3=w3-w7-w1]
		w3[i] = w3[i] - w1[i];
	}
	for (i = 0; i < sub_len; i++) {
		aT[i] = a_half[i];
		bT[i] = b_half[i];
	}
	ks1_read(ks1_res32);

	//---w6---
	ks1_write_int32(aT, bT);	//(a_half, b_half);
	toom4_saber_10bits_reduce_w(w7, ks1_res32, toom4_precision_mask);		//->w7
	for (i = 0; i < w_len; i++) {
		res_full[i] += w7[i];
	}
	//[step13.2,17.2,18,21,23]  w2_copy
	for (i = 0; i < w_len; i++) {
		//step13.2"w5=w5-64*w7;[w5=w5-w1-64*w7]
		w5[i] = w5[i] - (w7[i] << 6);
		//step17.2:w3=w3-w7; [w3=w3-w7-w1]
		w3[i] = w3[i] - w7[i];
		//step18: w2=w2+45*w3
		w2[i] = w2[i] + w3[i] * int45; //w2 <- w2+45*w3
		w2_copy[i] = w2[i]; //w2_copy
		//step21: w2=(w2+16*w4)/18
		w2[i] = w2[i] + (w4[i] << 4); //w2 <- w2+16*w4
		w2[i] = w2[i] * inv9;
		w2[i] = w2[i] >> 1;
		//step23: w4=-(w4+w2/2)
		w4[i] = -(w4[i] + w2[i]); //->w4_final
	}
	for (i = 0; i < w_len; i++) {
		res_full[sub_len3 + i] += w4[i];
	}
	ks1_read(ks1_res32);

	//====以下是不能并行的=====
	toom4_saber_10bits_reduce_w(w6, ks1_res32, toom4_precision_mask);		//->w6

	for (i = 0; i < w_len; i++) {
		//step11: w6=w6-w5
		w6[i] = w6[i] - w5_copy[i];	// w6 <- w6-w5

		//step15: w5=2*w5+w6
		w5[i] = w6[i] + (w5[i] << 1); //w5 <- 2*w5+w6

		//step19: w5=(w5-8*w3)/24
		w5[i] = w5[i] - (w3[i] << 3); //w5 <- w5-8*w3
		w5[i] = w5[i] * inv3;
		w5[i] = w5[i] >> 3;

		//step20: w6=w6+w2
		w6[i] = w2_copy[i] + w6[i]; //w6 <- w6+w2

		//step22. w3 <- w3-w5
		w3[i] = w3[i] - w5[i];

		//step24:  w6=(30*w2-w6)/60
		w6[i] = w2[i] * int30 - w6[i]; //w6 <- 30*w2-w6
		w6[i] = w6[i] * inv15; //*inv15 / 4
		w6[i] = w6[i] >> 2;

		//step25: w2=w2-w6
		w2[i] = w2[i] - w6[i];
	}

	for (i = 0; i < w_len; i++) { //w1,w7,w4已经计算过了
		res_full[sub_len + i] += w6[i];
		res_full[sub_len2 + i] += w5[i];
		res_full[sub_len4 + i] += w3[i];
		res_full[sub_len5 + i] += w2[i];
	}
	//---------------reduction-------
	const uint16_t result_mask = SABER_P - 1;
	//i从n到2 * n - 1;即从256到511
	for (i = 0; i < 255; i++) { //res[0 ~ (2n-2)-n=n-2=256-2=254]，剩下res[n-1]
		//result[i] = toom4_reduce_int32_full(res_full[i] - res_full[i + n]);
		result[i] = (res_full[i] - res_full[i + n]) & result_mask;
	}
	//result[n - 1] = toom4_reduce_int32_full(res_full[n - 1]);
	result[n - 1] = res_full[n - 1] & result_mask;
}

void toom4_saber_10bits_parallel_vector_bak(uint16_t* va0, uint16_t* va1, uint16_t* va2,
		uint16_t* vb0, uint16_t* vb1, uint16_t* vb2, uint16_t* acc) {
	//printf("---toom4_saber_10bits_parallel---\n");
	REG_WRITE_NOCHECK(RSA_MULT_MODE_REG, 15);
	const uint16_t toom4_precision_mask = SABER_P * 8 - 1;
	const uint16_t result_mask = SABER_P - 1;

	const int16_t n = 256;
	int16_t i;
	uint16_t sub_len = 64;	//64,n / 4
	uint16_t sub_len2 = 128;
	uint16_t sub_len3 = 192;
	uint16_t sub_len4 = 256;
	uint16_t sub_len5 = 320;
	uint16_t sub_len6 = 384;
	uint16_t w_len = 127;	//2 * sub_len - 1;

	uint16_t at1, at2;  //temp for w
	uint16_t bt1, bt2;

	int16_t inv3 = 43691;
	int16_t inv9 = 36409;
	int16_t inv15 = 61167;

	int16_t int45 = 45;
	int16_t int30 = 30;

	//====== va0*vb0 ===========

	uint16_t* a = &va0[0];
	uint16_t* b = &vb0[0];

	uint16_t* ainf = &a[sub_len3];	//A(inf),B(inf)->w1, inf:infinite
	uint16_t* binf = &b[sub_len3];

	//---w1---
	ks1_write_uint16(ainf, binf);

	memset(acc, 0, sizeof(acc[0]) * SABER_N);

	uint16_t a2[sub_len], b2[sub_len];	//A(2)*B(2)->w2
	uint16_t a1[sub_len], b1[sub_len];	//A(1)*B(1)->w3
	uint16_t a_1[sub_len], b_1[sub_len];	//A(-1)*B(-1)->w4
	uint16_t ahalf[sub_len], bhalf[sub_len];	//A(1/2)*B(1/2)->w5
	uint16_t a_half[sub_len], b_half[sub_len];	//A(-1/2)*B(-1/2)->w6
	uint16_t* a0 = &a[0];	//A0,B0->w7
	uint16_t* b0 = &b[0];

	uint16_t w1[w_len], w2[w_len], w3[w_len], w4[w_len], w5[w_len], w6[w_len], w7[w_len];			//w_len=128

	uint16_t w5_copy[w_len];	//备份w5给step10,11
	uint16_t w2_copy[w_len];	//在step21之前备份w5给step20

	uint32_t ks1_res32[127];	//存放从寄存器读取的ks1结果
	int32_t aT[sub_len]; //32-bit alignment
	int32_t bT[sub_len];

	for (i = 0; i < sub_len; i++) {
		//--- w3,w4=a(±1)*b(±1)=(A2+A0)±(A3+A1) ---
		at1 = a[sub_len2 + i] + a[i];	//A2 + A0
		bt1 = b[sub_len2 + i] + b[i];
		at2 = a[sub_len3 + i] + a[sub_len + i];	//A3 + A1
		bt2 = b[sub_len3 + i] + b[sub_len + i];
		a1[i] = at1 + at2; //(A2 + A0) + (A3 + A1) = A0 + A1 + A2 + A3
		b1[i] = bt1 + bt2;
		a_1[i] = at1 - at2;	//(A2 + A0) - (A3 + A1) = A0 - A1 + A2 - A3
		b_1[i] = bt1 - bt2;

		a1[i] = a1[i] & toom4_precision_mask;
		b1[i] = b1[i] & toom4_precision_mask;
		a_1[i] = a_1[i] & toom4_precision_mask;
		b_1[i] = b_1[i] & toom4_precision_mask;

		aT[i] = a1[i];
		bT[i] = b1[i];
	}
	ks1_read(ks1_res32);

	//---w3---
	ks1_write_int32(aT, bT);	//a1, b1;
	toom4_saber_10bits_reduce_w(w1, ks1_res32, toom4_precision_mask);	//->w1
	//poly_print_uint16(w1, "w1", 256, 0);
	for (i = 0; i < sub_len; i++) {
		//---w5,w6=a±half*b±half=A(±1/2)*B(±1/2)=(8A0+2A2)±(4A1+A3)---
		at1 = a0[i] << 2;	//4*A0
		bt1 = b0[i] << 2;
		at1 = at1 + a[sub_len2 + i];	//4A0+A2
		bt1 = bt1 + b[sub_len2 + i];
		at1 = at1 << 1;	//8A0 + 2A2
		bt1 = bt1 << 1;
		at2 = a[sub_len + i];	//A1
		bt2 = b[sub_len + i];
		at2 = at2 << 2;	//4A1
		bt2 = bt2 << 2;
		at2 = at2 + a[sub_len3 + i];	//4A1 + A3
		bt2 = bt2 + b[sub_len3 + i];
		ahalf[i] = at1 + at2;	//(8A0 + 2A2) + (4A1 + A3) = 8A0 + 4A1 + 2A2 + A3
		bhalf[i] = bt1 + bt2;
		a_half[i] = at1 - at2;	//(8A0 + 2A2) - (4A1 + A3) = 8A0 - 4A1 - 2A2 + A3
		b_half[i] = bt1 - bt2;

		ahalf[i] = ahalf[i] & toom4_precision_mask;
		bhalf[i] = bhalf[i] & toom4_precision_mask;
		a_half[i] = a_half[i] & toom4_precision_mask;
		b_half[i] = b_half[i] & toom4_precision_mask;

		aT[i] = a_1[i];
		bT[i] = b_1[i];
	}
	ks1_read(ks1_res32);

	//---w4---
	ks1_write_int32(aT, bT);	//(a_1, b_1);
	toom4_saber_10bits_reduce_w(w3, ks1_res32, toom4_precision_mask);	//->w3
	for (i = 0; i < sub_len; i++) {
		//--- w2=a2*b2=8A3+4A2+2A1+A0 ---
		a2[i] = a[sub_len3 + i] << 1;	//2A3
		b2[i] = b[sub_len3 + i] << 1;
		a2[i] += a[sub_len2 + i];	//2A3+A2
		b2[i] += b[sub_len2 + i];
		a2[i] <<= 1;	//4A3+2A2
		b2[i] <<= 1;
		a2[i] += a[sub_len + i];	//4A3+2A2+A1
		b2[i] += b[sub_len + i];
		a2[i] <<= 1;	//8A3+2A2+2A1
		b2[i] <<= 1;
		a2[i] += a[i];	//8A3+2A2+2A1+A0
		b2[i] += b[i];

		a2[i] = a2[i] & toom4_precision_mask;
		b2[i] = b2[i] & toom4_precision_mask;

		aT[i] = ahalf[i];
		bT[i] = bhalf[i];
	}
	const int16_t res_full_len = 2 * n - 1;
	uint16_t res_full[res_full_len];
	memset(res_full, 0, sizeof(res_full[0]) * res_full_len);
	for (i = 0; i < w_len; i++) {
		res_full[sub_len6 + i] += w1[i];
	}
	ks1_read(ks1_res32);

	//---w5---
	ks1_write_int32(aT, bT);	//(ahalf, bhalf);
	toom4_saber_10bits_reduce_w(w4, ks1_res32, toom4_precision_mask);	//->w4
	//[step12,14]
	for (i = 0; i < w_len; i++) {
		//step12: w4=(w4-w3)/2
		w4[i] = w4[i] - w3[i];
		w4[i] = w4[i] >> 1;
		//step14: w3=w3+w4
		w3[i] = w3[i] + w4[i]; //w3 <- w3+w4
	}
	for (i = 0; i < sub_len; i++) {
		aT[i] = a2[i];
		bT[i] = b2[i];
	}
	ks1_read(ks1_res32);

	//---w2---
	ks1_write_int32(aT, bT);	//(a2, b2);
	toom4_saber_10bits_reduce_w(w5, ks1_res32, toom4_precision_mask); //->w5
	//[step13.1,copy_w5]
	for (i = 0; i < w_len; i++) {
		w5_copy[i] = w5[i];	//copy_w5
		//step13.1:w5=w5-w1;[w5=w5-w1-64*w7]
		w5[i] = w5[i] - w1[i];
	}
	for (i = 0; i < sub_len; i++) {
		aT[i] = a0[i];
		bT[i] = b0[i];
	}
	ks1_read(ks1_res32);

	//--w7---
	ks1_write_int32(aT, bT);	//(a0, b0);
	toom4_saber_10bits_reduce_w(w2, ks1_res32, toom4_precision_mask);		//->w2
	//[step10,16,17.1]
	for (i = 0; i < w_len; i++) {
		//step10: w2=w2+w5
		w2[i] = w2[i] + w5_copy[i];	//w2 <- w2+w5
		//step16: w2=w2-65*w3
		w2[i] = w2[i] - (w3[i] << 6) - w3[i];
		//step17.1:w3=w3-w1; [w3=w3-w7-w1]
		w3[i] = w3[i] - w1[i];
	}
	for (i = 0; i < sub_len; i++) {
		aT[i] = a_half[i];
		bT[i] = b_half[i];
	}
	ks1_read(ks1_res32);

	//---w6---
	ks1_write_int32(aT, bT);	//(a_half, b_half);
	toom4_saber_10bits_reduce_w(w7, ks1_res32, toom4_precision_mask);		//->w7
	for (i = 0; i < w_len; i++) {
		res_full[i] += w7[i];
	}
	//[step13.2,17.2,18,21,23]  w2_copy
	for (i = 0; i < w_len; i++) {
		//step13.2"w5=w5-64*w7;[w5=w5-w1-64*w7]
		w5[i] = w5[i] - (w7[i] << 6);
		//step17.2:w3=w3-w7; [w3=w3-w7-w1]
		w3[i] = w3[i] - w7[i];
		//step18: w2=w2+45*w3
		w2[i] = w2[i] + w3[i] * int45; //w2 <- w2+45*w3
		w2_copy[i] = w2[i]; //w2_copy
		//step21: w2=(w2+16*w4)/18
		w2[i] = w2[i] + (w4[i] << 4); //w2 <- w2+16*w4
		w2[i] = w2[i] * inv9;
		w2[i] = w2[i] >> 1;
		//step23: w4=-(w4+w2/2)
		w4[i] = -(w4[i] + w2[i]); //->w4_final
	}
	for (i = 0; i < w_len; i++) {
		res_full[sub_len3 + i] += w4[i];
	}

	ks1_read(ks1_res32);

	//====== va1*vb1 ===========

	a = &va1[0];
	b = &vb1[0];

	ainf = &a[sub_len3];	//A(inf),B(inf)->w1, inf:infinite
	binf = &b[sub_len3];

	//---w1---
	ks1_write_uint16(ainf, binf);
	//----- previous final step(1 of 2)
	toom4_saber_10bits_reduce_w(w6, ks1_res32, toom4_precision_mask);		//->w6
	//----- previous final step(1 of 2)

	a0 = &a[0];	//A0,B0->w7
	b0 = &b[0];

	for (i = 0; i < sub_len; i++) {
		//--- w3,w4=a(±1)*b(±1)=(A2+A0)±(A3+A1) ---
		at1 = a[sub_len2 + i] + a[i];	//A2 + A0
		bt1 = b[sub_len2 + i] + b[i];
		at2 = a[sub_len3 + i] + a[sub_len + i];	//A3 + A1
		bt2 = b[sub_len3 + i] + b[sub_len + i];
		a1[i] = at1 + at2; //(A2 + A0) + (A3 + A1) = A0 + A1 + A2 + A3
		b1[i] = bt1 + bt2;
		a_1[i] = at1 - at2;	//(A2 + A0) - (A3 + A1) = A0 - A1 + A2 - A3
		b_1[i] = bt1 - bt2;

		a1[i] = a1[i] & toom4_precision_mask;
		b1[i] = b1[i] & toom4_precision_mask;
		a_1[i] = a_1[i] & toom4_precision_mask;
		b_1[i] = b_1[i] & toom4_precision_mask;

		aT[i] = a1[i];
		bT[i] = b1[i];
	}

	ks1_read(ks1_res32);

	//---w3---
	ks1_write_int32(aT, bT);	//a1, b1;
	//----- previous final step(2 of 3)
	for (i = 0; i < w_len; i++) {
		//step11: w6=w6-w5
		w6[i] = w6[i] - w5_copy[i];	// w6 <- w6-w5

		//step15: w5=2*w5+w6
		w5[i] = w6[i] + (w5[i] << 1); //w5 <- 2*w5+w6

		//step19: w5=(w5-8*w3)/24
		w5[i] = w5[i] - (w3[i] << 3); //w5 <- w5-8*w3
		w5[i] = w5[i] * inv3;
		w5[i] = w5[i] >> 3;

		//step20: w6=w6+w2
		w6[i] = w2_copy[i] + w6[i]; //w6 <- w6+w2

		//step22. w3 <- w3-w5
		w3[i] = w3[i] - w5[i];

		//step24:  w6=(30*w2-w6)/60
		w6[i] = w2[i] * int30 - w6[i]; //w6 <- 30*w2-w6
		w6[i] = w6[i] * inv15; //*inv15 / 4
		w6[i] = w6[i] >> 2;

		//step25: w2=w2-w6
		w2[i] = w2[i] - w6[i];
	}
	//----- previous final step(2 of 3)

	toom4_saber_10bits_reduce_w(w1, ks1_res32, toom4_precision_mask);	//->w1
	for (i = 0; i < sub_len; i++) {
		//---w5,w6=a±half*b±half=A(±1/2)*B(±1/2)=(8A0+2A2)±(4A1+A3)---
		at1 = a0[i] << 2;	//4*A0
		bt1 = b0[i] << 2;
		at1 = at1 + a[sub_len2 + i];	//4A0+A2
		bt1 = bt1 + b[sub_len2 + i];
		at1 = at1 << 1;	//8A0 + 2A2
		bt1 = bt1 << 1;
		at2 = a[sub_len + i];	//A1
		bt2 = b[sub_len + i];
		at2 = at2 << 2;	//4A1
		bt2 = bt2 << 2;
		at2 = at2 + a[sub_len3 + i];	//4A1 + A3
		bt2 = bt2 + b[sub_len3 + i];
		ahalf[i] = at1 + at2;	//(8A0 + 2A2) + (4A1 + A3) = 8A0 + 4A1 + 2A2 + A3
		bhalf[i] = bt1 + bt2;
		a_half[i] = at1 - at2;	//(8A0 + 2A2) - (4A1 + A3) = 8A0 - 4A1 - 2A2 + A3
		b_half[i] = bt1 - bt2;

		ahalf[i] = ahalf[i] & toom4_precision_mask;
		bhalf[i] = bhalf[i] & toom4_precision_mask;
		a_half[i] = a_half[i] & toom4_precision_mask;
		b_half[i] = b_half[i] & toom4_precision_mask;

		aT[i] = a_1[i];
		bT[i] = b_1[i];
	}
	ks1_read(ks1_res32);

	//---w4---
	ks1_write_int32(aT, bT);	//(a_1, b_1);

	//----- previous final step(3 of 3)
	for (i = 0; i < w_len; i++) { //w1,w7,w4已经计算过了
		res_full[sub_len + i] += w6[i];
		res_full[sub_len2 + i] += w5[i];
		res_full[sub_len4 + i] += w3[i];
		res_full[sub_len5 + i] += w2[i];
	}
	//---------------reduction-------
	for (i = 0; i < 255; i++) { //res[0 ~ (2n-2)-n=n-2=256-2=254]，剩下res[n-1]
		acc[i] = (acc[i] + res_full[i] - res_full[i + n]) & result_mask;
	}
	acc[n - 1] = (acc[n - 1] + res_full[n - 1]) & result_mask;
	//----- previous final step(3 of 3)

	toom4_saber_10bits_reduce_w(w3, ks1_res32, toom4_precision_mask);	//->w3
	for (i = 0; i < sub_len; i++) {
		//--- w2=a2*b2=8A3+4A2+2A1+A0 ---
		a2[i] = a[sub_len3 + i] << 1;	//2A3
		b2[i] = b[sub_len3 + i] << 1;
		a2[i] += a[sub_len2 + i];	//2A3+A2
		b2[i] += b[sub_len2 + i];
		a2[i] <<= 1;	//4A3+2A2
		b2[i] <<= 1;
		a2[i] += a[sub_len + i];	//4A3+2A2+A1
		b2[i] += b[sub_len + i];
		a2[i] <<= 1;	//8A3+2A2+2A1
		b2[i] <<= 1;
		a2[i] += a[i];	//8A3+2A2+2A1+A0
		b2[i] += b[i];

		a2[i] = a2[i] & toom4_precision_mask;
		b2[i] = b2[i] & toom4_precision_mask;

		aT[i] = ahalf[i];
		bT[i] = bhalf[i];
	}

	memset(res_full, 0, sizeof(res_full[0]) * res_full_len);
	for (i = 0; i < w_len; i++) {
		res_full[sub_len6 + i] += w1[i];
	}
	ks1_read(ks1_res32);

	//---w5---
	ks1_write_int32(aT, bT);	//(ahalf, bhalf);
	toom4_saber_10bits_reduce_w(w4, ks1_res32, toom4_precision_mask);	//->w4
	//[step12,14]
	for (i = 0; i < w_len; i++) {
		//step12: w4=(w4-w3)/2
		w4[i] = w4[i] - w3[i];
		w4[i] = w4[i] >> 1;
		//step14: w3=w3+w4
		w3[i] = w3[i] + w4[i]; //w3 <- w3+w4
	}
	for (i = 0; i < sub_len; i++) {
		aT[i] = a2[i];
		bT[i] = b2[i];
	}
	ks1_read(ks1_res32);

	//---w2---
	ks1_write_int32(aT, bT);	//(a2, b2);
	toom4_saber_10bits_reduce_w(w5, ks1_res32, toom4_precision_mask); //->w5
	//[step13.1,copy_w5]
	for (i = 0; i < w_len; i++) {
		w5_copy[i] = w5[i];	//copy_w5
		//step13.1:w5=w5-w1;[w5=w5-w1-64*w7]
		w5[i] = w5[i] - w1[i];
	}
	for (i = 0; i < sub_len; i++) {
		aT[i] = a0[i];
		bT[i] = b0[i];
	}
	ks1_read(ks1_res32);

	//--w7---
	ks1_write_int32(aT, bT);	//(a0, b0);
	toom4_saber_10bits_reduce_w(w2, ks1_res32, toom4_precision_mask);		//->w2
	//[step10,16,17.1]
	for (i = 0; i < w_len; i++) {
		//step10: w2=w2+w5
		w2[i] = w2[i] + w5_copy[i];	//w2 <- w2+w5
		//step16: w2=w2-65*w3
		w2[i] = w2[i] - (w3[i] << 6) - w3[i];
		//step17.1:w3=w3-w1; [w3=w3-w7-w1]
		w3[i] = w3[i] - w1[i];
	}
	for (i = 0; i < sub_len; i++) {
		aT[i] = a_half[i];
		bT[i] = b_half[i];
	}
	ks1_read(ks1_res32);

	//---w6---
	ks1_write_int32(aT, bT);	//(a_half, b_half);
	toom4_saber_10bits_reduce_w(w7, ks1_res32, toom4_precision_mask);		//->w7
	for (i = 0; i < w_len; i++) {
		res_full[i] += w7[i];
	}
	//[step13.2,17.2,18,21,23]  w2_copy
	for (i = 0; i < w_len; i++) {
		//step13.2"w5=w5-64*w7;[w5=w5-w1-64*w7]
		w5[i] = w5[i] - (w7[i] << 6);
		//step17.2:w3=w3-w7; [w3=w3-w7-w1]
		w3[i] = w3[i] - w7[i];
		//step18: w2=w2+45*w3
		w2[i] = w2[i] + w3[i] * int45; //w2 <- w2+45*w3
		w2_copy[i] = w2[i]; //w2_copy
		//step21: w2=(w2+16*w4)/18
		w2[i] = w2[i] + (w4[i] << 4); //w2 <- w2+16*w4
		w2[i] = w2[i] * inv9;
		w2[i] = w2[i] >> 1;
		//step23: w4=-(w4+w2/2)
		w4[i] = -(w4[i] + w2[i]); //->w4_final
	}
	for (i = 0; i < w_len; i++) {
		res_full[sub_len3 + i] += w4[i];
	}
	ks1_read(ks1_res32);

	//====== va2*vb2 ===========

	a = &va2[0];
	b = &vb2[0];

	ainf = &a[sub_len3];	//A(inf),B(inf)->w1, inf:infinite
	binf = &b[sub_len3];

	//---w1---
	ks1_write_uint16(ainf, binf);

	//----- previous final step(1 of 3)
	toom4_saber_10bits_reduce_w(w6, ks1_res32, toom4_precision_mask);		//->w6
	//----- previous final step(1 of 3)

	a0 = &a[0];	//A0,B0->w7
	b0 = &b[0];

	for (i = 0; i < sub_len; i++) {
		//--- w3,w4=a(±1)*b(±1)=(A2+A0)±(A3+A1) ---
		at1 = a[sub_len2 + i] + a[i];	//A2 + A0
		bt1 = b[sub_len2 + i] + b[i];
		at2 = a[sub_len3 + i] + a[sub_len + i];	//A3 + A1
		bt2 = b[sub_len3 + i] + b[sub_len + i];
		a1[i] = at1 + at2; //(A2 + A0) + (A3 + A1) = A0 + A1 + A2 + A3
		b1[i] = bt1 + bt2;
		a_1[i] = at1 - at2;	//(A2 + A0) - (A3 + A1) = A0 - A1 + A2 - A3
		b_1[i] = bt1 - bt2;

		a1[i] = a1[i] & toom4_precision_mask;
		b1[i] = b1[i] & toom4_precision_mask;
		a_1[i] = a_1[i] & toom4_precision_mask;
		b_1[i] = b_1[i] & toom4_precision_mask;

		aT[i] = a1[i];
		bT[i] = b1[i];
	}
	//toom4_reduce_a_b_unsigned(a1, b1);
	//toom4_reduce_a_b_signed(a_1, b_1);
	ks1_read(ks1_res32);

	//---w3---
	ks1_write_int32(aT, bT);	//a1, b1;
	//----- previous final step(2 of 3) begin
	for (i = 0; i < w_len; i++) {
		//step11: w6=w6-w5
		w6[i] = w6[i] - w5_copy[i];	// w6 <- w6-w5

		//step15: w5=2*w5+w6
		w5[i] = w6[i] + (w5[i] << 1); //w5 <- 2*w5+w6

		//step19: w5=(w5-8*w3)/24
		w5[i] = w5[i] - (w3[i] << 3); //w5 <- w5-8*w3
		w5[i] = w5[i] * inv3;
		w5[i] = w5[i] >> 3;

		//step20: w6=w6+w2
		w6[i] = w2_copy[i] + w6[i]; //w6 <- w6+w2

		//step22. w3 <- w3-w5
		w3[i] = w3[i] - w5[i];

		//step24:  w6=(30*w2-w6)/60
		w6[i] = w2[i] * int30 - w6[i]; //w6 <- 30*w2-w6
		w6[i] = w6[i] * inv15; //*inv15 / 4
		w6[i] = w6[i] >> 2;

		//step25: w2=w2-w6
		w2[i] = w2[i] - w6[i];
	}
	//----- previous final step(2 of 3) end

	toom4_saber_10bits_reduce_w(w1, ks1_res32, toom4_precision_mask);	//->w1
	//poly_print_uint16(w1, "w1", 256, 0);
	for (i = 0; i < sub_len; i++) {
		//---w5,w6=a±half*b±half=A(±1/2)*B(±1/2)=(8A0+2A2)±(4A1+A3)---
		at1 = a0[i] << 2;	//4*A0
		bt1 = b0[i] << 2;
		at1 = at1 + a[sub_len2 + i];	//4A0+A2
		bt1 = bt1 + b[sub_len2 + i];
		at1 = at1 << 1;	//8A0 + 2A2
		bt1 = bt1 << 1;
		at2 = a[sub_len + i];	//A1
		bt2 = b[sub_len + i];
		at2 = at2 << 2;	//4A1
		bt2 = bt2 << 2;
		at2 = at2 + a[sub_len3 + i];	//4A1 + A3
		bt2 = bt2 + b[sub_len3 + i];
		ahalf[i] = at1 + at2;	//(8A0 + 2A2) + (4A1 + A3) = 8A0 + 4A1 + 2A2 + A3
		bhalf[i] = bt1 + bt2;
		a_half[i] = at1 - at2;	//(8A0 + 2A2) - (4A1 + A3) = 8A0 - 4A1 - 2A2 + A3
		b_half[i] = bt1 - bt2;

		ahalf[i] = ahalf[i] & toom4_precision_mask;
		bhalf[i] = bhalf[i] & toom4_precision_mask;
		a_half[i] = a_half[i] & toom4_precision_mask;
		b_half[i] = b_half[i] & toom4_precision_mask;

		aT[i] = a_1[i];
		bT[i] = b_1[i];
	}
	ks1_read(ks1_res32);

	//---w4---
	ks1_write_int32(aT, bT);	//(a_1, b_1);

	//----- previous final step(3 of 3) begin
	for (i = 0; i < w_len; i++) { //w1,w7,w4已经计算过了
		res_full[sub_len + i] += w6[i];
		res_full[sub_len2 + i] += w5[i];
		res_full[sub_len4 + i] += w3[i];
		res_full[sub_len5 + i] += w2[i];
	}
	//---------------reduction-------
	for (i = 0; i < 255; i++) { //res[0 ~ (2n-2)-n=n-2=256-2=254]，剩下res[n-1]
		acc[i] = (acc[i] + res_full[i] - res_full[i + n]) & result_mask;
	}
	acc[n - 1] = (acc[n - 1] + res_full[n - 1]) & result_mask;
	//----- previous final step(3 of 3) end

	toom4_saber_10bits_reduce_w(w3, ks1_res32, toom4_precision_mask);	//->w3
	for (i = 0; i < sub_len; i++) {
		//--- w2=a2*b2=8A3+4A2+2A1+A0 ---
		a2[i] = a[sub_len3 + i] << 1;	//2A3
		b2[i] = b[sub_len3 + i] << 1;
		a2[i] += a[sub_len2 + i];	//2A3+A2
		b2[i] += b[sub_len2 + i];
		a2[i] <<= 1;	//4A3+2A2
		b2[i] <<= 1;
		a2[i] += a[sub_len + i];	//4A3+2A2+A1
		b2[i] += b[sub_len + i];
		a2[i] <<= 1;	//8A3+2A2+2A1
		b2[i] <<= 1;
		a2[i] += a[i];	//8A3+2A2+2A1+A0
		b2[i] += b[i];

		a2[i] = a2[i] & toom4_precision_mask;
		b2[i] = b2[i] & toom4_precision_mask;

		aT[i] = ahalf[i];
		bT[i] = bhalf[i];
	}
	memset(res_full, 0, sizeof(res_full[0]) * res_full_len);
	for (i = 0; i < w_len; i++) {
		res_full[sub_len6 + i] += w1[i];
	}
	ks1_read(ks1_res32);

	//---w5---
	ks1_write_int32(aT, bT);	//(ahalf, bhalf);
	toom4_saber_10bits_reduce_w(w4, ks1_res32, toom4_precision_mask);	//->w4
	//[step12,14]
	for (i = 0; i < w_len; i++) {
		//step12: w4=(w4-w3)/2
		w4[i] = w4[i] - w3[i];
		w4[i] = w4[i] >> 1;
		//step14: w3=w3+w4
		w3[i] = w3[i] + w4[i]; //w3 <- w3+w4
	}
	for (i = 0; i < sub_len; i++) {
		aT[i] = a2[i];
		bT[i] = b2[i];
	}
	ks1_read(ks1_res32);

	//---w2---
	ks1_write_int32(aT, bT);	//(a2, b2);
	toom4_saber_10bits_reduce_w(w5, ks1_res32, toom4_precision_mask); //->w5
	//[step13.1,copy_w5]
	for (i = 0; i < w_len; i++) {
		w5_copy[i] = w5[i];	//copy_w5
		//step13.1:w5=w5-w1;[w5=w5-w1-64*w7]
		w5[i] = w5[i] - w1[i];
	}
	for (i = 0; i < sub_len; i++) {
		aT[i] = a0[i];
		bT[i] = b0[i];
	}
	ks1_read(ks1_res32);

	//--w7---
	ks1_write_int32(aT, bT);	//(a0, b0);
	toom4_saber_10bits_reduce_w(w2, ks1_res32, toom4_precision_mask);		//->w2
	//[step10,16,17.1]
	for (i = 0; i < w_len; i++) {
		//step10: w2=w2+w5
		w2[i] = w2[i] + w5_copy[i];	//w2 <- w2+w5
		//step16: w2=w2-65*w3
		w2[i] = w2[i] - (w3[i] << 6) - w3[i];
		//step17.1:w3=w3-w1; [w3=w3-w7-w1]
		w3[i] = w3[i] - w1[i];
	}
	for (i = 0; i < sub_len; i++) {
		aT[i] = a_half[i];
		bT[i] = b_half[i];
	}
	ks1_read(ks1_res32);

	//---w6---
	ks1_write_int32(aT, bT);	//(a_half, b_half);
	toom4_saber_10bits_reduce_w(w7, ks1_res32, toom4_precision_mask);		//->w7
	for (i = 0; i < w_len; i++) {
		res_full[i] += w7[i];
	}
	//[step13.2,17.2,18,21,23]  w2_copy
	for (i = 0; i < w_len; i++) {
		//step13.2"w5=w5-64*w7;[w5=w5-w1-64*w7]
		w5[i] = w5[i] - (w7[i] << 6);
		//step17.2:w3=w3-w7; [w3=w3-w7-w1]
		w3[i] = w3[i] - w7[i];
		//step18: w2=w2+45*w3
		w2[i] = w2[i] + w3[i] * int45; //w2 <- w2+45*w3
		w2_copy[i] = w2[i]; //w2_copy
		//step21: w2=(w2+16*w4)/18
		w2[i] = w2[i] + (w4[i] << 4); //w2 <- w2+16*w4
		w2[i] = w2[i] * inv9;
		w2[i] = w2[i] >> 1;
		//step23: w4=-(w4+w2/2)
		w4[i] = -(w4[i] + w2[i]); //->w4_final
	}
	for (i = 0; i < w_len; i++) {
		res_full[sub_len3 + i] += w4[i];
	}
	ks1_read(ks1_res32);

	//====以下是不能并行的=====
	toom4_saber_10bits_reduce_w(w6, ks1_res32, toom4_precision_mask);		//->w6

	for (i = 0; i < w_len; i++) {
		//step11: w6=w6-w5
		w6[i] = w6[i] - w5_copy[i];	// w6 <- w6-w5

		//step15: w5=2*w5+w6
		w5[i] = w6[i] + (w5[i] << 1); //w5 <- 2*w5+w6

		//step19: w5=(w5-8*w3)/24
		w5[i] = w5[i] - (w3[i] << 3); //w5 <- w5-8*w3
		w5[i] = w5[i] * inv3;
		w5[i] = w5[i] >> 3;

		//step20: w6=w6+w2
		w6[i] = w2_copy[i] + w6[i]; //w6 <- w6+w2

		//step22. w3 <- w3-w5
		w3[i] = w3[i] - w5[i];

		//step24:  w6=(30*w2-w6)/60
		w6[i] = w2[i] * int30 - w6[i]; //w6 <- 30*w2-w6
		w6[i] = w6[i] * inv15; //*inv15 / 4
		w6[i] = w6[i] >> 2;

		//step25: w2=w2-w6
		w2[i] = w2[i] - w6[i];
	}

	for (i = 0; i < w_len; i++) { //w1,w7,w4已经计算过了
		res_full[sub_len + i] += w6[i];
		res_full[sub_len2 + i] += w5[i];
		res_full[sub_len4 + i] += w3[i];
		res_full[sub_len5 + i] += w2[i];
	}
	//---------------reduction-------
	for (i = 0; i < 255; i++) { //res[0 ~ (2n-2)-n=n-2=256-2=254]，剩下res[n-1]
		//result[i] = (res_full[i] - res_full[i + n]) & result_mask;
		acc[i] = (acc[i] + res_full[i] - res_full[i + n]) & result_mask;
	}
	//result[n - 1] = res_full[n - 1] & result_mask;
	acc[n - 1] = (acc[n - 1] + res_full[n - 1]) & result_mask;

}

void toom4_saber_10bits_parallel_vector(uint16_t* va0, uint16_t* va1, uint16_t* va2,
		uint16_t* vb0, uint16_t* vb1, uint16_t* vb2, uint16_t* acc) {
	//printf("---toom4_saber_10bits_parallel---\n");

	REG_WRITE_NOCHECK(RSA_MULT_MODE_REG, 15);
	const uint16_t toom4_precision_mask = SABER_P * 8 - 1;
	const uint16_t result_mask = SABER_P - 1;

	const int16_t n = 256;
	int16_t i;
	uint16_t sub_len = 64;	//64,n / 4
	uint16_t sub_len2 = 128;
	uint16_t sub_len3 = 192;
	uint16_t sub_len4 = 256;
	uint16_t sub_len5 = 320;
	uint16_t sub_len6 = 384;
	uint16_t w_len = 127;	//2 * sub_len - 1;

	uint16_t at1, at2;  //temp for w
	uint16_t bt1, bt2;

	int16_t inv3 = 43691;
	int16_t inv9 = 36409;
	int16_t inv15 = 61167;

	int16_t int45 = 45;
	int16_t int30 = 30;

	//====== va0*vb0 ===========

	uint16_t* a = &va0[0];
	uint16_t* b = &vb0[0];

	uint16_t* ainf = &a[sub_len3];	//A(inf),B(inf)->w1, inf:infinite
	uint16_t* binf = &b[sub_len3];

	//---w1---
	ks1_write_uint16(ainf, binf);

	memset(acc, 0, sizeof(acc[0]) * SABER_N);

	uint16_t a2[sub_len], b2[sub_len];	//A(2)*B(2)->w2
	uint16_t a1[sub_len], b1[sub_len];	//A(1)*B(1)->w3
	uint16_t a_1[sub_len], b_1[sub_len];	//A(-1)*B(-1)->w4
	uint16_t ahalf[sub_len], bhalf[sub_len];	//A(1/2)*B(1/2)->w5
	uint16_t a_half[sub_len], b_half[sub_len];	//A(-1/2)*B(-1/2)->w6
	uint16_t* a0 = &a[0];	//A0,B0->w7
	uint16_t* b0 = &b[0];

	uint16_t w1[w_len], w2[w_len], w3[w_len], w4[w_len], w5[w_len], w6[w_len], w7[w_len];			//w_len=128

	uint16_t w5_copy[w_len];	//备份w5给step10,11
	uint16_t w2_copy[w_len];	//在step21之前备份w5给step20

	uint32_t ks1_res32[127];	//存放从寄存器读取的ks1结果
	int32_t aT[sub_len]; //32-bit alignment
	int32_t bT[sub_len];

	for (i = 0; i < sub_len; i++) {
		//--- w3,w4=a(±1)*b(±1)=(A2+A0)±(A3+A1) ---
		at1 = a[sub_len2 + i] + a[i];	//A2 + A0
		bt1 = b[sub_len2 + i] + b[i];
		at2 = a[sub_len3 + i] + a[sub_len + i];	//A3 + A1
		bt2 = b[sub_len3 + i] + b[sub_len + i];
		a1[i] = at1 + at2; //(A2 + A0) + (A3 + A1) = A0 + A1 + A2 + A3
		b1[i] = bt1 + bt2;
		a_1[i] = at1 - at2;	//(A2 + A0) - (A3 + A1) = A0 - A1 + A2 - A3
		b_1[i] = bt1 - bt2;

		a1[i] = a1[i] & toom4_precision_mask;
		b1[i] = b1[i] & toom4_precision_mask;
		a_1[i] = a_1[i] & toom4_precision_mask;
		b_1[i] = b_1[i] & toom4_precision_mask;

		aT[i] = a1[i];
		bT[i] = b1[i];
	}
	ks1_read(ks1_res32);

	//---w3---
	ks1_write_int32(aT, bT);	//a1, b1;
	toom4_saber_10bits_reduce_w(w1, ks1_res32, toom4_precision_mask);	//->w1
	//poly_print_uint16(w1, "w1", 256, 0);
	for (i = 0; i < sub_len; i++) {
		//---w5,w6=a±half*b±half=A(±1/2)*B(±1/2)=(8A0+2A2)±(4A1+A3)---
		at1 = a0[i] << 2;	//4*A0
		bt1 = b0[i] << 2;
		at1 = at1 + a[sub_len2 + i];	//4A0+A2
		bt1 = bt1 + b[sub_len2 + i];
		at1 = at1 << 1;	//8A0 + 2A2
		bt1 = bt1 << 1;
		at2 = a[sub_len + i];	//A1
		bt2 = b[sub_len + i];
		at2 = at2 << 2;	//4A1
		bt2 = bt2 << 2;
		at2 = at2 + a[sub_len3 + i];	//4A1 + A3
		bt2 = bt2 + b[sub_len3 + i];
		ahalf[i] = at1 + at2;	//(8A0 + 2A2) + (4A1 + A3) = 8A0 + 4A1 + 2A2 + A3
		bhalf[i] = bt1 + bt2;
		a_half[i] = at1 - at2;	//(8A0 + 2A2) - (4A1 + A3) = 8A0 - 4A1 - 2A2 + A3
		b_half[i] = bt1 - bt2;

		ahalf[i] = ahalf[i] & toom4_precision_mask;
		bhalf[i] = bhalf[i] & toom4_precision_mask;
		a_half[i] = a_half[i] & toom4_precision_mask;
		b_half[i] = b_half[i] & toom4_precision_mask;

		aT[i] = a_1[i];
		bT[i] = b_1[i];
	}
	ks1_read(ks1_res32);

	//---w4---
	ks1_write_int32(aT, bT);	//(a_1, b_1);
	toom4_saber_10bits_reduce_w(w3, ks1_res32, toom4_precision_mask);	//->w3
	for (i = 0; i < sub_len; i++) {
		//--- w2=a2*b2=8A3+4A2+2A1+A0 ---
		a2[i] = a[sub_len3 + i] << 1;	//2A3
		b2[i] = b[sub_len3 + i] << 1;
		a2[i] += a[sub_len2 + i];	//2A3+A2
		b2[i] += b[sub_len2 + i];
		a2[i] <<= 1;	//4A3+2A2
		b2[i] <<= 1;
		a2[i] += a[sub_len + i];	//4A3+2A2+A1
		b2[i] += b[sub_len + i];
		a2[i] <<= 1;	//8A3+2A2+2A1
		b2[i] <<= 1;
		a2[i] += a[i];	//8A3+2A2+2A1+A0
		b2[i] += b[i];

		a2[i] = a2[i] & toom4_precision_mask;
		b2[i] = b2[i] & toom4_precision_mask;

		aT[i] = ahalf[i];
		bT[i] = bhalf[i];
	}
	const int16_t res_full_len = 2 * n - 1;
	uint16_t res_full[res_full_len];
	memset(res_full, 0, sizeof(res_full[0]) * res_full_len);
	for (i = 0; i < w_len; i++) {
		res_full[sub_len6 + i] += w1[i];
	}
	ks1_read(ks1_res32);

	//---w5---
	ks1_write_int32(aT, bT);	//(ahalf, bhalf);
	toom4_saber_10bits_reduce_w(w4, ks1_res32, toom4_precision_mask);	//->w4
	//[step12,14]
	for (i = 0; i < w_len; i++) {
		w4[i] = w4[i] - w3[i];//step12: w4=(w4-w3)/2
		w4[i] = w4[i] >> 1;
		w3[i] = w3[i] + w4[i]; //step14: w3=w3+w4
	}
	for (i = 0; i < sub_len; i++) {
		aT[i] = a2[i];
		bT[i] = b2[i];
	}
	ks1_read(ks1_res32);

	//---w2---
	ks1_write_int32(aT, bT);	//(a2, b2);
	toom4_saber_10bits_reduce_w(w5, ks1_res32, toom4_precision_mask); //->w5
	//[step13.1,copy_w5]
	for (i = 0; i < w_len; i++) {
		w5_copy[i] = w5[i];	//copy_w5
		w5[i] = w5[i] - w1[i];//step13.1:w5=w5-w1;[w5=w5-w1-64*w7]
	}
	for (i = 0; i < sub_len; i++) {
		aT[i] = a0[i];
		bT[i] = b0[i];
	}
	ks1_read(ks1_res32);

	//--w7---
	ks1_write_int32(aT, bT);	//(a0, b0);
	toom4_saber_10bits_reduce_w(w2, ks1_res32, toom4_precision_mask);		//->w2
	//[step10,16,17.1]
	for (i = 0; i < w_len; i++) {
		w2[i] = w2[i] + w5_copy[i];	//step10: w2=w2+w5
		w2[i] = w2[i] - (w3[i] << 6) - w3[i]; //step16: w2=w2-65*w3
		w3[i] = w3[i] - w1[i]; //step17.1:w3=w3-w1; [w3=w3-w7-w1]
	}
	for (i = 0; i < sub_len; i++) {
		aT[i] = a_half[i];
		bT[i] = b_half[i];
	}
	ks1_read(ks1_res32);

	//---w6---
	ks1_write_int32(aT, bT);	//(a_half, b_half);
	toom4_saber_10bits_reduce_w(w7, ks1_res32, toom4_precision_mask);		//->w7
	for (i = 0; i < w_len; i++) {
		res_full[i] += w7[i];
	}
	//[step13.2,17.2,18,21,23]  w2_copy
	for (i = 0; i < w_len; i++) {
		w5[i] = w5[i] - (w7[i] << 6);		//step13.2:w5=w5-64*w7;[w5=w5-w1-64*w7]
		w3[i] = w3[i] - w7[i];		//step17.2:w3=w3-w7; [w3=w3-w7-w1]
		w2[i] = w2[i] + w3[i] * int45; //step18: w2=w2+45*w3
		w2_copy[i] = w2[i]; //w2_copy
		w2[i] = w2[i] + (w4[i] << 4); //step21: w2=(w2+16*w4)/18
		w2[i] = w2[i] * inv9;
		w2[i] = w2[i] >> 1;
		w4[i] = -(w4[i] + w2[i]); ////step23: w4=-(w4+w2/2)->w4_final
	}
	for (i = 0; i < w_len; i++) {
		res_full[sub_len3 + i] += w4[i];
	}

	ks1_read(ks1_res32);

	//====== va1*vb1 ===========

	a = &va1[0];
	b = &vb1[0];

	ainf = &a[sub_len3];	//A(inf),B(inf)->w1, inf:infinite
	binf = &b[sub_len3];


	//---w1---
	ks1_write_uint16(ainf, binf);
	//----- previous final step(1 of 5)
	toom4_saber_10bits_reduce_w(w6, ks1_res32, toom4_precision_mask);		//->w6
	for (i = 0; i < w_len; i++) {
		w6[i] = w6[i] - w5_copy[i];	//step11: w6=w6-w5
		w5[i] = w6[i] + (w5[i] << 1); //step15: w5=2*w5+w6
		w5[i] = w5[i] - (w3[i] << 3); //step19: w5=(w5-8*w3)/24
		w5[i] = w5[i] * inv3;
		w5[i] = w5[i] >> 3;
	}
	//----- previous final step(1 of 5)

	a0 = &a[0];	//A0,B0->w7
	b0 = &b[0];

	for (i = 0; i < sub_len; i++) {
		//--- w3,w4=a(±1)*b(±1)=(A2+A0)±(A3+A1) ---
		at1 = a[sub_len2 + i] + a[i];	//A2 + A0
		bt1 = b[sub_len2 + i] + b[i];
		at2 = a[sub_len3 + i] + a[sub_len + i];	//A3 + A1
		bt2 = b[sub_len3 + i] + b[sub_len + i];
		a1[i] = at1 + at2; //(A2 + A0) + (A3 + A1) = A0 + A1 + A2 + A3
		b1[i] = bt1 + bt2;
		a_1[i] = at1 - at2;	//(A2 + A0) - (A3 + A1) = A0 - A1 + A2 - A3
		b_1[i] = bt1 - bt2;

		a1[i] = a1[i] & toom4_precision_mask;
		b1[i] = b1[i] & toom4_precision_mask;
		a_1[i] = a_1[i] & toom4_precision_mask;
		b_1[i] = b_1[i] & toom4_precision_mask;

		aT[i] = a1[i];
		bT[i] = b1[i];
	}
	ks1_read(ks1_res32);

	//---w3---
	ks1_write_int32(aT, bT);	//a1, b1;

	//----- previous final step(2 of 5) begin
	for (i = 0; i < w_len; i++) {
		w6[i] = w2_copy[i] + w6[i]; //step20: w6=w6+w2
		w3[i] = w3[i] - w5[i]; //step22. w3 <- w3-w5
	}
	//----- previous final step(2 of 5) end

	toom4_saber_10bits_reduce_w(w1, ks1_res32, toom4_precision_mask);	//->w1
	for (i = 0; i < sub_len; i++) {
		//---w5,w6=a±half*b±half=A(±1/2)*B(±1/2)=(8A0+2A2)±(4A1+A3)---
		at1 = a0[i] << 2;	//4*A0
		bt1 = b0[i] << 2;
		at1 = at1 + a[sub_len2 + i];	//4A0+A2
		bt1 = bt1 + b[sub_len2 + i];
		at1 = at1 << 1;	//8A0 + 2A2
		bt1 = bt1 << 1;
		at2 = a[sub_len + i];	//A1
		bt2 = b[sub_len + i];
		at2 = at2 << 2;	//4A1
		bt2 = bt2 << 2;
		at2 = at2 + a[sub_len3 + i];	//4A1 + A3
		bt2 = bt2 + b[sub_len3 + i];
		ahalf[i] = at1 + at2;	//(8A0 + 2A2) + (4A1 + A3) = 8A0 + 4A1 + 2A2 + A3
		bhalf[i] = bt1 + bt2;
		a_half[i] = at1 - at2;	//(8A0 + 2A2) - (4A1 + A3) = 8A0 - 4A1 - 2A2 + A3
		b_half[i] = bt1 - bt2;

		ahalf[i] = ahalf[i] & toom4_precision_mask;
		bhalf[i] = bhalf[i] & toom4_precision_mask;
		a_half[i] = a_half[i] & toom4_precision_mask;
		b_half[i] = b_half[i] & toom4_precision_mask;

		aT[i] = a_1[i];
		bT[i] = b_1[i];
	}
	ks1_read(ks1_res32);

	//---w4---
	ks1_write_int32(aT, bT);	//(a_1, b_1);
	//----- previous final step(3 of 5) end
	for (i = 0; i < w_len; i++) {
		w6[i] = w2[i] * int30 - w6[i]; //step24:  w6=(30*w2-w6)/60
		w6[i] = w6[i] * inv15; //*inv15 / 4
		w6[i] = w6[i] >> 2;
		w2[i] = w2[i] - w6[i]; //step25: w2=w2-w6
	}
	//----- previous final step(3 of 5) end

	//----- previous final step(4 of 5) begin
	for (i = 0; i < w_len; i++) { //w1,w7,w4已经计算过了
		res_full[sub_len + i] += w6[i];
		res_full[sub_len2 + i] += w5[i];
		res_full[sub_len4 + i] += w3[i];
		res_full[sub_len5 + i] += w2[i];
	}
	//----- previous final step(4 of 5) end

	toom4_saber_10bits_reduce_w(w3, ks1_res32, toom4_precision_mask);	//->w3
	for (i = 0; i < sub_len; i++) {
		aT[i] = ahalf[i];
		bT[i] = bhalf[i];
	}
	ks1_read(ks1_res32);

	//---w5---
	ks1_write_int32(aT, bT);	//(ahalf, bhalf);

	//----- previous final step(5 of 5) begin
	for (i = 0; i < 255; i++) { //res[0 ~ (2n-2)-n=n-2=256-2=254]，剩下res[n-1]
		acc[i] = (acc[i] + res_full[i] - res_full[i + n]); // & result_mask;
	}
	acc[n - 1] = (acc[n - 1] + res_full[n - 1]); // & result_mask;
	//----- previous final step(5 of 5) end

	toom4_saber_10bits_reduce_w(w4, ks1_res32, toom4_precision_mask);	//->w4

	for (i = 0; i < sub_len; i++) {
		//--- w2=a2*b2=8A3+4A2+2A1+A0 ---
		a2[i] = a[sub_len3 + i] << 1;	//2A3
		b2[i] = b[sub_len3 + i] << 1;
		a2[i] += a[sub_len2 + i];	//2A3+A2
		b2[i] += b[sub_len2 + i];
		a2[i] <<= 1;	//4A3+2A2
		b2[i] <<= 1;
		a2[i] += a[sub_len + i];	//4A3+2A2+A1
		b2[i] += b[sub_len + i];
		a2[i] <<= 1;	//8A3+2A2+2A1
		b2[i] <<= 1;
		a2[i] += a[i];	//8A3+2A2+2A1+A0
		b2[i] += b[i];

		a2[i] = a2[i] & toom4_precision_mask;
		b2[i] = b2[i] & toom4_precision_mask;
	}
	for (i = 0; i < sub_len; i++) {
		aT[i] = a2[i];
		bT[i] = b2[i];
	}
	ks1_read(ks1_res32);

	//---w2---
	ks1_write_int32(aT, bT);	//(a2, b2);
	toom4_saber_10bits_reduce_w(w5, ks1_res32, toom4_precision_mask); //->w5
	memset(res_full, 0, sizeof(res_full[0]) * res_full_len);
	for (i = 0; i < w_len; i++) {
		res_full[sub_len6 + i] += w1[i];
	}
	//[step12,14]
	for (i = 0; i < w_len; i++) {
		w4[i] = w4[i] - w3[i]; //step12: w4=(w4-w3)/2
		w4[i] = w4[i] >> 1;
		w3[i] = w3[i] + w4[i]; //step14: w3=w3+w4
	}
	//[step13.1,copy_w5]
	for (i = 0; i < w_len; i++) {
		w5_copy[i] = w5[i];	//copy_w5
		w5[i] = w5[i] - w1[i];	//step13.1:w5=w5-w1;[w5=w5-w1-64*w7]
	}
	for (i = 0; i < sub_len; i++) {
		aT[i] = a0[i];
		bT[i] = b0[i];
	}
	ks1_read(ks1_res32);

	//--w7---
	ks1_write_int32(aT, bT);	//(a0, b0);
	toom4_saber_10bits_reduce_w(w2, ks1_res32, toom4_precision_mask);		//->w2
	//[step10,16,17.1]
	for (i = 0; i < w_len; i++) {
		w2[i] = w2[i] + w5_copy[i];	//step10: w2=w2+w5
		w2[i] = w2[i] - (w3[i] << 6) - w3[i];	//step16: w2=w2-65*w3
		w3[i] = w3[i] - w1[i];	//step17.1:w3=w3-w1; [w3=w3-w7-w1]
	}
	for (i = 0; i < sub_len; i++) {
		aT[i] = a_half[i];
		bT[i] = b_half[i];
	}
	ks1_read(ks1_res32);

	//---w6---
	ks1_write_int32(aT, bT);	//(a_half, b_half);
	toom4_saber_10bits_reduce_w(w7, ks1_res32, toom4_precision_mask);		//->w7
	for (i = 0; i < w_len; i++) {
		res_full[i] += w7[i];
	}
	//[step13.2,17.2,18,21,23]  w2_copy
	for (i = 0; i < w_len; i++) {
		w5[i] = w5[i] - (w7[i] << 6);		//step13.2:w5=w5-64*w7;[w5=w5-w1-64*w7]
		w3[i] = w3[i] - w7[i];		//step17.2:w3=w3-w7; [w3=w3-w7-w1]
		w2[i] = w2[i] + w3[i] * int45;		//step18: w2=w2+45*w3
		w2_copy[i] = w2[i]; //w2_copy
		w2[i] = w2[i] + (w4[i] << 4); //step21: w2=(w2+16*w4)/18
		w2[i] = w2[i] * inv9;
		w2[i] = w2[i] >> 1;
		w4[i] = -(w4[i] + w2[i]); //step23: w4=-(w4+w2/2) w4 ok
	}
	for (i = 0; i < w_len; i++) {
		res_full[sub_len3 + i] += w4[i];
	}
	ks1_read(ks1_res32);

	//====== va2*vb2 ===========

	a = &va2[0];
	b = &vb2[0];

	ainf = &a[sub_len3];	//A(inf),B(inf)->w1, inf:infinite
	binf = &b[sub_len3];

	//---w1---
	ks1_write_uint16(ainf, binf);
	//----- previous final step(1 of 5)
	toom4_saber_10bits_reduce_w(w6, ks1_res32, toom4_precision_mask);		//->w6
	for (i = 0; i < w_len; i++) {
		w6[i] = w6[i] - w5_copy[i];	//step11: w6=w6-w5
		w5[i] = w6[i] + (w5[i] << 1); //step15: w5=2*w5+w6
		w5[i] = w5[i] - (w3[i] << 3); //step19: w5=(w5-8*w3)/24
		w5[i] = w5[i] * inv3;
		w5[i] = w5[i] >> 3;
	}
	//----- previous final step(1 of 5)

	a0 = &a[0];	//A0,B0->w7
	b0 = &b[0];

	for (i = 0; i < sub_len; i++) {
		//--- w3,w4=a(±1)*b(±1)=(A2+A0)±(A3+A1) ---
		at1 = a[sub_len2 + i] + a[i];	//A2 + A0
		bt1 = b[sub_len2 + i] + b[i];
		at2 = a[sub_len3 + i] + a[sub_len + i];	//A3 + A1
		bt2 = b[sub_len3 + i] + b[sub_len + i];
		a1[i] = at1 + at2; //(A2 + A0) + (A3 + A1) = A0 + A1 + A2 + A3
		b1[i] = bt1 + bt2;
		a_1[i] = at1 - at2;	//(A2 + A0) - (A3 + A1) = A0 - A1 + A2 - A3
		b_1[i] = bt1 - bt2;

		a1[i] = a1[i] & toom4_precision_mask;
		b1[i] = b1[i] & toom4_precision_mask;
		a_1[i] = a_1[i] & toom4_precision_mask;
		b_1[i] = b_1[i] & toom4_precision_mask;

		aT[i] = a1[i];
		bT[i] = b1[i];
	}
	ks1_read(ks1_res32);

	//---w3---
	ks1_write_int32(aT, bT);	//a1, b1;

	//----- previous final step(2 of 5) begin
	for (i = 0; i < w_len; i++) {
		w6[i] = w2_copy[i] + w6[i]; //step20: w6=w6+w2
		w3[i] = w3[i] - w5[i]; //step22. w3 <- w3-w5
	}
	//----- previous final step(2 of 5) end

	toom4_saber_10bits_reduce_w(w1, ks1_res32, toom4_precision_mask);	//->w1
	for (i = 0; i < sub_len; i++) {
		//---w5,w6=a±half*b±half=A(±1/2)*B(±1/2)=(8A0+2A2)±(4A1+A3)---
		at1 = a0[i] << 2;	//4*A0
		bt1 = b0[i] << 2;
		at1 = at1 + a[sub_len2 + i];	//4A0+A2
		bt1 = bt1 + b[sub_len2 + i];
		at1 = at1 << 1;	//8A0 + 2A2
		bt1 = bt1 << 1;
		at2 = a[sub_len + i];	//A1
		bt2 = b[sub_len + i];
		at2 = at2 << 2;	//4A1
		bt2 = bt2 << 2;
		at2 = at2 + a[sub_len3 + i];	//4A1 + A3
		bt2 = bt2 + b[sub_len3 + i];
		ahalf[i] = at1 + at2;	//(8A0 + 2A2) + (4A1 + A3) = 8A0 + 4A1 + 2A2 + A3
		bhalf[i] = bt1 + bt2;
		a_half[i] = at1 - at2;	//(8A0 + 2A2) - (4A1 + A3) = 8A0 - 4A1 - 2A2 + A3
		b_half[i] = bt1 - bt2;

		ahalf[i] = ahalf[i] & toom4_precision_mask;
		bhalf[i] = bhalf[i] & toom4_precision_mask;
		a_half[i] = a_half[i] & toom4_precision_mask;
		b_half[i] = b_half[i] & toom4_precision_mask;

		aT[i] = a_1[i];
		bT[i] = b_1[i];
	}
	ks1_read(ks1_res32);

	//---w4---
	ks1_write_int32(aT, bT);	//(a_1, b_1);
	//----- previous final step(3 of 5) end
	for (i = 0; i < w_len; i++) {
		w6[i] = w2[i] * int30 - w6[i]; //step24:  w6=(30*w2-w6)/60
		w6[i] = w6[i] * inv15; //*inv15 / 4
		w6[i] = w6[i] >> 2;
		w2[i] = w2[i] - w6[i]; //step25: w2=w2-w6
	}
	//----- previous final step(3 of 5) end

	//----- previous final step(4 of 5) begin
	for (i = 0; i < w_len; i++) { //w1,w7,w4已经计算过了
		res_full[sub_len + i] += w6[i];
		res_full[sub_len2 + i] += w5[i];
		res_full[sub_len4 + i] += w3[i];
		res_full[sub_len5 + i] += w2[i];
	}
	//----- previous final step(4 of 5) end

	toom4_saber_10bits_reduce_w(w3, ks1_res32, toom4_precision_mask);	//->w3
	for (i = 0; i < sub_len; i++) {
		aT[i] = ahalf[i];
		bT[i] = bhalf[i];
	}

	ks1_read(ks1_res32);

	//---w5---
	ks1_write_int32(aT, bT);	//(ahalf, bhalf);

	//----- previous final step(4 of 4) begin
	for (i = 0; i < 255; i++) { //res[0 ~ (2n-2)-n=n-2=256-2=254]，剩下res[n-1]
		acc[i] = (acc[i] + res_full[i] - res_full[i + n]); // & result_mask;
	}
	acc[n - 1] = (acc[n - 1] + res_full[n - 1]); // & result_mask;
	//----- previous final step(4 of 4) end

	toom4_saber_10bits_reduce_w(w4, ks1_res32, toom4_precision_mask);	//->w4

	for (i = 0; i < sub_len; i++) {
		//--- w2=a2*b2=8A3+4A2+2A1+A0 ---
		a2[i] = a[sub_len3 + i] << 1;	//2A3
		b2[i] = b[sub_len3 + i] << 1;
		a2[i] += a[sub_len2 + i];	//2A3+A2
		b2[i] += b[sub_len2 + i];
		a2[i] <<= 1;	//4A3+2A2
		b2[i] <<= 1;
		a2[i] += a[sub_len + i];	//4A3+2A2+A1
		b2[i] += b[sub_len + i];
		a2[i] <<= 1;	//8A3+2A2+2A1
		b2[i] <<= 1;
		a2[i] += a[i];	//8A3+2A2+2A1+A0
		b2[i] += b[i];

		a2[i] = a2[i] & toom4_precision_mask;
		b2[i] = b2[i] & toom4_precision_mask;
	}
	for (i = 0; i < sub_len; i++) {
		aT[i] = a2[i];
		bT[i] = b2[i];
	}
	ks1_read(ks1_res32);

	//---w2---
	ks1_write_int32(aT, bT);	//(a2, b2);
	toom4_saber_10bits_reduce_w(w5, ks1_res32, toom4_precision_mask); //->w5
	memset(res_full, 0, sizeof(res_full[0]) * res_full_len);
	for (i = 0; i < w_len; i++) {
		res_full[sub_len6 + i] += w1[i];
	}
	//[step12,14]
	for (i = 0; i < w_len; i++) {
		w4[i] = w4[i] - w3[i]; //step12: w4=(w4-w3)/2
		w4[i] = w4[i] >> 1;
		w3[i] = w3[i] + w4[i]; //step14: w3=w3+w4
	}
	//[step13.1,copy_w5]
	for (i = 0; i < w_len; i++) {
		w5_copy[i] = w5[i];	//copy_w5
		w5[i] = w5[i] - w1[i];	//step13.1:w5=w5-w1;[w5=w5-w1-64*w7]
	}
	for (i = 0; i < sub_len; i++) {
		aT[i] = a0[i];
		bT[i] = b0[i];
	}
	ks1_read(ks1_res32);

	//--w7---
	ks1_write_int32(aT, bT);	//(a0, b0);
	toom4_saber_10bits_reduce_w(w2, ks1_res32, toom4_precision_mask);		//->w2
	//[step10,16,17.1]
	for (i = 0; i < w_len; i++) {
		w2[i] = w2[i] + w5_copy[i];	//step10: w2=w2+w5
		w2[i] = w2[i] - (w3[i] << 6) - w3[i];	//step16: w2=w2-65*w3
		w3[i] = w3[i] - w1[i];	//step17.1:w3=w3-w1; [w3=w3-w7-w1]
	}
	for (i = 0; i < sub_len; i++) {
		aT[i] = a_half[i];
		bT[i] = b_half[i];
	}
	ks1_read(ks1_res32);

	//---w6---
	ks1_write_int32(aT, bT);	//(a_half, b_half);
	toom4_saber_10bits_reduce_w(w7, ks1_res32, toom4_precision_mask);		//->w7
	for (i = 0; i < w_len; i++) {
		res_full[i] += w7[i];
	}
	//[step13.2,17.2,18,21,23]  w2_copy
	for (i = 0; i < w_len; i++) {
		w5[i] = w5[i] - (w7[i] << 6);		//step13.2:w5=w5-64*w7;[w5=w5-w1-64*w7]
		w3[i] = w3[i] - w7[i];		//step17.2:w3=w3-w7; [w3=w3-w7-w1]
		w2[i] = w2[i] + w3[i] * int45; //step18: w2=w2+45*w3
		w2_copy[i] = w2[i]; //w2_copy
		w2[i] = w2[i] + (w4[i] << 4); //step21: w2=(w2+16*w4)/18
		w2[i] = w2[i] * inv9;
		w2[i] = w2[i] >> 1;
		w4[i] = -(w4[i] + w2[i]); //step23: w4=-(w4+w2/2)
	}
	for (i = 0; i < w_len; i++) {
		res_full[sub_len3 + i] += w4[i];
	}
	ks1_read(ks1_res32);

	//====以下是不能并行的=====
	toom4_saber_10bits_reduce_w(w6, ks1_res32, toom4_precision_mask);

	for (i = 0; i < w_len; i++) {
		w6[i] = w6[i] - w5_copy[i];	//step11: w6=w6-w5
		w5[i] = w6[i] + (w5[i] << 1); //step15: w5=2*w5+w6
		w5[i] = w5[i] - (w3[i] << 3); //step19: w5=(w5-8*w3)/24
		w5[i] = w5[i] * inv3;
		w5[i] = w5[i] >> 3;
		w6[i] = w2_copy[i] + w6[i]; //step20: w6=w6+w2
		w3[i] = w3[i] - w5[i]; //step22. w3 <- w3-w5
		w6[i] = w2[i] * int30 - w6[i]; //step24:  w6=(30*w2-w6)/60
		w6[i] = w6[i] * inv15; //*inv15 / 4
		w6[i] = w6[i] >> 2;
		w2[i] = w2[i] - w6[i]; //step25: w2=w2-w6
	}

	for (i = 0; i < w_len; i++) { //w1,w7,w4已经计算过了
		res_full[sub_len + i] += w6[i];
		res_full[sub_len2 + i] += w5[i];
		res_full[sub_len4 + i] += w3[i];
		res_full[sub_len5 + i] += w2[i];
	}
	//---------------reduction-------
	for (i = 0; i < 255; i++) { //res[0 ~ (2n-2)-n=n-2=256-2=254]，剩下res[n-1]
		//result[i] = (res_full[i] - res_full[i + n]) & result_mask;
		acc[i] = (acc[i] + res_full[i] - res_full[i + n]) & result_mask;
	}
	//result[n - 1] = res_full[n - 1] & result_mask;
	acc[n - 1] = (acc[n - 1] + res_full[n - 1]) & result_mask;
}

void toom4_saber_10bits_parallel_vector_speedtest(uint16_t* va0, uint16_t* va1, uint16_t* va2,
		uint16_t* vb0, uint16_t* vb1, uint16_t* vb2, uint16_t* acc) {
	//printf("---toom4_saber_10bits_parallel---\n");
	speedtest_print_title("toom4_saber_10bits_parallel_vector_speedtest", 1);
	speedtest_reset_startcpucycles();
	speedtest_print_splitline();
	speedtest_print_splitline();

	REG_WRITE_NOCHECK(RSA_MULT_MODE_REG, 15);
	const uint16_t toom4_precision_mask = SABER_P * 8 - 1;
	const uint16_t result_mask = SABER_P - 1;

	const int16_t n = 256;
	int16_t i;
	uint16_t sub_len = 64;	//64,n / 4
	uint16_t sub_len2 = 128;
	uint16_t sub_len3 = 192;
	uint16_t sub_len4 = 256;
	uint16_t sub_len5 = 320;
	uint16_t sub_len6 = 384;
	uint16_t w_len = 127;	//2 * sub_len - 1;

	uint16_t at1, at2;  //temp for w
	uint16_t bt1, bt2;

	int16_t inv3 = 43691;
	int16_t inv9 = 36409;
	int16_t inv15 = 61167;

	int16_t int45 = 45;
	int16_t int30 = 30;

	//====== va0*vb0 ===========

	uint16_t* a = &va0[0];
	uint16_t* b = &vb0[0];

	uint16_t* ainf = &a[sub_len3];	//A(inf),B(inf)->w1, inf:infinite
	uint16_t* binf = &b[sub_len3];

	speedtest_print_cpucycles("init_all_pointers");
	speedtest_print_splitline();

	//---w1---
	ks1_write_uint16(ainf, binf);
	speedtest_print_cpucycles("w1_ainf_binf_write");

	memset(acc, 0, sizeof(acc[0]) * SABER_N);

	uint16_t a2[sub_len], b2[sub_len];	//A(2)*B(2)->w2
	uint16_t a1[sub_len], b1[sub_len];	//A(1)*B(1)->w3
	uint16_t a_1[sub_len], b_1[sub_len];	//A(-1)*B(-1)->w4
	uint16_t ahalf[sub_len], bhalf[sub_len];	//A(1/2)*B(1/2)->w5
	uint16_t a_half[sub_len], b_half[sub_len];	//A(-1/2)*B(-1/2)->w6
	uint16_t* a0 = &a[0];	//A0,B0->w7
	uint16_t* b0 = &b[0];

	uint16_t w1[w_len], w2[w_len], w3[w_len], w4[w_len], w5[w_len], w6[w_len], w7[w_len];			//w_len=128

	uint16_t w5_copy[w_len];	//备份w5给step10,11
	uint16_t w2_copy[w_len];	//在step21之前备份w5给step20

	uint32_t ks1_res32[127];	//存放从寄存器读取的ks1结果
	int32_t aT[sub_len]; //32-bit alignment
	int32_t bT[sub_len];

	for (i = 0; i < sub_len; i++) {
		//--- w3,w4=a(±1)*b(±1)=(A2+A0)±(A3+A1) ---
		at1 = a[sub_len2 + i] + a[i];	//A2 + A0
		bt1 = b[sub_len2 + i] + b[i];
		at2 = a[sub_len3 + i] + a[sub_len + i];	//A3 + A1
		bt2 = b[sub_len3 + i] + b[sub_len + i];
		a1[i] = at1 + at2; //(A2 + A0) + (A3 + A1) = A0 + A1 + A2 + A3
		b1[i] = bt1 + bt2;
		a_1[i] = at1 - at2;	//(A2 + A0) - (A3 + A1) = A0 - A1 + A2 - A3
		b_1[i] = bt1 - bt2;

		a1[i] = a1[i] & toom4_precision_mask;
		b1[i] = b1[i] & toom4_precision_mask;
		a_1[i] = a_1[i] & toom4_precision_mask;
		b_1[i] = b_1[i] & toom4_precision_mask;

		aT[i] = a1[i];
		bT[i] = b1[i];
	}
	speedtest_print_cpucycles("w1_ainf_binf_wait_cpu_used");
	ks1_read(ks1_res32);
	speedtest_print_cpucycles("w1_ainf_binf_read");
	speedtest_print_splitline();

	//---w3---
	ks1_write_int32(aT, bT);	//a1, b1;
	speedtest_print_cpucycles("w3_a1_b1_write");
	toom4_saber_10bits_reduce_w(w1, ks1_res32, toom4_precision_mask);	//->w1
	//poly_print_uint16(w1, "w1", 256, 0);
	for (i = 0; i < sub_len; i++) {
		//---w5,w6=a±half*b±half=A(±1/2)*B(±1/2)=(8A0+2A2)±(4A1+A3)---
		at1 = a0[i] << 2;	//4*A0
		bt1 = b0[i] << 2;
		at1 = at1 + a[sub_len2 + i];	//4A0+A2
		bt1 = bt1 + b[sub_len2 + i];
		at1 = at1 << 1;	//8A0 + 2A2
		bt1 = bt1 << 1;
		at2 = a[sub_len + i];	//A1
		bt2 = b[sub_len + i];
		at2 = at2 << 2;	//4A1
		bt2 = bt2 << 2;
		at2 = at2 + a[sub_len3 + i];	//4A1 + A3
		bt2 = bt2 + b[sub_len3 + i];
		ahalf[i] = at1 + at2;	//(8A0 + 2A2) + (4A1 + A3) = 8A0 + 4A1 + 2A2 + A3
		bhalf[i] = bt1 + bt2;
		a_half[i] = at1 - at2;	//(8A0 + 2A2) - (4A1 + A3) = 8A0 - 4A1 - 2A2 + A3
		b_half[i] = bt1 - bt2;

		ahalf[i] = ahalf[i] & toom4_precision_mask;
		bhalf[i] = bhalf[i] & toom4_precision_mask;
		a_half[i] = a_half[i] & toom4_precision_mask;
		b_half[i] = b_half[i] & toom4_precision_mask;

		aT[i] = a_1[i];
		bT[i] = b_1[i];
	}
	speedtest_print_cpucycles("w3_a1_b1_wait_cpu_used");
	ks1_read(ks1_res32);
	speedtest_print_cpucycles("w3_a1_b1_read");
	speedtest_print_splitline();

	//---w4---
	ks1_write_int32(aT, bT);	//(a_1, b_1);
	speedtest_print_cpucycles("w4_a_1_b_1_write");
	toom4_saber_10bits_reduce_w(w3, ks1_res32, toom4_precision_mask);	//->w3
	for (i = 0; i < sub_len; i++) {
		//--- w2=a2*b2=8A3+4A2+2A1+A0 ---
		a2[i] = a[sub_len3 + i] << 1;	//2A3
		b2[i] = b[sub_len3 + i] << 1;
		a2[i] += a[sub_len2 + i];	//2A3+A2
		b2[i] += b[sub_len2 + i];
		a2[i] <<= 1;	//4A3+2A2
		b2[i] <<= 1;
		a2[i] += a[sub_len + i];	//4A3+2A2+A1
		b2[i] += b[sub_len + i];
		a2[i] <<= 1;	//8A3+2A2+2A1
		b2[i] <<= 1;
		a2[i] += a[i];	//8A3+2A2+2A1+A0
		b2[i] += b[i];

		a2[i] = a2[i] & toom4_precision_mask;
		b2[i] = b2[i] & toom4_precision_mask;

		aT[i] = ahalf[i];
		bT[i] = bhalf[i];
	}
	const int16_t res_full_len = 2 * n - 1;
	uint16_t res_full[res_full_len];
	memset(res_full, 0, sizeof(res_full[0]) * res_full_len);
	for (i = 0; i < w_len; i++) {
		res_full[sub_len6 + i] += w1[i];
	}
	speedtest_print_cpucycles("w4_a_1_b_1_wait_cpu_used");
	ks1_read(ks1_res32);
	speedtest_print_cpucycles("w4_a_1_b_1_read");
	speedtest_print_splitline();

	//---w5---
	ks1_write_int32(aT, bT);	//(ahalf, bhalf);
	speedtest_print_cpucycles("w5_ahalf_bhalf_write");
	toom4_saber_10bits_reduce_w(w4, ks1_res32, toom4_precision_mask);	//->w4
	//[step12,14]
	for (i = 0; i < w_len; i++) {
		w4[i] = w4[i] - w3[i];//step12: w4=(w4-w3)/2
		w4[i] = w4[i] >> 1;
		w3[i] = w3[i] + w4[i]; //step14: w3=w3+w4
	}
	for (i = 0; i < sub_len; i++) {
		aT[i] = a2[i];
		bT[i] = b2[i];
	}
	speedtest_print_cpucycles("w5_ahalf_bhalf_wait_cpu_used");
	ks1_read(ks1_res32);
	speedtest_print_cpucycles("w5_ahalf_bhalf_read");
	speedtest_print_splitline();

	//---w2---
	ks1_write_int32(aT, bT);	//(a2, b2);
	speedtest_print_cpucycles("w2_a2_b2_write");
	toom4_saber_10bits_reduce_w(w5, ks1_res32, toom4_precision_mask); //->w5
	//[step13.1,copy_w5]
	for (i = 0; i < w_len; i++) {
		w5_copy[i] = w5[i];	//copy_w5
		w5[i] = w5[i] - w1[i];//step13.1:w5=w5-w1;[w5=w5-w1-64*w7]
	}
	for (i = 0; i < sub_len; i++) {
		aT[i] = a0[i];
		bT[i] = b0[i];
	}
	speedtest_print_cpucycles("w2_a2_b2_wait_cpu_used");
	ks1_read(ks1_res32);
	speedtest_print_cpucycles("w2_a2_b2_read");
	speedtest_print_splitline();

	//--w7---
	ks1_write_int32(aT, bT);	//(a0, b0);
	speedtest_print_cpucycles("w7_a0_b0_write");
	toom4_saber_10bits_reduce_w(w2, ks1_res32, toom4_precision_mask);		//->w2
	//[step10,16,17.1]
	for (i = 0; i < w_len; i++) {
		w2[i] = w2[i] + w5_copy[i];	//step10: w2=w2+w5
		w2[i] = w2[i] - (w3[i] << 6) - w3[i]; //step16: w2=w2-65*w3
		w3[i] = w3[i] - w1[i]; //step17.1:w3=w3-w1; [w3=w3-w7-w1]
	}
	for (i = 0; i < sub_len; i++) {
		aT[i] = a_half[i];
		bT[i] = b_half[i];
	}
	speedtest_print_cpucycles("w7_a0_b0_wait_cpu_used");
	ks1_read(ks1_res32);
	speedtest_print_cpucycles("w7_a0_b0_read");
	speedtest_print_splitline();

	//---w6---
	ks1_write_int32(aT, bT);	//(a_half, b_half);
	speedtest_print_cpucycles("w6_a_half_b_half_write");
	toom4_saber_10bits_reduce_w(w7, ks1_res32, toom4_precision_mask);		//->w7
	for (i = 0; i < w_len; i++) {
		res_full[i] += w7[i];
	}
	//[step13.2,17.2,18,21,23]  w2_copy
	for (i = 0; i < w_len; i++) {
		w5[i] = w5[i] - (w7[i] << 6);		//step13.2:w5=w5-64*w7;[w5=w5-w1-64*w7]
		w3[i] = w3[i] - w7[i];		//step17.2:w3=w3-w7; [w3=w3-w7-w1]
		w2[i] = w2[i] + w3[i] * int45; //step18: w2=w2+45*w3
		w2_copy[i] = w2[i]; //w2_copy
		w2[i] = w2[i] + (w4[i] << 4); //step21: w2=(w2+16*w4)/18
		w2[i] = w2[i] * inv9;
		w2[i] = w2[i] >> 1;
		w4[i] = -(w4[i] + w2[i]); ////step23: w4=-(w4+w2/2)->w4_final
	}
	for (i = 0; i < w_len; i++) {
		res_full[sub_len3 + i] += w4[i];
	}

	speedtest_print_cpucycles("w6_a_half_b_half_wait_cpu_used");
	ks1_read(ks1_res32);
	speedtest_print_cpucycles("w6_a_half_b_half_read");
	speedtest_print_splitline();

	//====== va1*vb1 ===========

	a = &va1[0];
	b = &vb1[0];

	ainf = &a[sub_len3];	//A(inf),B(inf)->w1, inf:infinite
	binf = &b[sub_len3];

	speedtest_print_cpucycles("init_all_pointers");
	speedtest_print_splitline();

	//---w1---
	ks1_write_uint16(ainf, binf);
	speedtest_print_cpucycles("w1_ainf_binf_write");
	//----- previous final step(1 of 5)
	toom4_saber_10bits_reduce_w(w6, ks1_res32, toom4_precision_mask);		//->w6
	for (i = 0; i < w_len; i++) {
		w6[i] = w6[i] - w5_copy[i];	//step11: w6=w6-w5
		w5[i] = w6[i] + (w5[i] << 1); //step15: w5=2*w5+w6
		w5[i] = w5[i] - (w3[i] << 3); //step19: w5=(w5-8*w3)/24
		w5[i] = w5[i] * inv3;
		w5[i] = w5[i] >> 3;
	}
	//----- previous final step(1 of 5)

	a0 = &a[0];	//A0,B0->w7
	b0 = &b[0];

	for (i = 0; i < sub_len; i++) {
		//--- w3,w4=a(±1)*b(±1)=(A2+A0)±(A3+A1) ---
		at1 = a[sub_len2 + i] + a[i];	//A2 + A0
		bt1 = b[sub_len2 + i] + b[i];
		at2 = a[sub_len3 + i] + a[sub_len + i];	//A3 + A1
		bt2 = b[sub_len3 + i] + b[sub_len + i];
		a1[i] = at1 + at2; //(A2 + A0) + (A3 + A1) = A0 + A1 + A2 + A3
		b1[i] = bt1 + bt2;
		a_1[i] = at1 - at2;	//(A2 + A0) - (A3 + A1) = A0 - A1 + A2 - A3
		b_1[i] = bt1 - bt2;

		a1[i] = a1[i] & toom4_precision_mask;
		b1[i] = b1[i] & toom4_precision_mask;
		a_1[i] = a_1[i] & toom4_precision_mask;
		b_1[i] = b_1[i] & toom4_precision_mask;

		aT[i] = a1[i];
		bT[i] = b1[i];
	}
	speedtest_print_cpucycles("w1_ainf_binf_wait_cpu_used");
	ks1_read(ks1_res32);
	speedtest_print_cpucycles("w1_ainf_binf_read");
	speedtest_print_splitline();

	//---w3---
	ks1_write_int32(aT, bT);	//a1, b1;
	speedtest_print_cpucycles("w3_a1_b1_write");

	//----- previous final step(2 of 5) begin
	for (i = 0; i < w_len; i++) {
		w6[i] = w2_copy[i] + w6[i]; //step20: w6=w6+w2
		w3[i] = w3[i] - w5[i]; //step22. w3 <- w3-w5
	}
	//----- previous final step(2 of 5) end

	toom4_saber_10bits_reduce_w(w1, ks1_res32, toom4_precision_mask);	//->w1
	for (i = 0; i < sub_len; i++) {
		//---w5,w6=a±half*b±half=A(±1/2)*B(±1/2)=(8A0+2A2)±(4A1+A3)---
		at1 = a0[i] << 2;	//4*A0
		bt1 = b0[i] << 2;
		at1 = at1 + a[sub_len2 + i];	//4A0+A2
		bt1 = bt1 + b[sub_len2 + i];
		at1 = at1 << 1;	//8A0 + 2A2
		bt1 = bt1 << 1;
		at2 = a[sub_len + i];	//A1
		bt2 = b[sub_len + i];
		at2 = at2 << 2;	//4A1
		bt2 = bt2 << 2;
		at2 = at2 + a[sub_len3 + i];	//4A1 + A3
		bt2 = bt2 + b[sub_len3 + i];
		ahalf[i] = at1 + at2;	//(8A0 + 2A2) + (4A1 + A3) = 8A0 + 4A1 + 2A2 + A3
		bhalf[i] = bt1 + bt2;
		a_half[i] = at1 - at2;	//(8A0 + 2A2) - (4A1 + A3) = 8A0 - 4A1 - 2A2 + A3
		b_half[i] = bt1 - bt2;

		ahalf[i] = ahalf[i] & toom4_precision_mask;
		bhalf[i] = bhalf[i] & toom4_precision_mask;
		a_half[i] = a_half[i] & toom4_precision_mask;
		b_half[i] = b_half[i] & toom4_precision_mask;

		aT[i] = a_1[i];
		bT[i] = b_1[i];
	}
	speedtest_print_cpucycles("w3_a1_b1_wait_cpu_used");
	ks1_read(ks1_res32);
	speedtest_print_cpucycles("w3_a1_b1_read");
	speedtest_print_splitline();

	//---w4---
	ks1_write_int32(aT, bT);	//(a_1, b_1);
	speedtest_print_cpucycles("w4_a_1_b_1_write");
	//----- previous final step(3 of 5) end
	for (i = 0; i < w_len; i++) {
		w6[i] = w2[i] * int30 - w6[i]; //step24:  w6=(30*w2-w6)/60
		w6[i] = w6[i] * inv15; //*inv15 / 4
		w6[i] = w6[i] >> 2;
		w2[i] = w2[i] - w6[i]; //step25: w2=w2-w6
	}
	//----- previous final step(3 of 5) end

	//----- previous final step(4 of 5) begin
	for (i = 0; i < w_len; i++) { //w1,w7,w4已经计算过了
		res_full[sub_len + i] += w6[i];
		res_full[sub_len2 + i] += w5[i];
		res_full[sub_len4 + i] += w3[i];
		res_full[sub_len5 + i] += w2[i];
	}
	//----- previous final step(4 of 5) end

	toom4_saber_10bits_reduce_w(w3, ks1_res32, toom4_precision_mask);	//->w3
	for (i = 0; i < sub_len; i++) {
		aT[i] = ahalf[i];
		bT[i] = bhalf[i];
	}
	speedtest_print_cpucycles("w4_a_1_b_1_wait_cpu_used");
	ks1_read(ks1_res32);
	speedtest_print_cpucycles("w4_a_1_b_1_read");
	speedtest_print_splitline();

	//---w5---
	ks1_write_int32(aT, bT);	//(ahalf, bhalf);
	speedtest_print_cpucycles("w5_ahalf_bhalf_write");

	//----- previous final step(5 of 5) begin
	for (i = 0; i < 255; i++) { //res[0 ~ (2n-2)-n=n-2=256-2=254]，剩下res[n-1]
		acc[i] = (acc[i] + res_full[i] - res_full[i + n]); // & result_mask;
	}
	acc[n - 1] = (acc[n - 1] + res_full[n - 1]); // & result_mask;
	//----- previous final step(5 of 5) end

	toom4_saber_10bits_reduce_w(w4, ks1_res32, toom4_precision_mask);	//->w4

	for (i = 0; i < sub_len; i++) {
		//--- w2=a2*b2=8A3+4A2+2A1+A0 ---
		a2[i] = a[sub_len3 + i] << 1;	//2A3
		b2[i] = b[sub_len3 + i] << 1;
		a2[i] += a[sub_len2 + i];	//2A3+A2
		b2[i] += b[sub_len2 + i];
		a2[i] <<= 1;	//4A3+2A2
		b2[i] <<= 1;
		a2[i] += a[sub_len + i];	//4A3+2A2+A1
		b2[i] += b[sub_len + i];
		a2[i] <<= 1;	//8A3+2A2+2A1
		b2[i] <<= 1;
		a2[i] += a[i];	//8A3+2A2+2A1+A0
		b2[i] += b[i];

		a2[i] = a2[i] & toom4_precision_mask;
		b2[i] = b2[i] & toom4_precision_mask;
	}
	for (i = 0; i < sub_len; i++) {
		aT[i] = a2[i];
		bT[i] = b2[i];
	}
	speedtest_print_cpucycles("w5_ahalf_bhalf_wait_cpu_used");
	ks1_read(ks1_res32);
	speedtest_print_cpucycles("w5_ahalf_bhalf_read");
	speedtest_print_splitline();

	//---w2---
	ks1_write_int32(aT, bT);	//(a2, b2);
	speedtest_print_cpucycles("w2_a2_b2_write");
	toom4_saber_10bits_reduce_w(w5, ks1_res32, toom4_precision_mask); //->w5
	memset(res_full, 0, sizeof(res_full[0]) * res_full_len);
	for (i = 0; i < w_len; i++) {
		res_full[sub_len6 + i] += w1[i];
	}
	//[step12,14]
	for (i = 0; i < w_len; i++) {
		w4[i] = w4[i] - w3[i]; //step12: w4=(w4-w3)/2
		w4[i] = w4[i] >> 1;
		w3[i] = w3[i] + w4[i]; //step14: w3=w3+w4
	}
	//[step13.1,copy_w5]
	for (i = 0; i < w_len; i++) {
		w5_copy[i] = w5[i];	//copy_w5
		w5[i] = w5[i] - w1[i];	//step13.1:w5=w5-w1;[w5=w5-w1-64*w7]
	}
	for (i = 0; i < sub_len; i++) {
		aT[i] = a0[i];
		bT[i] = b0[i];
	}
	speedtest_print_cpucycles("w2_a2_b2_wait_cpu_used");
	ks1_read(ks1_res32);
	speedtest_print_cpucycles("w2_a2_b2_read");
	speedtest_print_splitline();

	//--w7---
	ks1_write_int32(aT, bT);	//(a0, b0);
	speedtest_print_cpucycles("w7_a0_b0_write");
	toom4_saber_10bits_reduce_w(w2, ks1_res32, toom4_precision_mask);		//->w2
	//[step10,16,17.1]
	for (i = 0; i < w_len; i++) {
		w2[i] = w2[i] + w5_copy[i];	//step10: w2=w2+w5
		w2[i] = w2[i] - (w3[i] << 6) - w3[i];	//step16: w2=w2-65*w3
		w3[i] = w3[i] - w1[i];	//step17.1:w3=w3-w1; [w3=w3-w7-w1]
	}
	for (i = 0; i < sub_len; i++) {
		aT[i] = a_half[i];
		bT[i] = b_half[i];
	}
	speedtest_print_cpucycles("w7_a0_b0_wait_cpu_used");
	ks1_read(ks1_res32);
	speedtest_print_cpucycles("w7_a0_b0_read");
	speedtest_print_splitline();

	//---w6---
	ks1_write_int32(aT, bT);	//(a_half, b_half);
	speedtest_print_cpucycles("w6_a_half_b_half_write");
	toom4_saber_10bits_reduce_w(w7, ks1_res32, toom4_precision_mask);		//->w7
	for (i = 0; i < w_len; i++) {
		res_full[i] += w7[i];
	}
	//[step13.2,17.2,18,21,23]  w2_copy
	for (i = 0; i < w_len; i++) {
		w5[i] = w5[i] - (w7[i] << 6);		//step13.2:w5=w5-64*w7;[w5=w5-w1-64*w7]
		w3[i] = w3[i] - w7[i];		//step17.2:w3=w3-w7; [w3=w3-w7-w1]
		w2[i] = w2[i] + w3[i] * int45;		//step18: w2=w2+45*w3
		w2_copy[i] = w2[i]; //w2_copy
		w2[i] = w2[i] + (w4[i] << 4); //step21: w2=(w2+16*w4)/18
		w2[i] = w2[i] * inv9;
		w2[i] = w2[i] >> 1;
		w4[i] = -(w4[i] + w2[i]); //step23: w4=-(w4+w2/2) w4 ok
	}
	for (i = 0; i < w_len; i++) {
		res_full[sub_len3 + i] += w4[i];
	}
	speedtest_print_cpucycles("w6_a_half_b_half_wait_cpu_used");
	ks1_read(ks1_res32);
	speedtest_print_cpucycles("w6_a_half_b_half_read");
	speedtest_print_splitline();

	//====== va2*vb2 ===========

	a = &va2[0];
	b = &vb2[0];

	ainf = &a[sub_len3];	//A(inf),B(inf)->w1, inf:infinite
	binf = &b[sub_len3];

	speedtest_print_cpucycles("init_all_pointers");
	speedtest_print_splitline();

	//---w1---
	ks1_write_uint16(ainf, binf);
	speedtest_print_cpucycles("w1_ainf_binf_write");
	//----- previous final step(1 of 5)
	toom4_saber_10bits_reduce_w(w6, ks1_res32, toom4_precision_mask);		//->w6
	for (i = 0; i < w_len; i++) {
		w6[i] = w6[i] - w5_copy[i];	//step11: w6=w6-w5
		w5[i] = w6[i] + (w5[i] << 1); //step15: w5=2*w5+w6
		w5[i] = w5[i] - (w3[i] << 3); //step19: w5=(w5-8*w3)/24
		w5[i] = w5[i] * inv3;
		w5[i] = w5[i] >> 3;
	}
	//----- previous final step(1 of 5)

	a0 = &a[0];	//A0,B0->w7
	b0 = &b[0];

	for (i = 0; i < sub_len; i++) {
		//--- w3,w4=a(±1)*b(±1)=(A2+A0)±(A3+A1) ---
		at1 = a[sub_len2 + i] + a[i];	//A2 + A0
		bt1 = b[sub_len2 + i] + b[i];
		at2 = a[sub_len3 + i] + a[sub_len + i];	//A3 + A1
		bt2 = b[sub_len3 + i] + b[sub_len + i];
		a1[i] = at1 + at2; //(A2 + A0) + (A3 + A1) = A0 + A1 + A2 + A3
		b1[i] = bt1 + bt2;
		a_1[i] = at1 - at2;	//(A2 + A0) - (A3 + A1) = A0 - A1 + A2 - A3
		b_1[i] = bt1 - bt2;

		a1[i] = a1[i] & toom4_precision_mask;
		b1[i] = b1[i] & toom4_precision_mask;
		a_1[i] = a_1[i] & toom4_precision_mask;
		b_1[i] = b_1[i] & toom4_precision_mask;

		aT[i] = a1[i];
		bT[i] = b1[i];
	}
	speedtest_print_cpucycles("w1_ainf_binf_wait_cpu_used");
	ks1_read(ks1_res32);
	speedtest_print_cpucycles("w1_ainf_binf_read");
	speedtest_print_splitline();

	//---w3---
	ks1_write_int32(aT, bT);	//a1, b1;
	speedtest_print_cpucycles("w3_a1_b1_write");

	//----- previous final step(2 of 5) begin
	for (i = 0; i < w_len; i++) {
		w6[i] = w2_copy[i] + w6[i]; //step20: w6=w6+w2
		w3[i] = w3[i] - w5[i]; //step22. w3 <- w3-w5
	}
	//----- previous final step(2 of 5) end

	toom4_saber_10bits_reduce_w(w1, ks1_res32, toom4_precision_mask);	//->w1
	for (i = 0; i < sub_len; i++) {
		//---w5,w6=a±half*b±half=A(±1/2)*B(±1/2)=(8A0+2A2)±(4A1+A3)---
		at1 = a0[i] << 2;	//4*A0
		bt1 = b0[i] << 2;
		at1 = at1 + a[sub_len2 + i];	//4A0+A2
		bt1 = bt1 + b[sub_len2 + i];
		at1 = at1 << 1;	//8A0 + 2A2
		bt1 = bt1 << 1;
		at2 = a[sub_len + i];	//A1
		bt2 = b[sub_len + i];
		at2 = at2 << 2;	//4A1
		bt2 = bt2 << 2;
		at2 = at2 + a[sub_len3 + i];	//4A1 + A3
		bt2 = bt2 + b[sub_len3 + i];
		ahalf[i] = at1 + at2;	//(8A0 + 2A2) + (4A1 + A3) = 8A0 + 4A1 + 2A2 + A3
		bhalf[i] = bt1 + bt2;
		a_half[i] = at1 - at2;	//(8A0 + 2A2) - (4A1 + A3) = 8A0 - 4A1 - 2A2 + A3
		b_half[i] = bt1 - bt2;

		ahalf[i] = ahalf[i] & toom4_precision_mask;
		bhalf[i] = bhalf[i] & toom4_precision_mask;
		a_half[i] = a_half[i] & toom4_precision_mask;
		b_half[i] = b_half[i] & toom4_precision_mask;

		aT[i] = a_1[i];
		bT[i] = b_1[i];
	}
	speedtest_print_cpucycles("w3_a1_b1_wait_cpu_used");
	ks1_read(ks1_res32);
	speedtest_print_cpucycles("w3_a1_b1_read");
	speedtest_print_splitline();

	//---w4---
	ks1_write_int32(aT, bT);	//(a_1, b_1);
	speedtest_print_cpucycles("w4_a_1_b_1_write");
	//----- previous final step(3 of 5) end
	for (i = 0; i < w_len; i++) {
		w6[i] = w2[i] * int30 - w6[i]; //step24:  w6=(30*w2-w6)/60
		w6[i] = w6[i] * inv15; //*inv15 / 4
		w6[i] = w6[i] >> 2;
		w2[i] = w2[i] - w6[i]; //step25: w2=w2-w6
	}
	//----- previous final step(3 of 5) end

	//----- previous final step(4 of 5) begin
	for (i = 0; i < w_len; i++) { //w1,w7,w4已经计算过了
		res_full[sub_len + i] += w6[i];
		res_full[sub_len2 + i] += w5[i];
		res_full[sub_len4 + i] += w3[i];
		res_full[sub_len5 + i] += w2[i];
	}
	//----- previous final step(4 of 5) end

	toom4_saber_10bits_reduce_w(w3, ks1_res32, toom4_precision_mask);	//->w3
	for (i = 0; i < sub_len; i++) {
		aT[i] = ahalf[i];
		bT[i] = bhalf[i];
	}

	speedtest_print_cpucycles("w4_a_1_b_1_wait_cpu_used");
	ks1_read(ks1_res32);
	speedtest_print_cpucycles("w4_a_1_b_1_read");
	speedtest_print_splitline();

	//---w5---
	ks1_write_int32(aT, bT);	//(ahalf, bhalf);
	speedtest_print_cpucycles("w5_ahalf_bhalf_write");

	//----- previous final step(4 of 4) begin
	for (i = 0; i < 255; i++) { //res[0 ~ (2n-2)-n=n-2=256-2=254]，剩下res[n-1]
		acc[i] = (acc[i] + res_full[i] - res_full[i + n]); // & result_mask;
	}
	acc[n - 1] = (acc[n - 1] + res_full[n - 1]); // & result_mask;
	//----- previous final step(4 of 4) end

	toom4_saber_10bits_reduce_w(w4, ks1_res32, toom4_precision_mask);	//->w4

	for (i = 0; i < sub_len; i++) {
		//--- w2=a2*b2=8A3+4A2+2A1+A0 ---
		a2[i] = a[sub_len3 + i] << 1;	//2A3
		b2[i] = b[sub_len3 + i] << 1;
		a2[i] += a[sub_len2 + i];	//2A3+A2
		b2[i] += b[sub_len2 + i];
		a2[i] <<= 1;	//4A3+2A2
		b2[i] <<= 1;
		a2[i] += a[sub_len + i];	//4A3+2A2+A1
		b2[i] += b[sub_len + i];
		a2[i] <<= 1;	//8A3+2A2+2A1
		b2[i] <<= 1;
		a2[i] += a[i];	//8A3+2A2+2A1+A0
		b2[i] += b[i];

		a2[i] = a2[i] & toom4_precision_mask;
		b2[i] = b2[i] & toom4_precision_mask;
	}
	for (i = 0; i < sub_len; i++) {
		aT[i] = a2[i];
		bT[i] = b2[i];
	}
	speedtest_print_cpucycles("w5_ahalf_bhalf_wait_cpu_used");
	ks1_read(ks1_res32);
	speedtest_print_cpucycles("w5_ahalf_bhalf_read");
	speedtest_print_splitline();

	//---w2---
	ks1_write_int32(aT, bT);	//(a2, b2);
	speedtest_print_cpucycles("w2_a2_b2_write");
	toom4_saber_10bits_reduce_w(w5, ks1_res32, toom4_precision_mask); //->w5
	memset(res_full, 0, sizeof(res_full[0]) * res_full_len);
	for (i = 0; i < w_len; i++) {
		res_full[sub_len6 + i] += w1[i];
	}
	//[step12,14]
	for (i = 0; i < w_len; i++) {
		w4[i] = w4[i] - w3[i]; //step12: w4=(w4-w3)/2
		w4[i] = w4[i] >> 1;
		w3[i] = w3[i] + w4[i]; //step14: w3=w3+w4
	}
	//[step13.1,copy_w5]
	for (i = 0; i < w_len; i++) {
		w5_copy[i] = w5[i];	//copy_w5
		w5[i] = w5[i] - w1[i];	//step13.1:w5=w5-w1;[w5=w5-w1-64*w7]
	}
	for (i = 0; i < sub_len; i++) {
		aT[i] = a0[i];
		bT[i] = b0[i];
	}
	speedtest_print_cpucycles("w2_a2_b2_wait_cpu_used");
	ks1_read(ks1_res32);
	speedtest_print_cpucycles("w2_a2_b2_read");
	speedtest_print_splitline();

	//--w7---
	ks1_write_int32(aT, bT);	//(a0, b0);
	speedtest_print_cpucycles("w7_a0_b0_write");
	toom4_saber_10bits_reduce_w(w2, ks1_res32, toom4_precision_mask);		//->w2
	//[step10,16,17.1]
	for (i = 0; i < w_len; i++) {
		w2[i] = w2[i] + w5_copy[i];	//step10: w2=w2+w5
		w2[i] = w2[i] - (w3[i] << 6) - w3[i];	//step16: w2=w2-65*w3
		w3[i] = w3[i] - w1[i];	//step17.1:w3=w3-w1; [w3=w3-w7-w1]
	}
	for (i = 0; i < sub_len; i++) {
		aT[i] = a_half[i];
		bT[i] = b_half[i];
	}
	speedtest_print_cpucycles("w7_a0_b0_wait_cpu_used");
	ks1_read(ks1_res32);
	speedtest_print_cpucycles("w7_a0_b0_read");
	speedtest_print_splitline();

	//---w6---
	ks1_write_int32(aT, bT);	//(a_half, b_half);
	speedtest_print_cpucycles("w6_a_half_b_half_write");
	toom4_saber_10bits_reduce_w(w7, ks1_res32, toom4_precision_mask);		//->w7
	for (i = 0; i < w_len; i++) {
		res_full[i] += w7[i];
	}
	//[step13.2,17.2,18,21,23]  w2_copy
	for (i = 0; i < w_len; i++) {
		w5[i] = w5[i] - (w7[i] << 6);		//step13.2:w5=w5-64*w7;[w5=w5-w1-64*w7]
		w3[i] = w3[i] - w7[i];		//step17.2:w3=w3-w7; [w3=w3-w7-w1]
		w2[i] = w2[i] + w3[i] * int45; //step18: w2=w2+45*w3
		w2_copy[i] = w2[i]; //w2_copy
		w2[i] = w2[i] + (w4[i] << 4); //step21: w2=(w2+16*w4)/18
		w2[i] = w2[i] * inv9;
		w2[i] = w2[i] >> 1;
		w4[i] = -(w4[i] + w2[i]); //step23: w4=-(w4+w2/2)
	}
	for (i = 0; i < w_len; i++) {
		res_full[sub_len3 + i] += w4[i];
	}
	speedtest_print_cpucycles("w6_a_half_b_half_wait_cpu_used");
	ks1_read(ks1_res32);
	speedtest_print_cpucycles("w6_a_half_b_half_read");
	speedtest_print_splitline();

	//====以下是不能并行的=====
	toom4_saber_10bits_reduce_w(w6, ks1_res32, toom4_precision_mask);
	speedtest_print_cpucycles("reduce_w6"); //->w6

	for (i = 0; i < w_len; i++) {
		w6[i] = w6[i] - w5_copy[i];	//step11: w6=w6-w5
		w5[i] = w6[i] + (w5[i] << 1); //step15: w5=2*w5+w6
		w5[i] = w5[i] - (w3[i] << 3); //step19: w5=(w5-8*w3)/24
		w5[i] = w5[i] * inv3;
		w5[i] = w5[i] >> 3;
		w6[i] = w2_copy[i] + w6[i]; //step20: w6=w6+w2
		w3[i] = w3[i] - w5[i]; //step22. w3 <- w3-w5
		w6[i] = w2[i] * int30 - w6[i]; //step24:  w6=(30*w2-w6)/60
		w6[i] = w6[i] * inv15; //*inv15 / 4
		w6[i] = w6[i] >> 2;
		w2[i] = w2[i] - w6[i]; //step25: w2=w2-w6
	}

	speedtest_print_cpucycles("rest_interpolation");

	for (i = 0; i < w_len; i++) { //w1,w7,w4已经计算过了
		res_full[sub_len + i] += w6[i];
		res_full[sub_len2 + i] += w5[i];
		res_full[sub_len4 + i] += w3[i];
		res_full[sub_len5 + i] += w2[i];
	}
	speedtest_print_cpucycles("res_full+=w2356");
	//---------------reduction-------
	for (i = 0; i < 255; i++) { //res[0 ~ (2n-2)-n=n-2=256-2=254]，剩下res[n-1]
		//result[i] = (res_full[i] - res_full[i + n]) & result_mask;
		acc[i] = (acc[i] + res_full[i] - res_full[i + n]) & result_mask;
	}
	//result[n - 1] = res_full[n - 1] & result_mask;
	acc[n - 1] = (acc[n - 1] + res_full[n - 1]) & result_mask;
	speedtest_print_cpucycles("reduction(X^N+1)");
	speedtest_print_title("toom4_saber_10bits_parallel_vector_speedtest", 0);
}
