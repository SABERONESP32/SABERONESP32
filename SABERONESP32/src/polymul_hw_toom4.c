#include<stdio.h>
#include<stdlib.h>
#include<stdint.h>
#include "bignum_lite.h"
#include "polymul_hw_speed.h"
#include "polymul_hw_common.h"

//2019.1.9
//尝试对于kyber_q=7681在运算过程中进行reduce
//注意：在进行a*b->w时，a和b必须先化为模q的正数
//a,b分别有13位，w的元素最多=13+13+log(64)=32位
//将w化为模q（其实都是正数），为13位
//但是w的变换中需要做乘以64以及乘以(inv60=7553)的计算、减法等，所以需要把w置为有符号的int32

//经测试，reduce_int32_mul_inv_to_int16和reduce_int32_to_uint16带上static inline:96670
//不带 static inline:103382。使用静态内联快很多
//对其他行数较多的函数做静态内联没有效果
//inv最大为2^13，
inline int16_t toom4_w_mul_inv(int32_t v, uint16_t inv, uint16_t q) {
	//这里必须用int64，暂时不知道哪里会超位数.哪怕把a*b->w的w也化为13位，这里也需要int64
	//return ((int64_t) v * inv) % q;//这一步非常耗时
	//return (v % q) * inv;//不对
	//return (((int16_t)(v % q)) * inv) % q;//0.571 ms
	//这里的v是有符号的，用barrett不划算
	return ((int32_t) (((int16_t) (v % q)) * inv)) % q;	//和上面的相同，应该是默认情况下int16相乘得出的结果为int32
	//return v;//直接这么写是  0.426，上面这种正确的是 0.460
//	v = v % REDUCE_Q + REDUCE_Q;
//	v = v * inv;
//	return v - (v >> 13) * REDUCE_Q;
}

//将带符号的int32约化为正数int16。约化w4,w6的带符号a,b; 约化res_full
inline uint16_t toom4_reduce_int32_full(int32_t v, uint16_t q) {
	int16_t t = v % q;
	if ((t >> 15) == 0) {	//比判断>0要快
		return t;
	}
	return t + q;
}

//a,b->[0,q)。适用于w4,w6,其中的元素可能有负数，a,b长度为64！
void toom4_reduce_a_b_signed(int32_t* a, int32_t* b, uint16_t q) {
	int i;
	int16_t t;
	for (i = 0; i < 64; i++) {
//		a[i] = toom4_reduce_int32_full(a[i]);
//		b[i] = toom4_reduce_int32_full(b[i]);
		t = a[i] % q;
		if (t < 0) {
			t += q;
		}
		a[i] = t;
		t = b[i] % q;
		if (t < 0) {
			t += q;
		}
		b[i] = t;
	}
}
//a,b->[0,q)，适用于w2,w3,w5,其中的元素均为正数
void toom4_reduce_a_b_unsigned(int32_t* a, int32_t* b, uint16_t q) {
	int i;
	for (i = 0; i < 64; i++) {
		a[i] = a[i] % q;
		b[i] = b[i] % q;
	}
}

//对硬件乘法器计算出来的uint32的w进行约化，以便后续使用int32来进行计算
void toom4_reduce_w(int32_t* w, uint32_t* ks1_res32, uint16_t q) {
	int i;
	for (i = 0; i < 127; i++) {	//2 * 64 - 1
		w[i] = ks1_res32[i] % q;
	}
}

void IRAM_ATTR ks1_full_uint16_int32_kyber(uint16_t *a, uint16_t *b, int32_t* res, uint16_t q) {
	int i = 0;
	//X低到高->64个X，64个0，
	//Y低到高->64个0，64个Y
	for (i = 0; i < 64; i++) {
		REG_WRITE_NOCHECK(((uint32_t *) REG_BASE_X) + i, (uint32_t ) a[i]);
		REG_WRITE_NOCHECK(((uint32_t *) REG_BASE_Y) + i + 64, (uint32_t ) b[i]);
	}
	bzero(((uint32_t *) REG_BASE_X) + 64, 256);		// 64 * 4);
	bzero(((uint32_t *) REG_BASE_Y), 256);		//64 * 4);	// 64个0
	REG_WRITE_NOCHECK(RSA_MULT_START_REG, 1);	//start
	while (REG_READ_NOCHECK(RSA_INTERRUPT_REG) != 1) {	//wait...
	}
	uint32_t t;
	for (i = 0; i < 127; i++) {	//注意res的长度是2*n-1
		t = REG_READ_NOCHECK(((uint32_t *) REG_BASE_Z) + i);
		res[i] = t % q;
	}
	REG_WRITE_NOCHECK(RSA_INTERRUPT_REG, 1);	//clear interrupt
}

void ks1_full_int32_int32_kyber(int32_t* a, int32_t* b, int32_t* res, uint16_t q) {
	int i = 0;
	//X低到高->64个X，64个0
	memcpy(((uint32_t *) REG_BASE_X), a, 64 * 4);
	bzero(((uint32_t *) REG_BASE_X) + 64, 64 * 4);
	//Y低到高->64个0，64个Y
	bzero(((uint32_t *) REG_BASE_Y), 64 * 4);	// 64个0
	memcpy(((uint32_t *) REG_BASE_Y) + 64, b, 64 * 4);
	REG_WRITE_NOCHECK(RSA_MULT_START_REG, 1);	//start
	while (REG_READ_NOCHECK(RSA_INTERRUPT_REG) != 1) {	//wait...
	}
	uint32_t t;
	for (i = 0; i < 127; i++) {	//注意res的长度是2*n-1
		t = REG_READ_NOCHECK(((uint32_t *) REG_BASE_Z) + i);
		res[i] = t % q;
	}
	REG_WRITE_NOCHECK(RSA_INTERRUPT_REG, 1);	//clear interrupt
}

//串行toom4
void toom4_kyber(uint16_t* a, uint16_t* b, uint16_t* result) {
	//printf("---toom4_kyber---\n");
	REG_WRITE_NOCHECK(RSA_MULT_MODE_REG, 15);
	const uint16_t n = 256;
	const uint16_t q = 7681;
	int16_t i;
	uint16_t sub_len = 64;	//64,n / 4
	uint16_t sub_len2 = 128;
	uint16_t sub_len3 = 192;
	uint16_t sub_len4 = 256;
	uint16_t sub_len5 = 320;
	uint16_t sub_len6 = 384;
	uint16_t w_len = 127;	//2 * sub_len - 1;

	//w需要做乘以8，64，30等运算，用int16不够
	int32_t w1[w_len], w2[w_len], w3[w_len], w4[w_len], w5[w_len], w6[w_len], w7[w_len];			//w_len=128

	int32_t at1, at2;  //temp for w
	int32_t bt1, bt2;

	uint16_t* ainf = &a[sub_len3];	//A(inf),B(inf)->w1, inf:infinite
	uint16_t* binf = &b[sub_len3];
	int32_t a2[sub_len], b2[sub_len];	//A(2)*B(2)->w2
	int32_t a1[sub_len], b1[sub_len];	//A(1)*B(1)->w3
	int32_t a_1[sub_len], b_1[sub_len];	//A(-1)*B(-1)->w4
	int32_t ahalf[sub_len], bhalf[sub_len];	//A(1/2)*B(1/2)->w5
	int32_t a_half[sub_len], b_half[sub_len];	//A(-1/2)*B(-1/2)->w6
	uint16_t* a0 = &a[0];	//A0,B0->w7
	uint16_t* b0 = &b[0];

	//kyber inverse intv3 * 3 = 1 mod 7681 * 8
	uint16_t inv2 = 3841;			//2^12	//=inv(2, 7681)
	uint16_t inv18 = 4694;			//2^13
	uint16_t inv24 = 7361;			//2^13
	uint16_t inv60 = 7553;			//2^13

	int16_t int45 = 45;
	int16_t int30 = 30;

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
		a_half[i] = at1 - at2;	//(8A0 + 2A2) - (4A1 + A3) = 8A0 - 4A1 + 2A2 - A3
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
	}

	ks1_full_uint16_int32_kyber(ainf, binf, w1, q);

	toom4_reduce_a_b_unsigned(a2, b2, q);
	ks1_full_int32_int32_kyber(a2, b2, w2, q);

	toom4_reduce_a_b_unsigned(a1, b1, q);
	ks1_full_int32_int32_kyber(a1, b1, w3, q);

	toom4_reduce_a_b_signed(a_1, b_1, q);
	ks1_full_int32_int32_kyber(a_1, b_1, w4, q);

	toom4_reduce_a_b_unsigned(ahalf, bhalf, q);
	ks1_full_int32_int32_kyber(ahalf, bhalf, w5, q);

	toom4_reduce_a_b_signed(a_half, b_half, q);
	ks1_full_int32_int32_kyber(a_half, b_half, w6, q);

	ks1_full_uint16_int32_kyber(a0, b0, w7, q);

	//----------
	for (i = 0; i < w_len; i++) {
		//step10: w2=w2+w5
		w2[i] = w2[i] + w5[i];	//w2 <- w2+w5

		//step11: w6=w6-w5
		w6[i] = w6[i] - w5[i];	// w6 <- w6-w5

		//step12: w4=(w4-w3)/2
		w4[i] = w4[i] - w3[i];
		w4[i] = toom4_w_mul_inv(w4[i], inv2, q);

		//step13: w5=w5-w1-64*w7
		w5[i] = w5[i] - w1[i] - (w7[i] << 6);	//temp1[i]; // w5 <- w5-64*w7

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
		w5[i] = toom4_w_mul_inv(w5[i], inv24, q);

		//step20: w6=w6+w2
		w6[i] = w2[i] + w6[i]; //w6 <- w6+w2

		//step21: w2=(w2+16*w4)/18
		w2[i] = w2[i] + (w4[i] << 4); //w2 <- w2+16*w4
		w2[i] = toom4_w_mul_inv(w2[i], inv18, q);

		//step22. w3 <- w3-w5
		w3[i] = w3[i] - w5[i];

		//step23: w4=-(w4+w2)
		w4[i] = -(w4[i] + w2[i]);

		//step24:  w6=(30*w2-w6)/60
		w6[i] = w2[i] * int30 - w6[i]; //w6 <- 30*w2-w6
		w6[i] = toom4_w_mul_inv(w6[i], inv60, q);

		//step25: w2=w2-w6
		w2[i] = w2[i] - w6[i];
	}

	int16_t res_full_len = 511; //2 * n - 1;
	int32_t res_full[res_full_len];
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
	memcpy(res_full, w7, w_len * 4);
	memcpy(res_full + sub_len2, w5, w_len * 4);
	memcpy(res_full + sub_len4, w3, w_len * 4);
	memcpy(res_full + sub_len6, w1, w_len * 4);
	res_full[sub_len2 - 1] = 0;
	res_full[sub_len4 - 1] = 0;
	res_full[sub_len6 - 1] = 0;
	for (i = 0; i < w_len; i++) {
		res_full[sub_len + i] += w6[i];
		res_full[sub_len3 + i] += w4[i];
		res_full[sub_len5 + i] += w2[i];
	}
	//---------------reduction-------
	//i从n到2 * n - 1;即从256到511
	for (i = 0; i < 255; i++) { //res[0 ~ (2n-2)-n=n-2=256-2=254]，剩下res[n-1]
		result[i] = toom4_reduce_int32_full(res_full[i] - res_full[i + n], q);
	}
	result[n - 1] = toom4_reduce_int32_full(res_full[n - 1], q);
}

//由toom4_kyber重新排列
void toom4_kyber_parallel(uint16_t* a, uint16_t* b, uint16_t* result) {
	//printf("---toom_cook_4way_kyber_reduce_v3_parallel---\n");
	REG_WRITE_NOCHECK(RSA_MULT_MODE_REG, 15);
	const uint16_t n = 256;
	const uint16_t q = 7681;
	int16_t i;
	uint16_t sub_len = 64;	//64,n / 4
	uint16_t sub_len2 = 128;
	uint16_t sub_len3 = 192;
	uint16_t sub_len4 = 256;
	uint16_t sub_len5 = 320;
	uint16_t sub_len6 = 384;
	uint16_t w_len = 127;	//2 * sub_len - 1;

	int32_t at1, at2;  //temp
	int32_t bt1, bt2;

	//kyber inverse intv3 * 3 = 1 mod 7681 * 8
	uint16_t inv2 = 3841;			//2^12	//=inv(2, 7681)
	uint16_t inv18 = 4694;			//2^13
	uint16_t inv24 = 7361;			//2^13
	uint16_t inv60 = 7553;			//2^13

	int16_t int45 = 45;
	//int16_t int30 = 30;

	uint16_t* ainf = &a[sub_len3];	//A(inf),B(inf)->w1, inf:infinite
	uint16_t* binf = &b[sub_len3];

	//---w1---
	ks1_write_uint16(ainf, binf);
	int16_t res_full_len = 511;	//2 * n - 1;
	int32_t res_full[res_full_len];
	//①memcpy(w7,w5,w3,w1);②res[127,255,383]=0;③res+=(w2,w4,w6)
	res_full[sub_len2 - 1] = 0;
	res_full[sub_len4 - 1] = 0;
	res_full[sub_len6 - 1] = 0;
	int32_t a2[sub_len], b2[sub_len];	//A(2)*B(2)->w2
	int32_t a1[sub_len], b1[sub_len];	//A(1)*B(1)->w3
	int32_t a_1[sub_len], b_1[sub_len];	//A(-1)*B(-1)->w4
	int32_t ahalf[sub_len], bhalf[sub_len];	//A(1/2)*B(1/2)->w5
	int32_t a_half[sub_len], b_half[sub_len];	//A(-1/2)*B(-1/2)->w6
	//uint16_t* a0 = &a[0];	//A0,B0->w7 提前做好int32的a0和b0,write更快
	//uint16_t* b0 = &b[0];
	//w需要做乘以8，64，30等运算，用int16不够
	int32_t w1[w_len], w2[w_len], w3[w_len], w4[w_len], w5[w_len], w6[w_len], w7[w_len];
	int32_t w5_copy[w_len];	//备份w5给step10,11
	int32_t w2_copy[w_len];	//在step21之前备份w5给step20
	uint32_t ks1_res32[w_len];	//存放从寄存器读取的ks1结果
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
	}
	toom4_reduce_a_b_unsigned(a1, b1, q);
	toom4_reduce_a_b_signed(a_1, b_1, q);
	ks1_read(ks1_res32);

	//---w3---
	ks1_write_int32(a1, b1);
	toom4_reduce_w(w1, ks1_res32, q);	//->w1
	//w3这里cpu耗时有点多，可以把w1->res_full放到下一个
	memcpy(res_full + sub_len6, w1, w_len * 4); //copy w1
	for (i = 0; i < sub_len; i++) {
		//---w5,w6=a±half*b±half=A(±1/2)*B(±1/2)=(8A0+2A2)±(4A1+A3)---
		at1 = a[i] << 2;	//4*A0
		bt1 = b[i] << 2;
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
		ahalf[i] = at1 + at2;	//(8A0+2A2)+(4A1+A3) = 8A0+4A1+2A2+A3
		bhalf[i] = bt1 + bt2;
		a_half[i] = at1 - at2;	//(8A0+2A2)-(4A1+A3) = 8A0-4A1-2A2+A3
		b_half[i] = bt1 - bt2;
	}
	//w3这里cpu消耗比较多，可以把对ahalf(w5)和a_half(w6)的约化放到下一个
	ks1_read(ks1_res32);

	//---w4---
	ks1_write_int32(a_1, b_1);
	toom4_reduce_w(w3, ks1_res32, q);	//->w3
	//这里约化ahalf(w5)和a_half(w6)的话CPU时间会超
	//w4之后是计算w5，这里必须先把ahalf约化了，对于a_half(w6)往后放
	toom4_reduce_a_b_unsigned(ahalf, bhalf, q);	//(for w5)
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
	}
	//这里约化a2(w2)时间会超，往后放
	ks1_read(ks1_res32);

	//---w5---
	ks1_write_int32(ahalf, bhalf);
	toom4_reduce_w(w4, ks1_res32, q);	//->w4
	toom4_reduce_a_b_unsigned(a2, b2, q);	// for w2
	//这里约化a_half(w6)时间会超，往后放
	//[step12,14]
	for (i = 0; i < w_len; i++) {
		//step12: w4=(w4-w3)/2
		w4[i] = w4[i] - w3[i];
		w4[i] = toom4_w_mul_inv(w4[i], inv2, q);
		//step14: w3=w3+w4
		w3[i] = w3[i] + w4[i];
	}
	ks1_read(ks1_res32);

	//---w2---
	ks1_write_int32(a2, b2);
	toom4_reduce_w(w5, ks1_res32, q); //->w5
	toom4_reduce_a_b_signed(a_half, b_half, q);	//(for w6)
	//[step13.1,copy_w5]
	for (i = 0; i < w_len; i++) {
		w5_copy[i] = w5[i];	//copy_w5
		//step13.1:w5=w5-w1;[w5=w5-w1-64*w7]
		w5[i] = w5[i] - w1[i];
	}
	for (i = 0; i < sub_len; i++) {
		a2[i] = a[i];		//a0 for w7
		b2[i] = b[i];		//b0
	}
	ks1_read(ks1_res32);

	//--w7---
	ks1_write_int32(a2, b2);		//int32 of (a0, b0);
	toom4_reduce_w(w2, ks1_res32, q);		//->w2
	//[step10,16,17.1]
	for (i = 0; i < w_len; i++) {
		//step10: w2=w2+w5
		w2[i] = w2[i] + w5_copy[i];
		//step16: w2=w2-65*w3
		w2[i] = w2[i] - (w3[i] << 6) - w3[i];
		//step17.1:w3=w3-w1; [w3=w3-w7-w1]
		w3[i] = w3[i] - w1[i];
	}
	ks1_read(ks1_res32);

	//---w6---
	ks1_write_int32(a_half, b_half);
	toom4_reduce_w(w7, ks1_res32, q);		//->w7
	memcpy(res_full, w7, w_len * 4);		//copy w7
	//[step13.2,17.2,18,21,23]  w2_copy
	for (i = 0; i < w_len; i++) {
		//step13.2:w5=w5-64*w7;[w5=w5-w1-64*w7]
		w5[i] = w5[i] - (w7[i] << 6);
		//step17.2:w3=w3-w7; [w3=w3-w7-w1]
		w3[i] = w3[i] - w7[i];
		//step18: w2=w2+45*w3
		w2[i] = w2[i] + w3[i] * int45;
		w2_copy[i] = w2[i]; //w2_copy
		//step21: w2=(w2+16*w4)/18
		w2[i] = w2[i] + (w4[i] << 4);
		w2[i] = toom4_w_mul_inv(w2[i], inv18, q);
		//step23: w4=-(w4+w2/2)
		w4[i] = -(w4[i] + w2[i]);
	}
	//	for (i = 0; i < w_len; i++) {
	//如果在这里把w4放进res_full，第一是这一步cpu耗时超长，
	//第二是影响w5和w3不能用memcpy进res_full
	//		res_full[sub_len3 + i] += w4[i];
	//	}
	ks1_read(ks1_res32);

	//====以下是不能并行的=====
	toom4_reduce_w(w6, ks1_res32, q);		//->w6

	for (i = 0; i < w_len; i++) {
		//step11: w6=w6-w5
		w6[i] = w6[i] - w5_copy[i];

		//step15: w5=2*w5+w6
		w5[i] = w6[i] + (w5[i] << 1);

		//step19: w5=(w5-8*w3)/24
		w5[i] = w5[i] - (w3[i] << 3);
		w5[i] = toom4_w_mul_inv(w5[i], inv24, q);

		//step20: w6=w6+w2
		w6[i] = w2_copy[i] + w6[i];

		//step22. w3 <- w3-w5
		w3[i] = w3[i] - w5[i];

		//step24:  w6=(30*w2-w6)/60
		w6[i] = (w2[i] << 5) - (w2[i] << 1) - w6[i];		//w2[i] * int30 - w6[i];
		w6[i] = toom4_w_mul_inv(w6[i], inv60, q);

		//step25: w2=w2-w6
		w2[i] = w2[i] - w6[i];
	}
	//w1,w7已经copy
	memcpy(res_full + sub_len2, w5, w_len * 4);
	memcpy(res_full + sub_len4, w3, w_len * 4);

	for (i = 0; i < w_len; i++) { //w1,w7,w4已经计算过了
		res_full[sub_len + i] += w6[i];
		res_full[sub_len3 + i] += w4[i];
		res_full[sub_len5 + i] += w2[i];
	}
	//---------------reduction-------
	//i从n到2 * n - 1;即从256到511
	for (i = 0; i < 255; i++) { //res[0 ~ (2n-2)-n=n-2=256-2=254]，剩下res[n-1]
		result[i] = toom4_reduce_int32_full(res_full[i] - res_full[i + n], q);
	}
	result[n - 1] = toom4_reduce_int32_full(res_full[n - 1], q);
}

void toom4_kyber_parallel_speedtest(uint16_t* a, uint16_t* b, uint16_t* result) {
//	speedtest_print_title("toom4_kyber_parallel_speedtest", 1);
//	printf("REDUCE_Q=%d\n", REDUCE_Q);
//	speedtest_reset_startcpucycles();
//	speedtest_print_splitline();
//	speedtest_print_splitline();
//	REG_WRITE_NOCHECK(RSA_MULT_MODE_REG, 15);
//	speedtest_print_cpucycles("RSA_MULT_MODE_REG");
//	int16_t n = 256;
//	int16_t i;
//	uint16_t sub_len = 64;	//64,n / 4
//	uint16_t sub_len2 = 128;
//	uint16_t sub_len3 = 192;
//	uint16_t sub_len4 = 256;
//	uint16_t sub_len5 = 320;
//	uint16_t sub_len6 = 384;
//	uint16_t w_len = 127;	//2 * sub_len - 1;
//
//	int32_t at1, at2;  //temp
//	int32_t bt1, bt2;
//
//	//kyber inverse intv3 * 3 = 1 mod 7681 * 8
//	uint16_t inv2 = 3841;			//2^12	//=inv(2, 7681)
//	uint16_t inv18 = 4694;			//2^13
//	uint16_t inv24 = 7361;			//2^13
//	uint16_t inv60 = 7553;			//2^13
//
//	int16_t int45 = 45;
//	//int16_t int30 = 30;
//
//	uint16_t* ainf = &a[sub_len3];	//A(inf),B(inf)->w1, inf:infinite
//	uint16_t* binf = &b[sub_len3];
//
//	speedtest_print_cpucycles("init_all_pointers");
//	speedtest_print_splitline();
//
//	//---w1---
//	ks1_write_uint16(ainf, binf);
//	speedtest_print_cpucycles("w1_ainf_binf_write");
//	int16_t res_full_len = 511;	//2 * n - 1;
//	int32_t res_full[res_full_len];
//	//①memcpy(w7,w5,w3,w1);②res[127,255,383]=0;③res+=(w2,w4,w6)
//	res_full[sub_len2 - 1] = 0;
//	res_full[sub_len4 - 1] = 0;
//	res_full[sub_len6 - 1] = 0;
//	int32_t a2[sub_len], b2[sub_len];	//A(2)*B(2)->w2
//	int32_t a1[sub_len], b1[sub_len];	//A(1)*B(1)->w3
//	int32_t a_1[sub_len], b_1[sub_len];	//A(-1)*B(-1)->w4
//	int32_t ahalf[sub_len], bhalf[sub_len];	//A(1/2)*B(1/2)->w5
//	int32_t a_half[sub_len], b_half[sub_len];	//A(-1/2)*B(-1/2)->w6
//	//uint16_t* a0 = &a[0];	//A0,B0->w7 提前做好int32的a0和b0,write更快
//	//uint16_t* b0 = &b[0];
//	//w需要做乘以8，64，30等运算，用int16不够
//	int32_t w1[w_len], w2[w_len], w3[w_len], w4[w_len], w5[w_len], w6[w_len], w7[w_len];
//	int32_t w5_copy[w_len];	//备份w5给step10,11
//	int32_t w2_copy[w_len];	//在step21之前备份w5给step20
//	uint32_t ks1_res32[w_len];	//存放从寄存器读取的ks1结果
//	for (i = 0; i < sub_len; i++) {
//		//--- w3,w4=a(±1)*b(±1)=(A2+A0)±(A3+A1) ---
//		at1 = a[sub_len2 + i] + a[i];	//A2 + A0
//		bt1 = b[sub_len2 + i] + b[i];
//		at2 = a[sub_len3 + i] + a[sub_len + i];	//A3 + A1
//		bt2 = b[sub_len3 + i] + b[sub_len + i];
//		a1[i] = at1 + at2; //(A2 + A0) + (A3 + A1) = A0 + A1 + A2 + A3
//		b1[i] = bt1 + bt2;
//		a_1[i] = at1 - at2;	//(A2 + A0) - (A3 + A1) = A0 - A1 + A2 - A3
//		b_1[i] = bt1 - bt2;
//	}
//	toom4_reduce_a_b_unsigned(a1, b1);
//	toom4_reduce_a_b_signed(a_1, b_1);
//	speedtest_print_cpucycles("w1_ainf_binf_wait_cpu_used");
//	ks1_read(ks1_res32);
//	speedtest_print_cpucycles("w1_ainf_binf_read");
//	speedtest_print_splitline();
//
//	//---w3---
//	ks1_write_int32(a1, b1);
//	speedtest_print_cpucycles("w3_a1_b1_write");
//	toom4_reduce_w(w1, ks1_res32);	//->w1
//	//w3这里cpu耗时有点多，可以把w1->res_full放到下一个
//	memcpy(res_full + sub_len6, w1, w_len * 4); //copy w1
//	for (i = 0; i < sub_len; i++) {
//		//---w5,w6=a±half*b±half=A(±1/2)*B(±1/2)=(8A0+2A2)±(4A1+A3)---
//		at1 = a[i] << 2;	//4*A0
//		bt1 = b[i] << 2;
//		at1 = at1 + a[sub_len2 + i];	//4A0+A2
//		bt1 = bt1 + b[sub_len2 + i];
//		at1 = at1 << 1;	//8A0 + 2A2
//		bt1 = bt1 << 1;
//		at2 = a[sub_len + i];	//A1
//		bt2 = b[sub_len + i];
//		at2 = at2 << 2;	//4A1
//		bt2 = bt2 << 2;
//		at2 = at2 + a[sub_len3 + i];	//4A1 + A3
//		bt2 = bt2 + b[sub_len3 + i];
//		ahalf[i] = at1 + at2;	//(8A0+2A2)+(4A1+A3) = 8A0+4A1+2A2+A3
//		bhalf[i] = bt1 + bt2;
//		a_half[i] = at1 - at2;	//(8A0+2A2)-(4A1+A3) = 8A0-4A1-2A2+A3
//		b_half[i] = bt1 - bt2;
//	}
//	//w3这里cpu消耗比较多，可以把对ahalf(w5)和a_half(w6)的约化放到下一个
//	speedtest_print_cpucycles("w3_a1_b1_wait_cpu_used");
//	ks1_read(ks1_res32);
//	speedtest_print_cpucycles("w3_a1_b1_read");
//	speedtest_print_splitline();
//
//	//---w4---
//	ks1_write_int32(a_1, b_1);
//	speedtest_print_cpucycles("w4_a_1_b_1_write");
//	toom4_reduce_w(w3, ks1_res32);	//->w3
//	//这里约化ahalf(w5)和a_half(w6)的话CPU时间会超
//	//w4之后是计算w5，这里必须先把ahalf约化了，对于a_half(w6)往后放
//	toom4_reduce_a_b_unsigned(ahalf, bhalf);	//(for w5)
//	for (i = 0; i < sub_len; i++) {
//		//--- w2=a2*b2=8A3+4A2+2A1+A0 ---
//		a2[i] = a[sub_len3 + i] << 1;	//2A3
//		b2[i] = b[sub_len3 + i] << 1;
//		a2[i] += a[sub_len2 + i];	//2A3+A2
//		b2[i] += b[sub_len2 + i];
//		a2[i] <<= 1;	//4A3+2A2
//		b2[i] <<= 1;
//		a2[i] += a[sub_len + i];	//4A3+2A2+A1
//		b2[i] += b[sub_len + i];
//		a2[i] <<= 1;	//8A3+2A2+2A1
//		b2[i] <<= 1;
//		a2[i] += a[i];	//8A3+2A2+2A1+A0
//		b2[i] += b[i];
//	}
//	//这里约化a2(w2)时间会超，往后放
//	speedtest_print_cpucycles("w4_a_1_b_1_wait_cpu_used");
//	ks1_read(ks1_res32);
//	speedtest_print_cpucycles("w4_a_1_b_1_read");
//	speedtest_print_splitline();
//
//	//---w5---
//	ks1_write_int32(ahalf, bhalf);
//	speedtest_print_cpucycles("w5_ahalf_bhalf_write");
//	toom4_reduce_w(w4, ks1_res32);	//->w4
//	toom4_reduce_a_b_unsigned(a2, b2);	// for w2
//	//这里约化a_half(w6)时间会超，往后放
//	//[step12,14]
//	for (i = 0; i < w_len; i++) {
//		//step12: w4=(w4-w3)/2
//		w4[i] = w4[i] - w3[i];
//		w4[i] = toom4_w_mul_inv(w4[i], inv2);
//		//step14: w3=w3+w4
//		w3[i] = w3[i] + w4[i];
//	}
//	speedtest_print_cpucycles("w5_ahalf_bhalf_wait_cpu_used");
//	ks1_read(ks1_res32);
//	speedtest_print_cpucycles("w5_ahalf_bhalf_read");
//	speedtest_print_splitline();
//
//	//---w2---
//	ks1_write_int32(a2, b2);
//	speedtest_print_cpucycles("w2_a2_b2_write");
//	toom4_reduce_w(w5, ks1_res32); //->w5
//	toom4_reduce_a_b_signed(a_half, b_half);	//(for w6)
//	//[step13.1,copy_w5]
//	for (i = 0; i < w_len; i++) {
//		w5_copy[i] = w5[i];	//copy_w5
//		//step13.1:w5=w5-w1;[w5=w5-w1-64*w7]
//		w5[i] = w5[i] - w1[i];
//	}
//	for (i = 0; i < sub_len; i++) {
//		a2[i] = a[i];		//a0 for w7
//		b2[i] = b[i];		//b0
//	}
//	speedtest_print_cpucycles("w2_a2_b2_wait_cpu_used");
//	ks1_read(ks1_res32);
//	speedtest_print_cpucycles("w2_a2_b2_read");
//	speedtest_print_splitline();
//
//	//--w7---
//	ks1_write_int32(a2, b2);		//int32 of (a0, b0);
//	speedtest_print_cpucycles("w7_a0_b0_write");
//	toom4_reduce_w(w2, ks1_res32);		//->w2
//	//[step10,16,17.1]
//	for (i = 0; i < w_len; i++) {
//		//step10: w2=w2+w5
//		w2[i] = w2[i] + w5_copy[i];
//		//step16: w2=w2-65*w3
//		w2[i] = w2[i] - (w3[i] << 6) - w3[i];
//		//step17.1:w3=w3-w1; [w3=w3-w7-w1]
//		w3[i] = w3[i] - w1[i];
//	}
//	speedtest_print_cpucycles("w7_a0_b0_wait_cpu_used");
//	ks1_read(ks1_res32);
//	speedtest_print_cpucycles("w7_a0_b0_read");
//	speedtest_print_splitline();
//
//	//---w6---
//	ks1_write_int32(a_half, b_half);
//	speedtest_print_cpucycles("w6_a_half_b_half_write");
//	toom4_reduce_w(w7, ks1_res32);		//->w7
//	memcpy(res_full, w7, w_len * 4);		//copy w7
//	//[step13.2,17.2,18,21,23]  w2_copy
//	for (i = 0; i < w_len; i++) {
//		//step13.2:w5=w5-64*w7;[w5=w5-w1-64*w7]
//		w5[i] = w5[i] - (w7[i] << 6);
//		//step17.2:w3=w3-w7; [w3=w3-w7-w1]
//		w3[i] = w3[i] - w7[i];
//		//step18: w2=w2+45*w3
//		w2[i] = w2[i] + w3[i] * int45;
//		w2_copy[i] = w2[i]; //w2_copy
//		//step21: w2=(w2+16*w4)/18
//		w2[i] = w2[i] + (w4[i] << 4);
//		w2[i] = toom4_w_mul_inv(w2[i], inv18);
//		//step23: w4=-(w4+w2/2)
//		w4[i] = -(w4[i] + w2[i]);
//	}
//	//	for (i = 0; i < w_len; i++) {
//	//如果在这里把w4放进res_full，第一是这一步cpu耗时超长，
//	//第二是影响w5和w3不能用memcpy进res_full
//	//		res_full[sub_len3 + i] += w4[i];
//	//	}
//	speedtest_print_cpucycles("w6_a_half_b_half_wait_cpu_used");
//	ks1_read(ks1_res32);
//	speedtest_print_cpucycles("w6_a_half_b_half_read");
//	speedtest_print_splitline();
//
//	//====以下是不能并行的=====
//	toom4_reduce_w(w6, ks1_res32);		//->w6
//	speedtest_print_cpucycles("reduce_w6");
//
//	for (i = 0; i < w_len; i++) {
//		//step11: w6=w6-w5
//		w6[i] = w6[i] - w5_copy[i];
//
//		//step15: w5=2*w5+w6
//		w5[i] = w6[i] + (w5[i] << 1);
//
//		//step19: w5=(w5-8*w3)/24
//		w5[i] = w5[i] - (w3[i] << 3);
//		w5[i] = toom4_w_mul_inv(w5[i], inv24);
//
//		//step20: w6=w6+w2
//		w6[i] = w2_copy[i] + w6[i];
//
//		//step22. w3 <- w3-w5
//		w3[i] = w3[i] - w5[i];
//
//		//step24:  w6=(30*w2-w6)/60
//		w6[i] = (w2[i] << 5) - (w2[i] << 1) - w6[i];		//w2[i] * int30 - w6[i];
//		w6[i] = toom4_w_mul_inv(w6[i], inv60);
//
//		//step25: w2=w2-w6
//		w2[i] = w2[i] - w6[i];
//	}
//	speedtest_print_cpucycles("rest_interpolation");
//	//w1,w7已经copy
//	memcpy(res_full + sub_len2, w5, w_len * 4);
//	memcpy(res_full + sub_len4, w3, w_len * 4);
//
//	for (i = 0; i < w_len; i++) { //w1,w7,w4已经计算过了
//		res_full[sub_len + i] += w6[i];
//		res_full[sub_len3 + i] += w4[i];
//		res_full[sub_len5 + i] += w2[i];
//	}
//	speedtest_print_cpucycles("res_full+=w23456");
//	//---------------reduction-------
//	//i从n到2 * n - 1;即从256到511
//	for (i = 0; i < 255; i++) { //res[0 ~ (2n-2)-n=n-2=256-2=254]，剩下res[n-1]
//		result[i] = toom4_reduce_int32_full(res_full[i] - res_full[i + n]);
//	}
//	result[n - 1] = toom4_reduce_int32_full(res_full[n - 1]);
//	speedtest_print_cpucycles("reduction(X^N+1)");
//	speedtest_print_title("toom4_kyber_parallel_speedtest", 0);

}

