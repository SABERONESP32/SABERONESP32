/*
 * polymul_hw_toom4_general.c
 *
 *  Created on: 2019年1月16日
 *      Author: Mai
 */
#include "polymul_hw_common.h"
#include "polymul_hw.h"
#include "polymul_hw_speed.h"
#include "bignum_lite.h"

mpi h0, h1;	//长度为2n
mpi h_add, h_sub;	//长度为2n

#define KS2_1536 1

void toom4_general_ks2_mpi_init() {
	int n2 = 48 * 2;
	printf("INIT KS2_1536\n");
	mpi_init(&h0);		//f0*g0
	mpi_init(&h1);		//f1*g1
	mpi_init(&h_add);
	mpi_init(&h_sub);
	h0.s = 1;
	h1.s = 1;
	mpi_grow(&h0, n2);
	mpi_grow(&h1, n2);
	mpi_grow(&h_add, n2);
	mpi_grow(&h_sub, n2);
	printf("#toom4_general_ks2_mpi_init#\n");

}

void toom4_general_ks2_mpi_free() {
	mpi_free(&h0);
	mpi_free(&h1);
	mpi_free(&h_add);
	mpi_free(&h_sub);
	printf("#toom4_general_ks2_mpi_free#\n");
}

//n=64
void toom4_general_fixborrow_f0g0(uint32_t *f0g0_fix, int32_t *a, int n) {
	int32_t x = 16777216; //2^24，带入x是24位的
	int32_t t;
	int borrow = 0;	//借位，+1表示借了一位
	//int sign;
	int32_t i;
	//求两个负点值，这一部分亲测几乎不耗时
	for (i = 0; i < n; i++) {		//f(-100),g(-100)
		//偶,系数为负, =-self-borrow//奇,系数为正, self=self-borrow
		//i%2==0:t = -a[i] - borrow;
		//i%2==1: t = a[i] - borrow;
		//sign = i % 2 == 0 ? -1 : 1;
		t = a[i] - borrow;		//(sign * a[i]) - borrow;
		if (t < 0) {
			borrow = 1;
			t += x;
		} else {
			borrow = 0;
		}
		f0g0_fix[i] = t;
	}
}

//使用fixborrow_f1g1处理借位，输入a，输出f1gi1_fix，然后使用insert_f1g1计算f1和g1
void toom4_general_fixborrow_f1g1(uint32_t *f1g1_fix, int32_t *a, int n) {
	int32_t x = 16777216; //2^24，带入x是24位的
	int32_t t;
	int borrow = 0;	//借位，+1表示借了一位
	int sign;
	int32_t i;
	//求两个负点值，这一部分亲测几乎不耗时
	for (i = 0; i < n; i++) {		//f(-100),g(-100)
		//偶,系数为负, =-self-borrow//奇,系数为正, self=self-borrow
		//i%2==0:t = -a[i] - borrow;
		//i%2==1: t = a[i] - borrow;
		sign = i % 2 == 0 ? -1 : 1;
		t = (sign * a[i]) - borrow;
		if (t < 0) {
			borrow = 1;
			t += x;
		} else {
			borrow = 0;
		}
		f1g1_fix[i] = t;
		//实测t都是正数
//		if (t < 0) {
//			printf("WARNING!! t<0 in fullbits_fixborrow_f1g1\n");
//		}
	}
}

//输入的a必须是已经处理过借位的a
//1536bit=48B，将64个uint32存入48B寄存器
void ks2_write_f0g0f1g1(uint32_t *a, uint32_t *reg_base) {
	int16_t i = 0;	//index of mpi
	int16_t j = 0;	//index of poly
	int mpi_n = 48;
	uint32_t t;
	while (i < mpi_n) {
		t = 0;
		t |= a[j + 1] & 0xFF;	//8位
		t = t << 24; //高8位
		t |= a[j] & 0xFFFFFF; //24位
		REG_WRITE(reg_base + i, t);
		//----------------------
		t = 0;
		t |= a[j + 2] & 0xFFFF; //16
		t = t << 16; //高16位为j+2
		t |= (a[j + 1] >> 8) & 0xFFFF; //低16位为j+1
		REG_WRITE(reg_base + i + 1, t);
		//----------------------
		t = 0;
		t |= a[j + 3] & 0xFFFFFF;		//24位
		t = t << 8;		//高24位
		t |= (a[j + 2] >> 16) & 0xFF;		//8低8位
		REG_WRITE(reg_base + i + 2, t);
		i += 3;
		j += 4;
	}
}

//a,b.len=64; res.len=127
void ks2_full(int32_t* a, int32_t* b, int64_t* res, int fixcarry) {	//len=64

	int n = 64;	//用mpi_ks2代替kara，共64个数据
	int mpi_n = 48;	//每个实际需要18*2+6=42位，对齐是48位，ks2折半是24位,
	int mpi_n2 = mpi_n * 2;
	int16_t i;
	//24bit * 64 / 32bit = mpi中需要48个int32
	uint32_t borrow_fix_temp[n];	//借位

	//REG_WRITE(RSA_MULT_MODE_REG, 13);

	//---------计算h0=f0*g0，两个正点值f0和g0的乘积--------------
	//写入正点值f0和g0
	toom4_general_fixborrow_f0g0(borrow_fix_temp, a, n);
	ks2_write_f0g0f1g1(borrow_fix_temp, ((uint32_t *) REG_BASE_X));
	toom4_general_fixborrow_f0g0(borrow_fix_temp, b, n);
	ks2_write_f0g0f1g1(borrow_fix_temp, ((uint32_t *) REG_BASE_Y) + 48);	//注意，从y的48字开始写数据(前48字填0)
	//reg写入0位
	//x低到高：48个X，48个0
	bzero(((uint32_t *) REG_BASE_X) + 48, 48 * 4);
	//y低到高：48个0，48个Y
	bzero(((uint32_t *) REG_BASE_Y), 48 * 4);
	//求正点值的乘积
	REG_WRITE_NOCHECK(RSA_MULT_START_REG, 1);	//start
	while (REG_READ_NOCHECK(RSA_INTERRUPT_REG) != 1) {
	}
	//读出h0,h0=f0*g0
	for (i = 0; i < mpi_n2; i++) {	//注意h0为2倍的mpi_n=48*2
		h0.p[i] = REG_READ_NOCHECK(((uint32_t *) REG_BASE_Z) + i);
	}
	REG_WRITE_NOCHECK(RSA_INTERRUPT_REG, 1);
	//esp_mpi_release_hardware();//over

	//---------计算h1=f1*g1，两个负点的乘积--------------

	//写入负点值f1g1
	toom4_general_fixborrow_f1g1(borrow_fix_temp, a, n);
	ks2_write_f0g0f1g1(borrow_fix_temp, ((uint32_t *) REG_BASE_X));
	toom4_general_fixborrow_f1g1(borrow_fix_temp, b, n);
	ks2_write_f0g0f1g1(borrow_fix_temp, ((uint32_t *) REG_BASE_Y) + 48);	//注意，从y的48字开始写数据(前48字填0)

	//reg写入0位
	//x低到高：48个X，48个0
	bzero(((uint32_t *) REG_BASE_X) + 48, 48 * 4);
	//y低到高：48个0，48个Y
	bzero(((uint32_t *) REG_BASE_Y), 48 * 4);

	REG_WRITE_NOCHECK(RSA_MULT_START_REG, 1);	//start
	while (REG_READ_NOCHECK(RSA_INTERRUPT_REG) != 1) {
	}

	//读出h1,h1=f1*g1
	for (i = 0; i < mpi_n2; i++) {	//注意h1为2倍的mpi_n
		h1.p[i] = REG_READ_NOCHECK(((uint32_t *) REG_BASE_Z) + i);
	}
	REG_WRITE_NOCHECK(RSA_INTERRUPT_REG, 1);
	//esp_mpi_release_hardware();

	//-----------h0和h1计算完毕--------------------

	mpi_add_abs_lite(&h_add, &h0, &h1);
	mpi_sub_abs_lite(&h_sub, &h0, &h1);

	mpi_shift_right_1bit_lite(&h_add);		//h0/2
	mpi_shift_right_1bit_lite(&h_sub);		//h1/2

	int32_t i_mpi = 0;
	int32_t i_poly = 0;
	//res中每个数为48位
	while (i_poly < 2 * n) {	//h_add，取出偶数项
		//i_mpi //0
		res[i_poly] = h_add.p[i_mpi + 1] & 0xFFFF;		//res[0]高16位
		res[i_poly] = res[i_poly] << 32;
		res[i_poly] |= h_add.p[i_mpi];		//res[0]低32位
		//--------------------------
		i_poly += 2;
		i_mpi += 1;	//1
		res[i_poly] = h_add.p[i_mpi + 1];	//res[2]高32位
		res[i_poly] = res[i_poly] << 16;
		res[i_poly] |= (h_add.p[i_mpi] >> 16) & 0xFFFF;	//res[2] 低16位
		i_mpi += 2;	//3
		i_poly += 2;
	}
	i_mpi = 0;
	i_poly = 1;	//从第1项开始
	while (i_poly < 2 * n) {	//h_sub，取出奇数项
		//i_mpi // 0
		res[i_poly] = h_sub.p[i_mpi + 2] & 0xFF;		//res[1] high 8
		res[i_poly] = res[i_poly] << 40;
		res[i_poly] |= ((int64_t) h_sub.p[i_mpi + 1]) << 8;		// res[1] middle 32
		res[i_poly] |= (h_sub.p[i_mpi] >> 24) & 0xFF;	//res[1] low 8
		//-------------------------------
		i_mpi += 2;	//2
		i_poly += 2;	//3
		res[i_poly] = h_sub.p[i_mpi + 1] & 0xFFFFFF;	//res[3] high 24
		res[i_poly] = res[i_poly] << 24;
		res[i_poly] |= (h_sub.p[i_mpi] >> 8) & 0xFFFFFF;	//res[3] low 24
		i_mpi += 1;	//3
		i_poly += 2;	//5
	}

	if (!fixcarry) {
		return;
	}

	//下面是处理进位的
	int carry0 = 0;		//偶数位进位
	int carry1 = 0;		//奇数位进位

	int64_t x = 281474976710656;		//x=2^48
	int64_t x_half = 140737488355328;		//x_half=x/2

	for (i = 0; i < (n * 2); i++) {
		if (i % 2 == 0) {		//h_add
			res[i] += carry0;
			//res[i] > x_half等价于carry0 = res[i] < 0 但不严谨
			if (res[i] > x_half) {		////carry1 = res[i] < 0;
				res[i] = res[i] - x;
				carry0 = 1;
			} else {
				carry0 = 0;
			}
		} else {  //h_sub
			res[i] += carry1;
			if (res[i] > x_half) {
				res[i] = res[i] - x;
				carry1 = 1;
			} else {
				carry1 = 0;
			}
		}
	}

}

void ks1_full_uint16_int64(uint16_t *a, uint16_t *b, int64_t* res) {
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
	for (i = 0; i < 127; i++) {	//注意res的长度是2*n-1
		res[i] = REG_READ_NOCHECK(((uint32_t *) REG_BASE_Z) + i);
	}
	REG_WRITE_NOCHECK(RSA_INTERRUPT_REG, 1);	//clear interrupt
}

inline uint16_t toom4_general_reduce_int64_to_uint16(int64_t v) {
	int16_t t = v % REDUCE_Q;
	if ((t >> 15) == 0) {	//比判断>0要快
		return t;
	}
	return t + REDUCE_Q;
}

//通用版的toom4，适用于所有的Q, w1和w7是KS1，w2~w6是KS2
void toom4_general(uint16_t* a, uint16_t* b, uint16_t* result) {
	//a,b最大为带符号的15a，即1+4+13=18位
	//w为64个a，b相乘，最大为6+18+18=42位
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

	//w需要做乘以8，64，30等运算，用int16不够
	int64_t w1[w_len], w2[w_len], w3[w_len], w4[w_len], w5[w_len], w6[w_len], w7[w_len];			//w_len=128

	int64_t at1, at2;  //temp for w
	int64_t bt1, bt2;

	uint16_t* ainf = &a[sub_len3];	//A(inf),B(inf)->w1, inf:infinite
	uint16_t* binf = &b[sub_len3];
	int32_t a2[sub_len], b2[sub_len];	//A(2)*B(2)->w2
	int32_t a1[sub_len], b1[sub_len];	//A(1)*B(1)->w3
	int32_t a_1[sub_len], b_1[sub_len];	//A(-1)*B(-1)->w4
	int32_t ahalf[sub_len], bhalf[sub_len];	//A(1/2)*B(1/2)->w5
	int32_t a_half[sub_len], b_half[sub_len];	//A(-1/2)*B(-1/2)->w6
	uint16_t* a0 = &a[0];	//A0,B0->w7
	uint16_t* b0 = &b[0];

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
	}

	REG_WRITE_NOCHECK(RSA_MULT_MODE_REG, 15);
	ks1_full_uint16_int64(ainf, binf, w1); //  w1=A(∞)*B(∞)  ks1 OK
	ks1_full_uint16_int64(a0, b0, w7); //  w7=A(0)*B(0)   ks1 OK

	REG_WRITE_NOCHECK(RSA_MULT_MODE_REG, 13);
	ks2_full(a2, b2, w2, 0); //          w2=A(2)*B(2)
	ks2_full(a1, b1, w3, 0); //  w3=A(1)*B(1)
	ks2_full(a_1, b_1, w4, 1); //  w4=A(-1)*B(-1) fixcarry
	ks2_full(ahalf, bhalf, w5, 0); //  w5=A(1/2)
	ks2_full(a_half, b_half, w6, 1); //  w6=A(-1/2) fixcarry? 但好像不fixcarry也对

	for (i = 0; i < w_len; i++) {
		//step10: w2=w2+w5
		w2[i] = w2[i] + w5[i];	//w2 <- w2+w5

		//step11: w6=w6-w5
		w6[i] = w6[i] - w5[i];	// w6 <- w6-w5

		//step12: w4=(w4-w3)/2
		w4[i] = w4[i] - w3[i];
		w4[i] = w4[i] / 2;

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
		w5[i] = w5[i] / 24;

		//step20: w6=w6+w2
		w6[i] = w2[i] + w6[i]; //w6 <- w6+w2

		//step21: w2=(w2+16*w4)/18
		w2[i] = w2[i] + (w4[i] << 4); //w2 <- w2+16*w4
		w2[i] = w2[i] / 18;

		//step22. w3 <- w3-w5
		w3[i] = w3[i] - w5[i];

		//step23: w4=-(w4+w2/2)
		w4[i] = -(w4[i] + w2[i]);

		//step24:  w6=(30*w2-w6)/60
		w6[i] = w2[i] * int30 - w6[i]; //w6 <- 30*w2-w6
		w6[i] = w6[i] / 60;

		//step25: w2=w2-w6
		w2[i] = w2[i] - w6[i];
	}

	int16_t res_full_len = 511; //2 * n - 1;
	int64_t res_full[res_full_len];
	memset(res_full, 0, sizeof(res_full[0]) * res_full_len);

	for (i = 0; i < w_len; i++) {
		res_full[i] += w7[i];
		res_full[sub_len + i] += w6[i];
		res_full[sub_len2 + i] += w5[i];
		res_full[sub_len3 + i] += w4[i];
		res_full[sub_len4 + i] += w3[i];
		res_full[sub_len5 + i] += w2[i];
		res_full[sub_len6 + i] += w1[i];
	}
	//---------------reduction-------
	for (i = 0; i < 255; i++) { //res[0 ~ (2n-2)-n=n-2=256-2=254]，剩下res[n-1]
		result[i] = toom4_general_reduce_int64_to_uint16(res_full[i] - res_full[i + n]);
	}
	result[n - 1] = toom4_general_reduce_int64_to_uint16(res_full[n - 1]);
}

//----speedtest-------
//a4,b4是w4的数据
int32_t ks2_a4[] = { 6538, 5954, 7038, 499, -2129, -5969, 8337, 5177, 1679, 2582, 2863, -709, 936, -599, 8023, -8800,
		3163,
		-5453, -1512, 2632, -3874, -4745, 2946, 217, 3316, 526, -2173, -3989, -254, 8613, 4507, 5513, -4345, 149,
		903, 6494, -8691, 4330, -1045, -213, -3202, -504, -2445, -9077, 4453, 3001, -5185, -766, 9791, 3006, 8351,
		-5465, 5473, 2450, 632, 1319, -5762, 2, -5359, -6913, -2083, -7212, -1606, 1262 };
int32_t ks2_b4[] = { -3, -7, 2, 5, 0, 4, -1, 2, -2, 4, -3, 0, 5, 1, 2, 2, -1, 1, -2, 1, -1, -3, -4, -1, 3, 7, 1, 0, -3,
		-1, -5, 1, -1, -1, 7, 1, 4, 5, -1, -1, -4, 3, -1, -3, 2, 0, -1, -1, 5, -2, 1, -5, -1, -6, 1, -1, 3, 1, 1, 1, 2,
		-2, 3, 1 };

//a2,b2是w2的数据
int32_t ks2_a2[] = { 28393, 25700, 63216, 88444, 82528, 54304, 47028, 44690, 20519, 83147, 24505, 79772, 32058, 64516,
		37918, 79826, 54733, 79825, 107427, 36490, 52931, 51274, 69708, 91039, 20824, 41659, 37058, 74476, 80923, 71562,
		52978, 9923, 79841, 52847, 89589, 31781, 62670, 20053, 60314, 80085, 62150, 44394, 43809, 64789, 53416, 28213,
		84926, 55571, 41219, 55821, 37442, 29896, 46075, 90896, 21278, 54248, 65614, 61520, 69692, 53495, 87824, 78744,
		83852, 73472 };
int32_t ks2_b2[] = { 115197, 115229, 115220, 115217, 115221, 115228, 115202, 115220, 115201, 115201, 115221, 115182,
		115208, 115219, 115208, 115214, 115208, 115225, 115219, 115219, 115223, 115233, 115205, 115214, 115197, 115207,
		115198, 115227, 115221, 115208, 115243, 115210, 115211, 115226, 115213, 115219, 115216, 115211, 115226, 115196,
		115238, 115221, 115214, 115224, 115205, 115215, 115232, 115214, 115205, 115201, 115225, 115219, 115193, 115203,
		115207, 115220, 115218, 115216, 115201, 115207, 115217, 115210, 115218, 115195 };

//a,b.len=64; res.len=127
void ks2_speedtest_sample(int32_t *a, int32_t *b, int64_t* res) {
	speedtest_print_title("ks2_speedtest_sample", 1);
	const int n = 64;	//用mpi_ks2代替kara，共64个数据
	int mpi_n = 48;	//每个实际需要18*2+6=42位，对齐是48位，ks2折半是24位,
	int mpi_n2 = mpi_n * 2;
	int16_t i;
	//24bit * 64 / 32bit = mpi中需要48个int32
	uint32_t borrow_fix_temp[n];	//借位

	REG_WRITE_NOCHECK(RSA_MULT_MODE_REG, 13);
	speedtest_reset_startcpucycles();

	//---------计算h0=f0*g0，两个正点值f0和g0的乘积--------------
	//写入正点值f0和g0
	toom4_general_fixborrow_f0g0(borrow_fix_temp, a, n);
	speedtest_print_cpucycles("fixborrow_f0g0");
	ks2_write_f0g0f1g1(borrow_fix_temp, ((uint32_t *) REG_BASE_X));
	speedtest_print_cpucycles("write_f0g0f1g1");
	toom4_general_fixborrow_f0g0(borrow_fix_temp, b, n);
	speedtest_print_cpucycles("fixborrow_f0g0");
	ks2_write_f0g0f1g1(borrow_fix_temp, ((uint32_t *) REG_BASE_Y) + 48);	//注意，从y的48字开始写数据(前48字填0)
	speedtest_print_cpucycles("write_f0g0f1g1");
	//reg写入0位
	//x低到高：48个X，48个0
	bzero(((uint32_t *) REG_BASE_X) + 48, 48 * 4);
	//y低到高：48个0，48个Y
	bzero(((uint32_t *) REG_BASE_Y), 48 * 4);
	speedtest_print_cpucycles("write_zero");
	//求正点值的乘积
	REG_WRITE_NOCHECK(RSA_MULT_START_REG, 1);	//start
	//speedtest_print_cpucycles("RSA_MULT_START_REG");
	while (REG_READ_NOCHECK(RSA_INTERRUPT_REG) != 1) {
	}
	speedtest_print_cpucycles("cpu_wait");
	//读出h0,h0=f0*g0
	for (i = 0; i < mpi_n2; i++) {	//注意h0为2倍的mpi_n=48*2
		h0.p[i] = REG_READ_NOCHECK(((uint32_t *) REG_BASE_Z) + i);
	}
	speedtest_print_cpucycles("read_h0");
	REG_WRITE_NOCHECK(RSA_INTERRUPT_REG, 1);
	speedtest_print_cpucycles("RSA_INTERRUPT_REG");
	//esp_mpi_release_hardware();//over

	//---------计算h1=f1*g1，两个负点的乘积--------------

	//写入负点值f1g1
	toom4_general_fixborrow_f1g1(borrow_fix_temp, a, n);
	speedtest_print_cpucycles("fixborrow_f1g1");
	ks2_write_f0g0f1g1(borrow_fix_temp, ((uint32_t *) REG_BASE_X));
	speedtest_print_cpucycles("write_f0g0f1g1");
	toom4_general_fixborrow_f1g1(borrow_fix_temp, b, n);
	speedtest_print_cpucycles("fixborrow_f1g1");
	ks2_write_f0g0f1g1(borrow_fix_temp, ((uint32_t *) REG_BASE_Y) + 48);	//注意，从y的48字开始写数据(前48字填0)
	speedtest_print_cpucycles("write_f0g0f1g1");

	//reg写入0位
	//x低到高：48个X，48个0
	bzero(((uint32_t *) REG_BASE_X) + 48, 48 * 4);
	//y低到高：48个0，48个Y
	bzero(((uint32_t *) REG_BASE_Y), 48 * 4);
	speedtest_print_cpucycles("write_zero");

	REG_WRITE_NOCHECK(RSA_MULT_START_REG, 1);	//start
	//speedtest_print_cpucycles("RSA_MULT_START_REG");
	while (REG_READ_NOCHECK(RSA_INTERRUPT_REG) != 1) {
	}
	speedtest_print_cpucycles("cpu_wait");

	//读出h1,h1=f1*g1
	for (i = 0; i < mpi_n2; i++) {	//注意h1为2倍的mpi_n
		h1.p[i] = REG_READ(((uint32_t *) REG_BASE_Z) + i);
	}
	speedtest_print_cpucycles("read_h1");
	REG_WRITE_NOCHECK(RSA_INTERRUPT_REG, 1);
	speedtest_print_cpucycles("RSA_INTERRUPT_REG");
	//esp_mpi_release_hardware();

	//-----------h0和h1计算完毕--------------------

	mpi_add_abs_lite(&h_add, &h0, &h1);
	speedtest_print_cpucycles("mpi_add_abs_lite");
	mpi_sub_abs_lite(&h_sub, &h0, &h1);
	speedtest_print_cpucycles("mpi_sub_abs_lite");

	mpi_shift_right_1bit_lite(&h_add);		//h0/2
	speedtest_print_cpucycles("mpi_shift_right_1bit_lite");
	mpi_shift_right_1bit_lite(&h_sub);		//h1/2
	speedtest_print_cpucycles("mpi_shift_right_1bit_lite");

	int32_t i_mpi = 0;
	int32_t i_poly = 0;
	//res中每个数为48位
	while (i_poly < 2 * n) {	//h_add，取出偶数项
		//i_mpi //0
		res[i_poly] = h_add.p[i_mpi + 1] & 0xFFFF;		//res[0]高16位
		res[i_poly] = res[i_poly] << 32;
		res[i_poly] |= h_add.p[i_mpi];		//res[0]低32位
		//--------------------------
		i_poly += 2;
		i_mpi += 1;	//1
		res[i_poly] = h_add.p[i_mpi + 1];	//res[2]高32位
		res[i_poly] = res[i_poly] << 16;
		res[i_poly] |= (h_add.p[i_mpi] >> 16) & 0xFFFF;	//res[2] 低16位
		i_mpi += 2;	//3
		i_poly += 2;
	}
	speedtest_print_cpucycles("unpack_h_add");
	i_mpi = 0;
	i_poly = 1;	//从第1项开始
	while (i_poly < 2 * n) {	//h_sub，取出奇数项
		//i_mpi // 0
		res[i_poly] = h_sub.p[i_mpi + 2] & 0xFF;		//res[1] high 8
		res[i_poly] = res[i_poly] << 40;
		res[i_poly] |= ((int64_t) h_sub.p[i_mpi + 1]) << 8;		// res[1] middle 32
		res[i_poly] |= (h_sub.p[i_mpi] >> 24) & 0xFF;	//res[1] low 8
		//-------------------------------
		i_mpi += 2;	//2
		i_poly += 2;	//3
		res[i_poly] = h_sub.p[i_mpi + 1] & 0xFFFFFF;	//res[3] high 24
		res[i_poly] = res[i_poly] << 24;
		res[i_poly] |= (h_sub.p[i_mpi] >> 8) & 0xFFFFFF;	//res[3] low 24
		i_mpi += 1;	//3
		i_poly += 2;	//5
	}
	speedtest_print_cpucycles("unpack_h_sub");

//	if (!fixcarry) {
//		return;
//	}
	//下面是处理进位的
	int carry0 = 0;		//偶数位进位
	int carry1 = 0;		//奇数位进位

	int64_t x = 281474976710656;		//x=2^48
	int64_t x_half = 140737488355328;		//x_half=x/2

	for (i = 0; i < (n * 2); i++) {
		//printf("fix_carry->%d at %d\n", i, n*2);
		if (i % 2 == 0) {		//h_add
			res[i] += carry0;
			//res[i] > x_half等价于carry0 = res[i] < 0 但不严谨
			if (res[i] > x_half) {		////carry1 = res[i] < 0;
				res[i] = res[i] - x;
				carry0 = 1;
			} else {
				carry0 = 0;
			}
		} else {  //h_sub
			res[i] += carry1;
			if (res[i] > x_half) {
				res[i] = res[i] - x;
				carry1 = 1;
			} else {
				carry1 = 0;
			}
		}
	}
	//卧槽，如果后面对res不使用的话，编译器直接就把fix_carry和unpack_h的两个循环直接去掉了
	speedtest_print_cpucycles("fix_carry");

	speedtest_print_splitline();
	ks2_full(a, b, res, 0);
	speedtest_print_cpucycles("ks2_full");
	ks2_full(a, b, res, 1);
	speedtest_print_cpucycles("ks2_full_fix_carry");
	speedtest_print_title("ks2_speedtest_sample", 0);
	//poly_print_int64(res, "res", 127, 0);
}
//ks2性能测试，含有两个sample
void ks2_speedtest() {
	int64_t res[127];
	ks2_speedtest_sample(ks2_a4, ks2_b4, res);
	ks2_speedtest_sample(ks2_a2, ks2_b2, res);
}

