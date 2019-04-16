#include "polymul_hw_common.h"

void schoolbook_mul(uint16_t* a, uint16_t* b, uint16_t* res, uint16_t q) {
	int n = 256;
	uint32_t i;
	uint32_t j, mask = 2 * n;
	//-------------------normal multiplication-----------------
	int64_t c[2 * n];
	for (i = 0; i < mask; i++) {
		c[i] = 0;
	}
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			c[i + j] += (int32_t) (a[i] * b[j]);
		}
	}
	//---------------reduction-------
	int64_t t;
	for (i = n; i < 2 * n; i++) {
		//c[i - n] = (c[i - n] - c[i]);	// % KYBER_Q;//& (p - 1);
		t = (c[i - n] - c[i]) % q;
		if (t < 0) {
			t += q;
		}
		res[i - n] = t;
	}
	//c[n] = 0; //256th coefficient=0;
	//----------------------copy result back----------------
	//for (i = 0; i < n; i++) {
	//res[i] = c[i];
	//}
}

static _lock_t rsa_accelerator_lock;

//---------------------
static inline uint16_t reduce_q_uint32_uint16(uint32_t v) {
#if REDUCE_Q == SABER_Q
	return v & SABER_MOD_MASK;
#else
	return (uint16_t) (v % REDUCE_Q);
#endif
}

void esp32_acquire_hardware(void) {
	/* newlib locks lazy initialize on ESP-IDF */
	_lock_acquire(&rsa_accelerator_lock);
	REG_SET_BIT(DPORT_PERI_CLK_EN_REG, DPORT_PERI_EN_RSA);
	/* also clear reset on digital signature, otherwise RSA is held in reset */
	REG_CLR_BIT(DPORT_PERI_RST_EN_REG, DPORT_PERI_EN_RSA | DPORT_PERI_EN_DIGITAL_SIGNATURE);

	REG_CLR_BIT(DPORT_RSA_PD_CTRL_REG, DPORT_RSA_PD);

	while (REG_READ(RSA_CLEAN_REG) != 1) {
	}

}
void esp32_release_hardware(void) {
	REG_SET_BIT(DPORT_RSA_PD_CTRL_REG, DPORT_RSA_PD);

	/* don't reset digital signature unit, as this resets AES also */
	REG_SET_BIT(DPORT_PERI_RST_EN_REG, DPORT_PERI_EN_RSA);
	REG_CLR_BIT(DPORT_PERI_CLK_EN_REG, DPORT_PERI_EN_RSA);

	_lock_release(&rsa_accelerator_lock);
}

//------------------------------

//a，b长度为64（每个占32位，32*64=2048），res长度为64*2-1=127,计算a*b=res，a和b必须都是正数
//传入的a，b必须为正数，且不多于13位
void IRAM_ATTR ks1_full_int32_int32(int32_t* a, int32_t* b, int32_t* res) {
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
		res[i] = t % REDUCE_Q;
	}
	REG_WRITE_NOCHECK(RSA_INTERRUPT_REG, 1);	//clear interrupt
}

void IRAM_ATTR ks1_full_uint16_int32(uint16_t *a, uint16_t *b, int32_t* res) {
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
		res[i] = t % REDUCE_Q;
	}
	REG_WRITE_NOCHECK(RSA_INTERRUPT_REG, 1);	//clear interrupt
}

void IRAM_ATTR ks1_full_uint16_uint16(uint16_t *a, uint16_t *b, uint16_t* res) {
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
		res[i] = reduce_q_uint32_uint16(t);
	}
	REG_WRITE_NOCHECK(RSA_INTERRUPT_REG, 1);	//clear interrupt
}

void IRAM_ATTR ks1_write_uint16(uint16_t* a, uint16_t* b) {
	int i = 0;
	//void *memcpy(void *dest, const void *src, size_t n);
	//由于a,b是int16，寄存器是32位的，二者不一样大，所以不能用内存拷贝
	//X低到高->64个X，64个0，
	//Y低到高->64个0，64个Y
	for (i = 0; i < 64; i++) {
		REG_WRITE_NOCHECK(((uint32_t *) REG_BASE_X) + i, (uint32_t ) a[i]);
		REG_WRITE_NOCHECK(((uint32_t *) REG_BASE_Y) + i + 64, (uint32_t ) b[i]);
	}
	bzero(((uint32_t *) REG_BASE_X) + 64, 256);		// 64 * 4);
	bzero(((uint32_t *) REG_BASE_Y), 256);		//64 * 4);	// 64个0
	REG_WRITE_NOCHECK(RSA_MULT_START_REG, 1);	//start
}

//传入的a，b必须为正数，且不多于13位
//传入的n位a,b的长度，输出的res长度为2n-1.这里的n=64
void IRAM_ATTR ks1_write_int32(int32_t* a, int32_t* b) {
	//void *memcpy(void *dest, const void *src, size_t n);//n为字节数
	//X低到高->64个X，64个0
	memcpy(((uint32_t *) REG_BASE_X), a, 64 * 4);
	bzero(((uint32_t *) REG_BASE_X) + 64, 64 * 4);
	//Y低到高->64个0，64个Y
	bzero(((uint32_t *) REG_BASE_Y), 64 * 4);	// 64个0
	memcpy(((uint32_t *) REG_BASE_Y) + 64, b, 64 * 4);
	REG_WRITE_NOCHECK(RSA_MULT_START_REG, 1);	//start
}

//从寄存器读出结果，长度为127，res的长度为2*64-1=127
void IRAM_ATTR ks1_read(uint32_t* ks1_res32) {
	while (REG_READ_NOCHECK(RSA_INTERRUPT_REG) != 1) {
	}
	memcpy(ks1_res32, ((uint32_t *) REG_BASE_Z), 127 * 4);	//pbase_z->ks1_res32
	REG_WRITE_NOCHECK(RSA_INTERRUPT_REG, 1);
}

void IRAM_ATTR ks1_saber_uint16_uint16(uint16_t *a, uint16_t *b, uint16_t* res, uint16_t mask) {
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
		res[i] = t & mask;
	}
	REG_WRITE_NOCHECK(RSA_INTERRUPT_REG, 1);	//clear interrupt
}

//------------------------------

void poly_print_uint16(uint16_t *p, char* name, int n, int full) {
	int i = 0;
#define FORMATE "%5u"
	printf("%20s-> { ", name);
	if (full) {
		for (i = 0; i < n; i++) {
			printf(FORMATE, p[i]);
			if (i < (n - 1)) {
				printf(",");
			}
		}
	} else {
		int max = POLY_PRINT_MAX;
		for (i = 0; i < max; i++) {
			printf(FORMATE, p[i]);
			printf(",");
		}
		printf(" ... ");
		for (i = n - max; i < n; i++) {
			printf(FORMATE, p[i]);
			if (i < (n - 1)) {
				printf(",");
			}
		}
	}
#undef FORMATE
	printf(" }\n");
}

void poly_print_int32(int32_t *p, char* name, int n, int full) {
	int i = 0;
#define FORMATE "%13d"
	printf("%16s-> { ", name);
	if (full) {
		for (i = 0; i < n; i++) {
			printf(FORMATE, p[i]);
			if (i < (n - 1)) {
				printf(",");
			}
		}
	} else {
		int max = POLY_PRINT_MAX;
		for (i = 0; i < max; i++) {
			printf(FORMATE, p[i]);
			printf(",");
		}
		printf(" ... ");
		for (i = n - max; i < n; i++) {
			printf(FORMATE, p[i]);
			if (i < (n - 1)) {
				printf(",");
			}
		}
	}
#undef FORMATE
	printf(" }\n");
}

void poly_print_int64(int64_t *p, char* name, int n, int full) {
	int i = 0;
#define FORMATE  "%13lld"//  I64d "%13lu"//
	printf("%16s-> { ", name);
	if (full) {
		for (i = 0; i < n; i++) {
			printf(FORMATE, p[i]);
			if (i < (n - 1)) {
				printf(",");
			}
		}
	} else {
		int max = POLY_PRINT_MAX;
		for (i = 0; i < max; i++) {
			printf(FORMATE, p[i]);
			printf(",");
		}
		printf(" ... ");
		for (i = n - max; i < n; i++) {
			printf(FORMATE, p[i]);
			if (i < (n - 1)) {
				printf(",");
			}
		}
	}
#undef FORMATE
	printf(" }\n");
}

void check_equal_int16(int16_t* a, int16_t* b, uint32_t n, char* name, int w_index) {
	int i;
	int equal = 1;
	for (i = 0; i < n; i++) {
		if (a[i] != b[i]) {
			equal = 0;
			//printf("not equal -> a[%d]=%d, b[%d]=%d\n", i, a[i], i,  b[i]);
			//break;
		}
	}
	printf("%10s-> [%d] ", name, w_index);
	if (equal) {
		printf("equal~~\n");
	} else {
		printf("NOT equal!!!!!\n");
	}
}

void check_equal_uint16(uint16_t* a, uint16_t* b, uint32_t n, char* name) {
	int i;
	int equal = 1;
	for (i = 0; i < n; i++) {
		if (a[i] != b[i]) {
			equal = 0;
			//printf("[%d]->  %5u vs. %5u \n", i, a[i], b[i]);
			//break;
		}
	}
	printf("%17s-> ", name);
	if (equal) {
		printf("equal~~\n");
	} else {
		printf("NOT equal!!!!!\n");
	}

}

void check_equal_int64(int64_t* a, int64_t* b, uint32_t n, char* name) {
	int i;
	int equal = 1;
	for (i = 0; i < n; i++) {
		if (a[i] != b[i]) {
			equal = 0;
			//printf("[%d]->  %5u vs. %5u \n", i, a[i], b[i]);
			break;
		}
	}
	printf("%17s-> ", name);
	if (equal) {
		printf("equal~~\n");
	} else {
		printf("NOT equal!!!!!\n");
	}

}

