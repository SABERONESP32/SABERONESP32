//
//  rng.c
//
//  Created by Bassham, Lawrence E (Fed) on 8/29/17.
//  Copyright © 2017 Bassham, Lawrence E (Fed). All rights reserved.
//  Modified by saber on esp32
//

#include <string.h>
#include "rng.h"
#include "esp_attr.h"
#include "esp_system.h"
#include "polymul_hw_speed.h"
//#include <openssl/conf.h>
//#include <openssl/evp.h>
//#include <openssl/err.h>

#define RNG_DATA_REG   0x3FF75144

#define REG_WRITE_NOCHECK(_r, _v) ({                                                                                           \
            (*(volatile uint32_t *)(_r)) = (_v);                                                                       \
        })

//read value from register
#define REG_READ_NOCHECK(_r) ({                                                                                                \
            (*(volatile uint32_t *)(_r));                                                                              \
        })

uint32_t inline hw_random(){
	return REG_READ_NOCHECK(RNG_DATA_REG);
}

int IRAM_ATTR randombytes(unsigned char *x,  unsigned long long xlen) {
	union {
		unsigned char aschar[4];
		uint32_t asint;
	} random;


	while (xlen > 0) {
		//speedtest_print_cpucycles("random");
		random.asint = hw_random();//esp_random();
		//speedtest_reset_startcpucycles();
		//the time cost by unpacking is longer than 48 cycles.
		*x++ = random.aschar[0];
		*x++ = random.aschar[1];
		*x++ = random.aschar[2];
		*x++ = random.aschar[3];
		xlen -= 4;
		//如果xlen是ull,数据处理需要50时钟周期
		//如果xlen是ul或uint整个数据处理需要43cpucycles,而我们需要48个时钟周期
	}

	return 0;
}






