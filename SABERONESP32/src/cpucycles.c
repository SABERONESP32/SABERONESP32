#include "cpucycles.h"
#include "stdint.h"


uint32_t cpucycles(void) {
	//Copy the codes in Esp.cpp for compatible with standard C language.
	uint32_t ccount;
	__asm__ __volatile__("esync; rsr %0,ccount":"=a" (ccount));
	return ccount;
	//#include "Esp.h"
	//return ESP.getCpuFreqMHz();
}
