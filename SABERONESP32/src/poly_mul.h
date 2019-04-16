#ifndef POLY_MUL_H
#define POLY_MUL_H

#include"SABER_params.h"

#ifdef __cplusplus
extern "C" {
#endif

void print_poly2(int16_t *a, int32_t n, uint64_t p);
void pol_mul(uint16_t* a, uint16_t* b, uint16_t* res, uint16_t p, uint16_t n);

void pol_mul_sb(int16_t* a, int16_t* b, int16_t* res, uint16_t p, uint32_t n,uint32_t start);
void toom_cook_4way(uint16_t* a1,uint16_t* b1, uint16_t* result, uint32_t p_mod, uint16_t n);




#ifdef __cplusplus
}
#endif

#endif
