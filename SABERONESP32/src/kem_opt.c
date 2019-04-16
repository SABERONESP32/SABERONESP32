/*
 * kem_opt.c
 *
 *  Created on: 2019年3月8日
 *      Author: Mai
 */

#include "kem_opt.h"
#include "fips202.h"

void byte_bank2pol_part(unsigned char *bytes, uint16_t pol_part[],
		uint16_t pol_part_start_index, uint16_t num_8coeff) {
	uint32_t j;
	uint32_t offset_data = 0, offset_byte = 0;

	offset_byte = 0;

	for (j = 0; j < num_8coeff; j++) {
		offset_byte = 13 * j;
		offset_data = pol_part_start_index + 8 * j;
		pol_part[offset_data + 0] = (bytes[offset_byte + 0] & (0xff)) | ((bytes[offset_byte + 1] & 0x1f) << 8);
		pol_part[offset_data + 1] = (bytes[offset_byte + 1] >> 5 & (0x07)) | ((bytes[offset_byte + 2] & 0xff) << 3)
				| ((bytes[offset_byte + 3] & 0x03) << 11);
		pol_part[offset_data + 2] = (bytes[offset_byte + 3] >> 2 & (0x3f)) | ((bytes[offset_byte + 4] & 0x7f) << 6);
		pol_part[offset_data + 3] = (bytes[offset_byte + 4] >> 7 & (0x01)) | ((bytes[offset_byte + 5] & 0xff) << 1)
				| ((bytes[offset_byte + 6] & 0x0f) << 9);
		pol_part[offset_data + 4] = (bytes[offset_byte + 6] >> 4 & (0x0f)) | ((bytes[offset_byte + 7] & 0xff) << 4)
				| ((bytes[offset_byte + 8] & 0x01) << 12);
		pol_part[offset_data + 5] = (bytes[offset_byte + 8] >> 1 & (0x7f)) | ((bytes[offset_byte + 9] & 0x3f) << 7);
		pol_part[offset_data + 6] = (bytes[offset_byte + 9] >> 6 & (0x03)) | ((bytes[offset_byte + 10] & 0xff) << 2)
				| ((bytes[offset_byte + 11] & 0x07) << 10);
		pol_part[offset_data + 7] = (bytes[offset_byte + 11] >> 3 & (0x1f)) | ((bytes[offset_byte + 12] & 0xff) << 5);
	}
}

void GenMatrix_opt_justintime(polyvec *a, const unsigned char *seed, uint16_t transpose) {
	int j;	//,k;
	unsigned char shake_op_buf[SHAKE128_RATE + 112];// there can be at most 112 bytes left over from previous shake call

	uint64_t s[25];
	memset(s, 0, sizeof(s[0]) * 25);

	uint16_t pol_part_start_index, num_8coeff, num_8coeff_final, left_over_bytes, total_bytes;
	uint16_t row, column, num_polynomial;

	keccak_absorb(s, SHAKE128_RATE, seed, SABER_SEEDBYTES, 0x1F);

	pol_part_start_index = 0;
	num_8coeff = 0;
	left_over_bytes = 0;
	total_bytes = 0;
	num_polynomial = 0;

	for (;;) {	//while (num_polynomial != 9)

		keccak_squeezeblocks(shake_op_buf + left_over_bytes, 1, s, SHAKE128_RATE);

		total_bytes = left_over_bytes + SHAKE128_RATE;

		num_8coeff = total_bytes / 13;

		if ((num_8coeff * 8 + pol_part_start_index) > 255) {
			num_8coeff_final = 32 - pol_part_start_index / 8;
		} else {
			num_8coeff_final = num_8coeff;
		}

		row = num_polynomial / 3;
		column = num_polynomial % 3;
		byte_bank2pol_part(shake_op_buf, a[row].vec[column].coeffs, pol_part_start_index, num_8coeff_final);

		pol_part_start_index = pol_part_start_index + num_8coeff_final * 8;	// this will be >256 when the polynomial is complete.
		if (pol_part_start_index > 255) {
			pol_part_start_index = 0;
			if (num_polynomial == 8) {	//from 0 to 8, in total 9;
				break;
			}
			num_polynomial++;
		}
		left_over_bytes = total_bytes - num_8coeff_final * 13;
		for (j = 0; j < left_over_bytes; j++) {	// bring the leftover in the begining of the buffer.
			shake_op_buf[j] = shake_op_buf[num_8coeff_final * 13 + j];
		}
	}
}

void dualcore_GenMatrix_opt_justintime(polyvec *a, const unsigned char *seed, uint16_t transpose) {
	int j;	//,k;
	unsigned char shake_op_buf[SHAKE128_RATE + 112];// there can be at most 112 bytes left over from previous shake call

	uint64_t s[25];
	memset(s, 0, sizeof(s[0]) * 25);

	uint16_t pol_part_start_index, num_8coeff, num_8coeff_final, left_over_bytes, total_bytes;
	uint16_t row, column, num_polynomial;

	keccak_absorb(s, SHAKE128_RATE, seed, SABER_SEEDBYTES, 0x1F);

	pol_part_start_index = 0;
	num_8coeff = 0;
	left_over_bytes = 0;
	total_bytes = 0;
	num_polynomial = 0;

	for (;;) {	//while (num_polynomial != 9)

		keccak_squeezeblocks(shake_op_buf + left_over_bytes, 1, s, SHAKE128_RATE);

		total_bytes = left_over_bytes + SHAKE128_RATE;

		num_8coeff = total_bytes / 13;

		if ((num_8coeff * 8 + pol_part_start_index) > 255) {
			num_8coeff_final = 32 - pol_part_start_index / 8;
		} else {
			num_8coeff_final = num_8coeff;
		}

		row = num_polynomial / 3;
		column = num_polynomial % 3;
		byte_bank2pol_part(shake_op_buf, a[row].vec[column].coeffs, pol_part_start_index, num_8coeff_final);

		pol_part_start_index = pol_part_start_index + num_8coeff_final * 8;	// this will be >256 when the polynomial is complete.
		if (pol_part_start_index > 255) {
			pol_part_start_index = 0;
			//poly finish
			//printf("poly finish - %d\n", num_polynomial);
			dualcoreGiveSemaphore(semaphore_ext);			//give one poly
			if (num_polynomial == 8) {	//from 0 to 8, in total 9;
				//printf("Gen A finish\n");
				break;
			}
			num_polynomial++;
		}
		left_over_bytes = total_bytes - num_8coeff_final * 13;
		for (j = 0; j < left_over_bytes; j++) {	// bring the leftover in the begining of the buffer.
			shake_op_buf[j] = shake_op_buf[num_8coeff_final * 13 + j];
		}
	}
}

//void MatrixVectorMulRouding_justintime(polyvec *a, polyvec *skpv, polyvec* res,
//		uint16_t mod, int16_t transpose) {
////	saber_karatsuba_parallel_justintime(a, b, root_res_full, acc_previous,
////			acc_previous_need_rounding);
//	memset(res->vec[0].coeffs, 0, sizeof(res->vec[0].coeffs[0]) * SABER_N);
//	memset(res->vec[1].coeffs, 0, sizeof(res->vec[0].coeffs[0]) * SABER_N);
//	memset(res->vec[2].coeffs, 0, sizeof(res->vec[0].coeffs[0]) * SABER_N);
//	const int root_res_full_len = 511;	// 511;
//	uint16_t root_res_full[root_res_full_len];
//	//karatsuba result is saved in root_res_full
//	//reduce root_res_full to acc and rouding acc in next step
//	//if there is no acc_previous, acc_previous is NULL
//	//Gen A is always in row order
//	if (transpose == 0) {
//
//		//res_index = a_row_index
//		//---- A row 0 ----
//		saber_karatsuba_parallel_justintime(a[0].vec[0].coeffs, skpv->vec[0].coeffs, root_res_full,
//		NULL, 0);
//		saber_karatsuba_parallel_justintime(a[0].vec[1].coeffs, skpv->vec[1].coeffs, root_res_full,
//				res->vec[0].coeffs, 0);
//		saber_karatsuba_parallel_justintime(a[0].vec[2].coeffs, skpv->vec[2].coeffs, root_res_full,
//				res->vec[0].coeffs, 0);
//
//		//---- A row 1 ----
//		saber_karatsuba_parallel_justintime(a[1].vec[0].coeffs, skpv->vec[0].coeffs, root_res_full,
//				res->vec[0].coeffs, 1);
//		saber_karatsuba_parallel_justintime(a[1].vec[1].coeffs, skpv->vec[1].coeffs, root_res_full,
//				res->vec[1].coeffs, 0);
//		saber_karatsuba_parallel_justintime(a[1].vec[2].coeffs, skpv->vec[2].coeffs, root_res_full,
//				res->vec[1].coeffs, 0);
//
//		//---- A row 2 ----
//		saber_karatsuba_parallel_justintime(a[2].vec[0].coeffs, skpv->vec[0].coeffs, root_res_full,
//				res->vec[1].coeffs, 1);
//		saber_karatsuba_parallel_justintime(a[2].vec[1].coeffs, skpv->vec[1].coeffs, root_res_full,
//				res->vec[2].coeffs, 0);
//		saber_karatsuba_parallel_justintime(a[2].vec[2].coeffs, skpv->vec[2].coeffs, root_res_full,
//				res->vec[2].coeffs, 0);
//		reduce_root_res_full_to_acc_and_rounding(root_res_full, res->vec[2].coeffs);
//
//	} else {
//		//res_index = a_column_index
//		saber_karatsuba_parallel_justintime(a[0].vec[0].coeffs, skpv->vec[0].coeffs, root_res_full,
//		NULL, 0);
//		saber_karatsuba_parallel_justintime(a[0].vec[1].coeffs, skpv->vec[0].coeffs, root_res_full,
//				res->vec[0].coeffs, 0);
//		saber_karatsuba_parallel_justintime(a[0].vec[2].coeffs, skpv->vec[0].coeffs, root_res_full,
//				res->vec[1].coeffs, 0);
//
//		saber_karatsuba_parallel_justintime(a[1].vec[0].coeffs, skpv->vec[1].coeffs, root_res_full,
//				res->vec[2].coeffs, 0);
//		saber_karatsuba_parallel_justintime(a[1].vec[1].coeffs, skpv->vec[1].coeffs, root_res_full,
//				res->vec[0].coeffs, 0);
//		saber_karatsuba_parallel_justintime(a[1].vec[2].coeffs, skpv->vec[1].coeffs, root_res_full,
//				res->vec[1].coeffs, 0);
//
//		saber_karatsuba_parallel_justintime(a[2].vec[0].coeffs, skpv->vec[2].coeffs, root_res_full,
//				res->vec[2].coeffs, 0);
//		saber_karatsuba_parallel_justintime(a[2].vec[1].coeffs, skpv->vec[2].coeffs, root_res_full,
//				res->vec[0].coeffs, 1);
//		saber_karatsuba_parallel_justintime(a[2].vec[2].coeffs, skpv->vec[2].coeffs, root_res_full,
//				res->vec[1].coeffs, 1);
//
//		reduce_root_res_full_to_acc_and_rounding(root_res_full, res->vec[2].coeffs);
//	}
//}

void dualcore_MatrixVectorMulRouding_justintime(polyvec *a, polyvec *skpv, polyvec* res,
		uint16_t mod, int16_t transpose) {
//	saber_karatsuba_parallel_justintime(a, b, root_res_full, acc_previous,
//			acc_previous_need_rounding);
	memset(res->vec[0].coeffs, 0, sizeof(res->vec[0].coeffs[0]) * SABER_N);
	memset(res->vec[1].coeffs, 0, sizeof(res->vec[0].coeffs[0]) * SABER_N);
	memset(res->vec[2].coeffs, 0, sizeof(res->vec[0].coeffs[0]) * SABER_N);
	const int root_res_full_len = 511;		// 511;
	uint16_t root_res_full[root_res_full_len];
	//karatsuba result is saved in root_res_full
	//reduce root_res_full to acc and rouding acc in next step
	//if there is no acc_previous, acc_previous is NULL
	//Gen A is always in row order
	if (transpose == 0) {

		//res_index = a_row_index
		//---- A row 0 ----
		dualcoreTakeSemaphore(semaphore_ext);
		saber_karatsuba_parallel_justintime(a[0].vec[0].coeffs, skpv->vec[0].coeffs, root_res_full,
		NULL, 0);
		dualcoreTakeSemaphore(semaphore_ext);
		saber_karatsuba_parallel_justintime(a[0].vec[1].coeffs, skpv->vec[1].coeffs, root_res_full,
				res->vec[0].coeffs, 0);
		dualcoreTakeSemaphore(semaphore_ext);
		saber_karatsuba_parallel_justintime(a[0].vec[2].coeffs, skpv->vec[2].coeffs, root_res_full,
				res->vec[0].coeffs, 0);

		//---- A row 1 ----
		dualcoreTakeSemaphore(semaphore_ext);
		saber_karatsuba_parallel_justintime(a[1].vec[0].coeffs, skpv->vec[0].coeffs, root_res_full,
				res->vec[0].coeffs, 1);
		dualcoreTakeSemaphore(semaphore_ext);
		saber_karatsuba_parallel_justintime(a[1].vec[1].coeffs, skpv->vec[1].coeffs, root_res_full,
				res->vec[1].coeffs, 0);
		dualcoreTakeSemaphore(semaphore_ext);
		saber_karatsuba_parallel_justintime(a[1].vec[2].coeffs, skpv->vec[2].coeffs, root_res_full,
				res->vec[1].coeffs, 0);

		//---- A row 2 ----
		dualcoreTakeSemaphore(semaphore_ext);
		saber_karatsuba_parallel_justintime(a[2].vec[0].coeffs, skpv->vec[0].coeffs, root_res_full,
				res->vec[1].coeffs, 1);
		dualcoreTakeSemaphore(semaphore_ext);
		saber_karatsuba_parallel_justintime(a[2].vec[1].coeffs, skpv->vec[1].coeffs, root_res_full,
				res->vec[2].coeffs, 0);
		dualcoreTakeSemaphore(semaphore_ext);
		saber_karatsuba_parallel_justintime(a[2].vec[2].coeffs, skpv->vec[2].coeffs, root_res_full,
				res->vec[2].coeffs, 0);
		reduce_root_res_full_to_acc_and_rounding(root_res_full, res->vec[2].coeffs);

	} else {
		//res_index = a_column_index
		dualcoreTakeSemaphore(semaphore_ext);
		//printf("poly mul start - (0,0)\n" );
		saber_karatsuba_parallel_justintime(a[0].vec[0].coeffs, skpv->vec[0].coeffs, root_res_full,
		NULL, 0);
		dualcoreTakeSemaphore(semaphore_ext);
		//printf("poly mul start - (0,1)\n" );
		saber_karatsuba_parallel_justintime(a[0].vec[1].coeffs, skpv->vec[0].coeffs, root_res_full,
				res->vec[0].coeffs, 0);
		//printf("poly mul start - (0,2)\n" );
		dualcoreTakeSemaphore(semaphore_ext);
		saber_karatsuba_parallel_justintime(a[0].vec[2].coeffs, skpv->vec[0].coeffs, root_res_full,
				res->vec[1].coeffs, 0);

		dualcoreTakeSemaphore(semaphore_ext);
		//printf("poly mul start - (1,0)\n" );
		saber_karatsuba_parallel_justintime(a[1].vec[0].coeffs, skpv->vec[1].coeffs, root_res_full,
				res->vec[2].coeffs, 0);
		dualcoreTakeSemaphore(semaphore_ext);
		//printf("poly mul start - (1,1)\n" );
		saber_karatsuba_parallel_justintime(a[1].vec[1].coeffs, skpv->vec[1].coeffs, root_res_full,
				res->vec[0].coeffs, 0);
		dualcoreTakeSemaphore(semaphore_ext);
		//printf("poly mul start - (1,2)\n" );
		saber_karatsuba_parallel_justintime(a[1].vec[2].coeffs, skpv->vec[1].coeffs, root_res_full,
				res->vec[1].coeffs, 0);

		dualcoreTakeSemaphore(semaphore_ext);
		//printf("poly mul start - (2,0)\n" );
		saber_karatsuba_parallel_justintime(a[2].vec[0].coeffs, skpv->vec[2].coeffs, root_res_full,
				res->vec[2].coeffs, 0);
		dualcoreTakeSemaphore(semaphore_ext);
		//printf("poly mul start - (2,1)\n" );
		saber_karatsuba_parallel_justintime(a[2].vec[1].coeffs, skpv->vec[2].coeffs, root_res_full,
				res->vec[0].coeffs, 1);
		dualcoreTakeSemaphore(semaphore_ext);
		//printf("poly mul start - (2,2)\n" );
		saber_karatsuba_parallel_justintime(a[2].vec[2].coeffs, skpv->vec[2].coeffs, root_res_full,
				res->vec[1].coeffs, 1);

		reduce_root_res_full_to_acc_and_rounding(root_res_full, res->vec[2].coeffs);
	}
}

void dualcore_GenMatrix_opt(polyvec *a, const unsigned char *seed, const int column_first) {
	unsigned int one_vector = 13 * SABER_N / 8;
	unsigned int byte_bank_length = SABER_K * SABER_K * one_vector;
	unsigned char buf[byte_bank_length];

	uint16_t temp_ar[SABER_N];

	int i, j, k;
	const uint16_t mod = SABER_Q_MASK; //(SABER_Q-1);

	shake128(buf, byte_bank_length, seed, SABER_SEEDBYTES);
	if (column_first == 0) {
		for (i = 0; i < SABER_K; i++) {
			for (j = 0; j < SABER_K; j++) {
				BS2POLq(buf + (i * SABER_K + j) * one_vector, temp_ar);
				for (k = 0; k < SABER_N; k++) {
					a[i].vec[j].coeffs[k] = (temp_ar[k]) & mod;
				}
			}
			dualcoreGiveSemaphore(semaphore_ext);
		}
	} else { //column_first == 1
		for (j = 0; j < SABER_K; j++) {
			for (i = 0; i < SABER_K; i++) {
				BS2POLq(buf + (i * SABER_K + j) * one_vector, temp_ar);
				for (k = 0; k < SABER_N; k++) {
					a[i].vec[j].coeffs[k] = (temp_ar[k]) & mod;
				}
			}
			dualcoreGiveSemaphore(semaphore_ext);
		}
	}
}

void GenMatrix_opt(polyvec *a, const unsigned char *seed, const int column_first) {
	unsigned int one_vector = 13 * SABER_N / 8;
	unsigned int byte_bank_length = SABER_K * SABER_K * one_vector;
	unsigned char buf[byte_bank_length];

	uint16_t temp_ar[SABER_N];

	int i, j, k;
	const uint16_t mod = SABER_Q_MASK; //(SABER_Q-1);
	shake128(buf, byte_bank_length, seed, SABER_SEEDBYTES);
	if (column_first == 0) {
		for (i = 0; i < SABER_K; i++) {
			for (j = 0; j < SABER_K; j++) {
				BS2POLq(buf + (i * SABER_K + j) * one_vector, temp_ar);
				for (k = 0; k < SABER_N; k++) {
					a[i].vec[j].coeffs[k] = (temp_ar[k]) & mod;
				}
			}
		}
	} else {
		for (j = 0; j < SABER_K; j++) {
			for (i = 0; i < SABER_K; i++) {
				BS2POLq(buf + (i * SABER_K + j) * one_vector, temp_ar);
				for (k = 0; k < SABER_N; k++) {
					a[i].vec[j].coeffs[k] = (temp_ar[k]) & mod;
				}
			}
		}
	}
}

void GenSecret_opt(polyvec* r, const unsigned char *seed) {
	//r[SABER_K][SABER_N]
	//要传polyvec的话，必须传地址才对，不传地址的话win运行正确，esp32运行错误
	uint16_t i;
	int32_t buf_size = SABER_MU * SABER_N * SABER_K / 8;
	uint8_t buf[buf_size];
	shake128(buf, buf_size, seed, SABER_NOISESEEDBYTES);
	for (i = 0; i < SABER_K; i++) {
		cbd((r->vec[i].coeffs), buf + i * SABER_MU * SABER_N / 8);
	}
}

//void GenSecret_opt(polyvec* r, const unsigned char *seed) {
//	//r[SABER_K][SABER_N]
//	//要传polyvec的话，必须传地址才对，不传地址的话win运行正确，esp32运行错误
//	uint16_t i;
//	int32_t buf_size = SABER_MU * SABER_N * SABER_K / 8;
//	uint8_t buf[buf_size];
//	shake128(buf, buf_size, seed, SABER_NOISESEEDBYTES);
//	for (i = 0; i < SABER_K; i++) {
//		cbd((r->vec[i].coeffs), buf + i * SABER_MU * SABER_N / 8);
//	}
//}

//message_dec_unpacked->message_dec
void POL2MSG_opt(uint16_t *message_dec_unpacked, unsigned char *message_dec) {
	int32_t i, j;
	for (j = 0; j < SABER_KEYBYTES; j++) {
		message_dec[j] = 0;
		for (i = 0; i < 8; i++)
			message_dec[j] = message_dec[j] | (message_dec_unpacked[j * 8 + i] << i);
	}
}

//void VectorMul_opt(polyvec* pkcl, polyvec* skpv,
//		uint16_t mod, uint16_t *res) {
//	saber_karatsuba_parallel_vector(pkcl->vec[0].coeffs, pkcl->vec[1].coeffs,
//			pkcl->vec[2].coeffs, skpv->vec[0].coeffs, skpv->vec[1].coeffs,
//			skpv->vec[2].coeffs, res);
//}

void VectorMul_opt(polyvec* pkcl, polyvec* skpv, uint16_t mod, uint16_t *res) {
//	toom4_saber_10bits_parallel_vector_speedtest(pkcl->vec[0].coeffs, pkcl->vec[1].coeffs,
//			pkcl->vec[2].coeffs, skpv->vec[0].coeffs, skpv->vec[1].coeffs,
//			skpv->vec[2].coeffs, res);
	toom4_saber_10bits_parallel_vector(pkcl->vec[0].coeffs, pkcl->vec[1].coeffs,
			pkcl->vec[2].coeffs, skpv->vec[0].coeffs, skpv->vec[1].coeffs,
			skpv->vec[2].coeffs, res);
	return;
//	uint16_t j, k;
//	uint16_t acc[SABER_N];
//	memset(res, 0, sizeof(res[0]) * SABER_N);
//
//	// vector-vector scalar multiplication with mod p
//	for (j = 0; j < SABER_K; j++) {
//		toom4_saber_10bits_parallel(pkcl->vec[j].coeffs, skpv->vec[j].coeffs, acc);
//		//toom4_saber_10bits(pkcl->vec[j].coeffs, skpv->vec[j].coeffs, acc);
//		//pol_mul(pkcl[j], skpv[j], acc, SABER_P, SABER_N, 0);
//		for (k = 0; k < SABER_N; k++) {
//			res[k] = res[k] + acc[k];
//			res[k] = res[k] & mod; //reduction
//			//acc[k] = 0; //clear the accumulator
//		}
//	}

}

void POLVECq2BS_opt(uint8_t *bytes, polyvec * data) {
	//uint16_t data[SABER_K][SABER_N]
	uint32_t i, j;
	uint32_t offset_data = 0, offset_byte = 0, offset_byte1 = 0;

	offset_byte = 0;
	for (i = 0; i < SABER_K; i++) {
		offset_byte1 = i * (SABER_N * 13) / 8;
		for (j = 0; j < SABER_N / 8; j++) {
			offset_byte = offset_byte1 + 13 * j;
			offset_data = 8 * j;
			bytes[offset_byte + 0] = (data->vec[i].coeffs[offset_data + 0] & (0xff));

			bytes[offset_byte + 1] = ((data->vec[i].coeffs[offset_data + 0] >> 8) & 0x1f)
					| ((data->vec[i].coeffs[offset_data + 1] & 0x07) << 5);

			bytes[offset_byte + 2] = ((data->vec[i].coeffs[offset_data + 1] >> 3) & 0xff);

			bytes[offset_byte + 3] = ((data->vec[i].coeffs[offset_data + 1] >> 11) & 0x03)
					| ((data->vec[i].coeffs[offset_data + 2] & 0x3f) << 2);

			bytes[offset_byte + 4] = ((data->vec[i].coeffs[offset_data + 2] >> 6) & 0x7f)
					| ((data->vec[i].coeffs[offset_data + 3] & 0x01) << 7);

			bytes[offset_byte + 5] = ((data->vec[i].coeffs[offset_data + 3] >> 1) & 0xff);

			bytes[offset_byte + 6] = ((data->vec[i].coeffs[offset_data + 3] >> 9) & 0x0f)
					| ((data->vec[i].coeffs[offset_data + 4] & 0x0f) << 4);

			bytes[offset_byte + 7] = ((data->vec[i].coeffs[offset_data + 4] >> 4) & 0xff);

			bytes[offset_byte + 8] = ((data->vec[i].coeffs[offset_data + 4] >> 12) & 0x01)
					| ((data->vec[i].coeffs[offset_data + 5] & 0x7f) << 1);

			bytes[offset_byte + 9] = ((data->vec[i].coeffs[offset_data + 5] >> 7) & 0x3f)
					| ((data->vec[i].coeffs[offset_data + 6] & 0x03) << 6);

			bytes[offset_byte + 10] = ((data->vec[i].coeffs[offset_data + 6] >> 2) & 0xff);

			bytes[offset_byte + 11] = ((data->vec[i].coeffs[offset_data + 6] >> 10) & 0x07)
					| ((data->vec[i].coeffs[offset_data + 7] & 0x1f) << 3);

			bytes[offset_byte + 12] = ((data->vec[i].coeffs[offset_data + 7] >> 5) & 0xff);

		}
	}
}

void POLVECp2BS_opt(uint8_t *bytes, polyvec* data) {
	//uint16_t data[SABER_K][SABER_N]
	uint32_t i, j;
	uint32_t offset_data = 0, offset_byte = 0, offset_byte1 = 0;

	offset_byte = 0;
	for (i = 0; i < SABER_K; i++) {
		offset_byte1 = i * (SABER_N * 10) / 8;
		for (j = 0; j < SABER_N / 4; j++) {
			offset_byte = offset_byte1 + 5 * j;
			offset_data = 4 * j;
			bytes[offset_byte + 0] = (data->vec[i].coeffs[offset_data + 0] & (0xff));

			bytes[offset_byte + 1] = ((data->vec[i].coeffs[offset_data + 0] >> 8) & 0x03)
					| ((data->vec[i].coeffs[offset_data + 1] & 0x3f) << 2);

			bytes[offset_byte + 2] = ((data->vec[i].coeffs[offset_data + 1] >> 6) & 0x0f)
					| ((data->vec[i].coeffs[offset_data + 2] & 0x0f) << 4);

			bytes[offset_byte + 3] = ((data->vec[i].coeffs[offset_data + 2] >> 4) & 0x3f)
					| ((data->vec[i].coeffs[offset_data + 3] & 0x03) << 6);

			bytes[offset_byte + 4] = ((data->vec[i].coeffs[offset_data + 3] >> 2) & 0xff);
		}
	}

}

void BS2POLVECp_opt(const unsigned char *bytes, polyvec *data) {

	uint32_t i, j;
	uint32_t offset_data = 0, offset_byte = 0, offset_byte1 = 0;

	offset_byte = 0;
	for (i = 0; i < SABER_K; i++) {
		offset_byte1 = i * (SABER_N * 10) / 8;
		for (j = 0; j < SABER_N / 4; j++) {
			offset_byte = offset_byte1 + 5 * j;
			offset_data = 4 * j;
			data->vec[i].coeffs[offset_data + 0] = (bytes[offset_byte + 0] & (0xff))
					| ((bytes[offset_byte + 1] & 0x03) << 8);
			data->vec[i].coeffs[offset_data + 1] = ((bytes[offset_byte + 1] >> 2) & (0x3f))
					| ((bytes[offset_byte + 2] & 0x0f) << 6);
			data->vec[i].coeffs[offset_data + 2] = ((bytes[offset_byte + 2] >> 4) & (0x0f))
					| ((bytes[offset_byte + 3] & 0x3f) << 4);
			data->vec[i].coeffs[offset_data + 3] = ((bytes[offset_byte + 3] >> 6) & (0x03))
					| ((bytes[offset_byte + 4] & 0xff) << 2);
		}
	}

}
//bytes->polyvec
void BS2POLVECq_opt(const unsigned char *bytes, polyvec *data) {

	uint32_t i, j;
	uint32_t offset_data = 0, offset_byte = 0, offset_byte1 = 0;

	offset_byte = 0;
	for (i = 0; i < SABER_K; i++) {
		offset_byte1 = i * (SABER_N * 13) / 8;
		for (j = 0; j < SABER_N / 8; j++) {
			offset_byte = offset_byte1 + 13 * j;
			offset_data = 8 * j;
			data->vec[i].coeffs[offset_data + 0] = (bytes[offset_byte + 0] & (0xff))
					| ((bytes[offset_byte + 1] & 0x1f) << 8);
			data->vec[i].coeffs[offset_data + 1] = (bytes[offset_byte + 1] >> 5 & (0x07))
					| ((bytes[offset_byte + 2] & 0xff) << 3) | ((bytes[offset_byte + 3] & 0x03) << 11);
			data->vec[i].coeffs[offset_data + 2] = (bytes[offset_byte + 3] >> 2 & (0x3f))
					| ((bytes[offset_byte + 4] & 0x7f) << 6);
			data->vec[i].coeffs[offset_data + 3] = (bytes[offset_byte + 4] >> 7 & (0x01))
					| ((bytes[offset_byte + 5] & 0xff) << 1) | ((bytes[offset_byte + 6] & 0x0f) << 9);
			data->vec[i].coeffs[offset_data + 4] = (bytes[offset_byte + 6] >> 4 & (0x0f))
					| ((bytes[offset_byte + 7] & 0xff) << 4) | ((bytes[offset_byte + 8] & 0x01) << 12);
			data->vec[i].coeffs[offset_data + 5] = (bytes[offset_byte + 8] >> 1 & (0x7f))
					| ((bytes[offset_byte + 9] & 0x3f) << 7);
			data->vec[i].coeffs[offset_data + 6] = (bytes[offset_byte + 9] >> 6 & (0x03))
					| ((bytes[offset_byte + 10] & 0xff) << 2) | ((bytes[offset_byte + 11] & 0x07) << 10);
			data->vec[i].coeffs[offset_data + 7] = (bytes[offset_byte + 11] >> 3 & (0x1f))
					| ((bytes[offset_byte + 12] & 0xff) << 5);
		}
	}

}

int crypto_kem_keypair_opt(unsigned char *pk, unsigned char *sk) {
	polyvec a[SABER_K];
	polyvec skpv; //[SABER_K][SABER_N];
	polyvec res; //uint16_t res[SABER_K][SABER_N]; res = A * skpv;
	unsigned char seed[SABER_SEEDBYTES]; //(32B)
	unsigned char noiseseed[SABER_COINBYTES]; //(32B)

	randombytes(seed, SABER_SEEDBYTES); //(random 32B)
	shake128(seed, SABER_SEEDBYTES, seed, SABER_SEEDBYTES); // for not revealing system RNG state
	randombytes(noiseseed, SABER_COINBYTES); //(random 32B)
	GenMatrix_opt(a, seed, 0); //
	GenSecret_opt(&skpv, noiseseed); //
	MatrixVectorMul_rounding_parallel(a, &skpv, &res, SABER_Q_MASK, 0); //matrix vector multiplication and rounding
	POLVECq2BS_opt(sk, &skpv); //

	POLVECp2BS_opt(pk, &res); //

	memcpy(pk + SABER_POLYVECCOMPRESSEDBYTES, seed, SABER_SEEDBYTES);
	memcpy(sk + SABER_INDCPA_SECRETKEYBYTES, pk, SABER_INDCPA_PUBLICKEYBYTES);
	sha3_256(sk + SABER_SECRETKEYBYTES - 64, pk, SABER_INDCPA_PUBLICKEYBYTES); // Then hash(pk) is appended.
	randombytes(sk + SABER_SECRETKEYBYTES - SABER_KEYBYTES, SABER_KEYBYTES); // Remaining part of sk contains a pseudo-random number.
	return (0);
}

int crypto_kem_enc_opt(unsigned char *c, unsigned char *k,
		const unsigned char *pk) {
	unsigned char kr[64];   // Will contain key, coins (k || random-noise)
	unsigned char buf[64];
	randombytes(buf, 32);
	sha3_256(buf, buf, 32); // BUF[0:31] <-- random message (will be used as the key for client) Note: hash doesnot release system RNG output
	sha3_256(buf + 32, pk, SABER_INDCPA_PUBLICKEYBYTES); // BUF[32:63] <-- Hash(public key);  Multitarget countermeasure for coins + contributory KEM
	sha3_512(kr, buf, 64);   // kr[0:63] <-- Hash(buf[0:63]);// K^ <-- kr[0:31]// noiseseed (r) <-- kr[32:63];
	//indcpa_kem_enc(buf, kr + 32, pk, c);// buf[0:31] contains message; kr[32:63] contains randomness r;
	//void indcpa_kem_enc(unsigned char *message_received, unsigned char *noiseseed,
	//	const unsigned char *pk, unsigned char *ciphertext)
	//========== INDCPA.ENC ===========
	unsigned char *message_received = buf;
	unsigned char *noiseseed = kr + 32;
	//unsigned char *ciphertext = c;
	uint32_t i, j;
	polyvec a[SABER_K];	// skpv;
	//unsigned char seed[SABER_SEEDBYTES];
	polyvec pkcl;	//uint16_t pkcl[SABER_K][SABER_N];		//public key of received by the client
	polyvec skpv1;	//uint16_t skpv1[SABER_K][SABER_N];
	uint16_t message[SABER_KEYBYTES * 8];
	polyvec res;	//uint16_t res[SABER_K][SABER_N];
	uint16_t mod_p = SABER_P - 1;
	uint16_t vprime[SABER_N];

	//unsigned char rec_c[SABER_RECONBYTES_KEM];

	//	for (i = 0; i < SABER_SEEDBYTES; i++) { // extract the seedbytes from Public Key.
	//		seed[i] = pk[ SABER_POLYVECCOMPRESSEDBYTES + i];
	//	}
	//memcpy(seed, pk + SABER_POLYVECCOMPRESSEDBYTES, SABER_SEEDBYTES);
	//GenMatrix_opt(a, seed, 1);
	GenMatrix_opt(a, pk + SABER_POLYVECCOMPRESSEDBYTES, 1);

	GenSecret_opt(&skpv1, noiseseed);	//generate secret from constant-time binomial distribution
	MatrixVectorMul_rounding_parallel(a, &skpv1, &res, SABER_Q_MASK, 1);

	//POLVECp2BS_opt(ciphertext, &res);
	POLVECp2BS_opt(c, &res);

	//------now calculate the v'
	BS2POLVECp_opt(pk, &pkcl);	//unpack the public_key, pkcl is the b in the protocol
	//skpv1[i][j] = skpv1[i][j] & (mod_p);//可以省略
	VectorMul_opt(&pkcl, &skpv1, mod_p, vprime);

	for (j = 0; j < SABER_KEYBYTES; j++) {	// unpack message_received;
		for (i = 0; i < 8; i++) {
			message[8 * j + i] = ((message_received[j] >> i) & 0x01);
		}
	}
	for (i = 0; i < SABER_N; i++) {
		message[i] = (message[i] << 9);	// message encoding
		vprime[i] = vprime[i] + message[i];
	}
	ReconDataGen(vprime, c + SABER_POLYVECCOMPRESSEDBYTES);
	//============== INDCPA.ENC END ================

	sha3_256(kr + 32, c, SABER_BYTES_CCA_DEC);	//1088 B
	sha3_256(k, kr, 64);
	return (0);
}

int crypto_kem_dec_opt(unsigned char *k, const unsigned char *c,
		const unsigned char *sk) {
	int i, fail;
	unsigned char cmp[SABER_BYTES_CCA_DEC];
	unsigned char buf[64];	//buf[0:31] is message; buf[32:63] is hash(pk)
	unsigned char kr[64];	//Will contain key, coins (kr=sha512(message||hash(pk)))
	//kr[0:31] is Kprime(first 32B of hash(mssage||pk) or last 32B random in sk)
	//kr[32:63] is noise(last 32B of hash(mssage||pk) for enc, hash(c) for create K)
	const unsigned char *pk = sk + SABER_INDCPA_SECRETKEYBYTES;

	//==== INDCPA.dec begin =====
	//indcpa_kem_dec(sk, c, buf);	// dec(sk, c) -> message_dec(in buf[0:31])
	//indcpa_kem_dec(*sk, *ciphertext, *message_dec)
	const unsigned char *ciphertext = c;
	polyvec sksv;	//uint16_t sksv[SABER_K][SABER_N]; //secret key of the server
	polyvec pksv;	//uint16_t pksv[SABER_K][SABER_N];
	//uint8_t recon_ar[SABER_RECONBYTES_KEM];
	uint16_t message_dec_unpacked[SABER_KEYBYTES * 8];	// one element containes on decrypted bit;
	uint16_t mod_p = SABER_P - 1;
	uint16_t v[SABER_N];

	BS2POLVECq_opt(sk, &sksv);	//sksv is the secret-key
	BS2POLVECp_opt(ciphertext, &pksv);	//pksv is the ciphertext
	//sksv[i][j] = sksv[i][j] & (mod_p);//省略
	VectorMul_opt(&pksv, &sksv, mod_p, v);	//vector-vector scalar multiplication with mod p

	//recon_ar = ciphertext + SABER_POLYVECCOMPRESSEDBYTES
	Recon(v, ciphertext + SABER_POLYVECCOMPRESSEDBYTES, message_dec_unpacked);
	POL2MSG_opt(message_dec_unpacked, buf);	//pack(message_dec)->buf[0:31]
	//=== INDCPA.dec end =====

	// hash(pk) is in sk last[-64:-32]
	memcpy(buf + 32, sk + SABER_SECRETKEYBYTES - 64, 32);	//hash(pk) -> buf[32:63];
	sha3_512(kr, buf, 64);	//sha512(buf[0:63])->kr[0:63]; kr=sha512(message||pk)

	//(buf[0:31] is message, kr+32 is noise, pk) -> cmp (cmp is the ciphertext to be compared)

	//========== INDCPA.enc begin ======
	//indcpa_kem_enc(buf, kr + 32, pk, cmp);
	//indcpa_kem_enc(*message_received, *noiseseed, *pk, *ciphertext)
	//unsigned char *message_received = buf;//buf[0:31]
	unsigned char *noiseseed = kr + 32;
	//unsigned char *ciphertext = c;
	polyvec a[SABER_K];	// skpv;
	//unsigned char seed[SABER_SEEDBYTES];
	polyvec pkcl;	//uint16_t pkcl[SABER_K][SABER_N];		//public key of received by the client
	polyvec skpv1;	//uint16_t skpv1[SABER_K][SABER_N];
	//uint16_t message[SABER_KEYBYTES * 8];// message is just message_dec_unpacked
	uint16_t *message = message_dec_unpacked;
	polyvec res;	//uint16_t res[SABER_K][SABER_N];
	//uint16_t mod_p = SABER_P - 1;
	uint16_t vprime[SABER_N];

	//unsigned char rec_c[SABER_RECONBYTES_KEM];

	//seed = pk + SABER_POLYVECCOMPRESSEDBYTES
	GenMatrix_opt(a, pk + SABER_POLYVECCOMPRESSEDBYTES, 1);

	GenSecret_opt(&skpv1, noiseseed);	//generate secret from constant-time binomial distribution

	MatrixVectorMul_rounding_parallel(a, &skpv1, &res, SABER_Q_MASK, 1);

	POLVECp2BS_opt(cmp, &res);

	//------now calculate the v'
	BS2POLVECp_opt(pk, &pkcl);	//unpack the public_key, pkcl is the b in the protocol
	//skpv1[i][j] = skpv1[i][j] & (mod_p);//可以省略
	VectorMul_opt(&pkcl, &skpv1, mod_p, vprime);

//	for (j = 0; j < SABER_KEYBYTES; j++) {	// unpack message_received;
//		for (i = 0; i < 8; i++) {
//			message[8 * j + i] = ((message_received[j] >> i) & 0x01);
//		}
//	}// message is just message_dec_unpacked, no need to unpack

	for (i = 0; i < SABER_N; i++) {
		message[i] = (message[i] << 9);	// message encoding
		vprime[i] = vprime[i] + message[i];
	}

	ReconDataGen(vprime, cmp + SABER_POLYVECCOMPRESSEDBYTES);
	//============== INDCPA.ENC END ================

	//if c==cmp then k = hash(kr_r,c) else k = hash(z,c)  (z is in sk last 32B)
	fail = verify(c, cmp, SABER_BYTES_CCA_DEC);	//0 for equal strings, 1 for non-equal strings
	//fail = memcmp(c, cmp, SABER_BYTES_CCA_DEC);//a little slower than verify

	sha3_256(kr + 32, c, SABER_BYTES_CCA_DEC);	// overwrite coins in kr with h(c)

	//cmov(kr, sk + SABER_SECRETKEYBYTES - SABER_KEYBYTES, SABER_KEYBYTES, fail);
	if (fail == 0) {   //success
	} else {
		//if not success, random in sk last 32B -> kr[0:31]
		memcpy(kr, sk + SABER_SECRETKEYBYTES - SABER_KEYBYTES, SABER_KEYBYTES);
	}

	sha3_256(k, kr, 64);    // hash concatenation of pre-k and h(c) to k
	return (0);
}

//dualcore: run on EXT_CORE
void GenMatrix_opt_dualcore(polyvec *a, const unsigned char *seed, const int column_first) {
	unsigned int one_vector = 13 * SABER_N / 8;
	unsigned int byte_bank_length = SABER_K * SABER_K * one_vector;
	unsigned char buf[byte_bank_length];

	uint16_t temp_ar[SABER_N];

	int i, j, k;
	const uint16_t mod = SABER_Q_MASK; //(SABER_Q-1);

	shake128(buf, byte_bank_length, seed, SABER_SEEDBYTES);
	if (column_first == 0) {
		for (i = 0; i < SABER_K; i++) {
#if DUALCORE_DEBUG
			speedtest_reset_startcpucycles();
#endif
			for (j = 0; j < SABER_K; j++) {
				BS2POLq(buf + (i * SABER_K + j) * one_vector, temp_ar);
				for (k = 0; k < SABER_N; k++) {
					a[i].vec[j].coeffs[k] = (temp_ar[k]) & mod;
				}
			}
#if DUALCORE_DEBUG
			speedtest_print_cpucycles("Gen A[i]");
			dualcoreDebug("ext", "Gen A[i] over");
#endif
			dualcoreGiveSemaphore(semaphore_ext);
		}
	} else { //column_first == 1
		for (j = 0; j < SABER_K; j++) {
			for (i = 0; i < SABER_K; i++) {
				BS2POLq(buf + (i * SABER_K + j) * one_vector, temp_ar);
				for (k = 0; k < SABER_N; k++) {
					a[i].vec[j].coeffs[k] = (temp_ar[k]) & mod;
				}
			}
			//vTaskDelay(1000);
#if DUALCORE_DEBUG
			speedtest_print_cpucycles("Gen A[i]");
			dualcoreDebug("ext", "Gen A[i] over");
#endif
			dualcoreGiveSemaphore(semaphore_ext);
			//printf("A[i] ok, column_first\n");
		}
//		printf("A[0,1,2] ok, column_first\n");
		//poly_print_uint16(a[0].vec[0].coeffs, "a_test in ext", SABER_N, 0);

	}
}

//===========================
// dual core
//===========================


//semaphore_main is the semaphore for "main core"
//semaphore_ext is the semaphore for "secondary core"

//====== CCA.GEN ======

//keypair running on "secondary core"
void dualcore_crypto_kem_keypair_opt_ext_task(void * param) {
	PARAMS_5_T p = *((PARAMS_5_T*) param);//param pointers from main core
	unsigned char *sk = (unsigned char*) p.p0;
	unsigned char *pk = (unsigned char*) p.p1;
	unsigned char *seed = (unsigned char*) p.p2;
	polyvec *skpv = (polyvec*) p.p4;//vector sk
	polyvec *a = (polyvec*) p.p4; //matrix A

	shake128(seed, SABER_SEEDBYTES, seed, SABER_SEEDBYTES);
	//GenMatrix will call "GiveSemaphore" internally, here use the "justintime" version
	dualcore_GenMatrix_opt_justintime(a, seed, 0);
	dualcoreDeleteTask();
}

//keypair running on "main core"
int dualcore_crypto_kem_keypair_opt(unsigned char *pk, unsigned char *sk) {
	polyvec a[SABER_K];
	polyvec skpv;
	polyvec res; //res = A * skpv;
	unsigned char seed[SABER_SEEDBYTES]; //(32B)
	unsigned char noiseseed[SABER_COINBYTES]; //(32B)

	randombytes(seed, SABER_SEEDBYTES); //(random 32B)

	PARAMS_5_T p = { sk, pk, seed, &skpv, a };
	//start the task running on secondary core
	dualcoreCreateTask(dualcore_crypto_kem_keypair_opt_ext_task, &p);

	randombytes(noiseseed, SABER_COINBYTES); //(random 32B)
	GenSecret_opt(&skpv, noiseseed);

	//MatrixVectorMul_rounding will call "TakeSemaphore" internally
	dualcore_MatrixVectorMulRouding_justintime(a, &skpv, &res, SABER_Q_MASK, 0);//matrix vector multiplication and rounding

	POLVECq2BS_opt(sk, &skpv);

	POLVECp2BS_opt(pk, &res);

	memcpy(pk + SABER_POLYVECCOMPRESSEDBYTES, seed, SABER_SEEDBYTES);
	memcpy(sk + SABER_INDCPA_SECRETKEYBYTES, pk, SABER_INDCPA_PUBLICKEYBYTES);
	sha3_256(sk + SABER_SECRETKEYBYTES - 64, pk, SABER_INDCPA_PUBLICKEYBYTES); // Then hash(pk) is appended.
	randombytes(sk + SABER_SECRETKEYBYTES - SABER_KEYBYTES, SABER_KEYBYTES); // Remaining part of sk contains a pseudo-random number.
	return (0);
}

//====== CCA.ENC ======

//enc running on "secondary core"
void dualcore_crypto_kem_enc_opt_ext_task(void * param) {

	PARAMS_7_T p = *((PARAMS_7_T*) param);
	unsigned char *buf = (unsigned char*) p.p0;
	unsigned char *pk = (unsigned char*) p.p1;
	polyvec *pkcl = (polyvec*) p.p2;
	unsigned char *c = (unsigned char*) p.p3;
	uint16_t *vprime = (uint16_t*) p.p4;
	polyvec *a = (polyvec*) p.p5;
	uint16_t *message = (uint16_t*) p.p6;
	unsigned char *message_received = buf;

	randombytes(buf, 32);
	sha3_256(buf, buf, 32);	// BUF[0:31] <-- random message (used as the key for client)
	dualcoreGiveSemaphore(semaphore_ext);	//give buf[0:31]

	BS2POLVECp_opt(pk, pkcl);	//unpack the public_key

	dualcoreGiveSemaphore(semaphore_ext);	//give pkcl(pk_polyvec)

	//GenMatrix will call "GiveSemaphore" internally
	dualcore_GenMatrix_opt(a, pk + SABER_POLYVECCOMPRESSEDBYTES, 1);	//gen(seed in pk)->A

	uint16_t i, j;
	for (j = 0; j < SABER_KEYBYTES; j++) {	// unpack message_received;
		for (i = 0; i < 8; i++) {
			message[8 * j + i] = ((message_received[j] >> i) & 0x01);
		}
	}
	for (i = 0; i < SABER_N; i++) {
		message[i] = (message[i] << 9);	// message encoding
	}
	dualcoreTakeSemaphore(semaphore_main);	//take vprime

	for (i = 0; i < SABER_N; i++) {
		vprime[i] = vprime[i] + message[i];
	}
	ReconDataGen(vprime, c + SABER_POLYVECCOMPRESSEDBYTES);// ReconDataGen

	dualcoreGiveSemaphore(semaphore_ext);	//give ReconData (stored in c)

	dualcoreDeleteTask();

}

//enc running on "main core"
int dualcore_crypto_kem_enc_opt(unsigned char *c, unsigned char *k,
		const unsigned char *pk) {
	unsigned char kr[64];	// Will contain key, coins (k || random-noise)
	unsigned char buf[64];
	polyvec a[SABER_K];	// skpv;
	polyvec pkcl;	//public key of received by the client
	polyvec skpv1;
	uint16_t message[SABER_KEYBYTES * 8];
	polyvec res;
	uint16_t mod_p = SABER_P - 1;
	uint16_t vprime[SABER_N];
	unsigned char *noiseseed = kr + 32;

	//start the task running on secondary core
	PARAMS_7_T params = { buf, pk, &pkcl, c, vprime, a, message };
	dualcoreCreateTask(dualcore_crypto_kem_enc_opt_ext_task, &params);
	sha3_256(buf + 32, pk, SABER_INDCPA_PUBLICKEYBYTES); // BUF[32:63] <-- Hash(public key);  Multitarget countermeasure for coins + contributory KEM

	dualcoreTakeSemaphore(semaphore_ext); //take buf[0:31]
	sha3_512(kr, buf, 64); // kr[0:63] <-- Hash(buf[0:63]);// K^ <-- kr[0:31]// noiseseed (r) <-- kr[32:63];

	GenSecret_opt(&skpv1, noiseseed); //generate secret from constant-time binomial distribution

	dualcoreTakeSemaphore(semaphore_ext); //take pkcl(pk_polyvec)
	VectorMul_opt(&pkcl, &skpv1, mod_p, vprime);
	dualcoreGiveSemaphore(semaphore_main); //give vprime
	//MatrixVectorMul_rounding will call "TakeSemaphore" internally
	dualcore_MatrixVectorMul_rounding_parallel(a, &skpv1, &res, SABER_Q_MASK, 1);
	POLVECp2BS_opt(c, &res);

	dualcoreTakeSemaphore(semaphore_ext); //take ReconData (stored in c)

	sha3_256(kr + 32, c, SABER_BYTES_CCA_DEC); //1088B
	sha3_256(k, kr, 64);

	return (0);
}

//====== CCA.DEC ======

//dec running on "secondary core"
void dualcore_crypto_kem_dec_opt_ext_task(void * param) {
	PARAMS_9_T p = *((PARAMS_9_T*) param);
	const unsigned char *c = (unsigned char*) p.p0;
	const unsigned char *sk = (unsigned char*) p.p1;
	uint16_t *vprime = (uint16_t*) p.p2;
	uint16_t *message = (uint16_t*) p.p3;
	unsigned char *cmp = (unsigned char*) p.p4;
	polyvec *a = (polyvec*) p.p5;
	polyvec *pksv = (polyvec*) p.p6;
	polyvec *pkcl = (polyvec*) p.p7;
	unsigned char *kr = (unsigned char*) p.p8;

	const unsigned char *ciphertext = c;
	const unsigned char *pk = sk + SABER_INDCPA_SECRETKEYBYTES;

	BS2POLVECp_opt(ciphertext, pksv); //pksv is the ciphertext

	dualcoreGiveSemaphore(semaphore_ext); //give pksv

	BS2POLVECp_opt(pk, pkcl);	//unpack the public_key, pkcl is the b in the protocol
	dualcoreGiveSemaphore(semaphore_ext);	//give pkcl

	//GenMatrix will call "GiveSemaphore" internally
	dualcore_GenMatrix_opt(a, pk + SABER_POLYVECCOMPRESSEDBYTES, 1);	//give semaphore already in GenMatrix

	//no need to unpack message from buff, we already have it.
	uint16_t i;

	dualcoreTakeSemaphore(semaphore_main);	//take vprime
	for (i = 0; i < SABER_N; i++) {
		message[i] = (message[i] << 9);	// message encoding
		vprime[i] = vprime[i] + message[i];
	}
	ReconDataGen(vprime, cmp + SABER_POLYVECCOMPRESSEDBYTES);
	dualcoreGiveSemaphore(semaphore_ext);	//give cmp(with recon_data from vprime)
	//c=((3 * 320) + ((3+1)*256/8))=1088B
	sha3_256(kr + 32, c, SABER_BYTES_CCA_DEC);  // overwrite coins in kr with h(c)
	dualcoreGiveSemaphore(semaphore_ext); //give kr(with sha256(c) in kr[32:63])
	dualcoreDeleteTask();
}

//dec running on "main core"
int dualcore_crypto_kem_dec_opt(unsigned char *k, const unsigned char *c,
		const unsigned char *sk) {
	int fail;
	unsigned char cmp[SABER_BYTES_CCA_DEC];
	unsigned char buf[64];	//buf[0:31]->message; buf[32:43]->hash(pk)
	unsigned char kr[64];	// Will contain key, coins (kr=sha512(message||hash(pk)))
	//kr[0:31] is Kprime(first 32B of hash(mssage||pk) or last 32B random in sk)
	//kr[32:63] is noise(last 32B of hash(mssage||pk) for enc, hash(c) for create K)
	//const unsigned char *pk = sk + SABER_INDCPA_SECRETKEYBYTES;

	const unsigned char *ciphertext = c;
	polyvec sksv;	//secret key of the server
	polyvec pksv;
	uint16_t message_dec_unpacked[SABER_KEYBYTES * 8]; // one element containes on decrypted bit;
	uint16_t mod_p = SABER_P - 1;
	uint16_t v[SABER_N];


	unsigned char *noiseseed = kr + 32;
	polyvec a[SABER_K];

	polyvec pkcl;	//public key of received by the client
	polyvec skpv1;

	polyvec res;

	uint16_t vprime[SABER_N];
	//message is just message_dec_unpacked
	PARAMS_9_T params = { c, sk, vprime, message_dec_unpacked, cmp, a, &pksv, &pkcl, kr };
	dualcoreCreateTask(dualcore_crypto_kem_dec_opt_ext_task, &params);

	BS2POLVECq_opt(sk, &sksv); //sksv is the secret-key
	dualcoreTakeSemaphore(semaphore_ext); //take pksv

	VectorMul_opt(&pksv, &sksv, mod_p, v); //vector-vector scalar multiplication with mod p
	//recon_ar = ciphertext + SABER_POLYVECCOMPRESSEDBYTES
	Recon(v, ciphertext + SABER_POLYVECCOMPRESSEDBYTES, message_dec_unpacked);
	POL2MSG_opt(message_dec_unpacked, buf);	//pack(message_dec)->msg in buf[0:31]


	// hash(pk) is in sk last[-64:-32]
	memcpy(buf + 32, sk + SABER_SECRETKEYBYTES - 64, 32);	//hash(pk) -> buf[32:63];
	sha3_512(kr, buf, 64);	//sha512(buf[0:63])->kr[0:63]; kr=sha512(message||pk)
	//(buf[0:31] is message, kr+32 is noise, pk) -> cmp (cmp is the ciphertext to be compared)

	GenSecret_opt(&skpv1, noiseseed);	//generate secret from constant-time binomial distribution

	dualcoreTakeSemaphore(semaphore_ext);	//take pkcl
	//skpv1[i][j] = skpv1[i][j] & (mod_p);//可以省略
	VectorMul_opt(&pkcl, &skpv1, mod_p, vprime);
	dualcoreGiveSemaphore(semaphore_main);	//give vprime
	//MatrixVectorMul_rounding will call "TakeSemaphore" internally
	dualcore_MatrixVectorMul_rounding_parallel(a, &skpv1, &res, SABER_Q_MASK, 1);

	POLVECp2BS_opt(cmp, &res);

	dualcoreTakeSemaphore(semaphore_ext);	//take cmp(with recon_data from vprime)

	//if c==cmp then k = hash(kr_r,c) else k = hash(z,c)  (z is in sk last 32B)
	fail = verify(c, cmp, SABER_BYTES_CCA_DEC);//0 for equal strings, 1 for non-equal strings

	if (fail == 0) {   //success
	} else {   //if not success, copy the random in sk last 32B -> kr[0:31]
		memcpy(kr, sk + SABER_SECRETKEYBYTES - SABER_KEYBYTES, SABER_KEYBYTES);
	}

	dualcoreTakeSemaphore(semaphore_ext);   //take kr(with sha256(c) in kr[32:63])
	sha3_256(k, kr, 64);    // hash concatenation of pre-k and h(c) to k
	return (0);
}


//==========================================
// SpeedTest
//==========================================
void GenMatrix_opt_speedtest(polyvec *a, const unsigned char *seed, const int column_first) {
	unsigned int one_vector = 13 * SABER_N / 8; //416
	unsigned int byte_bank_length = SABER_K * SABER_K * one_vector;//3*3*416=3744
	unsigned char buf[byte_bank_length];

	uint16_t temp_ar[SABER_N];

	int i, j, k;
	const uint16_t mod = SABER_Q_MASK;//(SABER_Q-1);
	speedtest_print_cpucycles("Gen(A) prepare");
	shake128(buf, byte_bank_length, seed, SABER_SEEDBYTES);
	speedtest_print_cpucycles("Gen(A) shake128(32->3744)");
	if (column_first == 0) {
		for (i = 0; i < SABER_K; i++) {
			for (j = 0; j < SABER_K; j++) {
				BS2POLq(buf + (i * SABER_K + j) * one_vector, temp_ar);
				for (k = 0; k < SABER_N; k++) {
					a[i].vec[j].coeffs[k] = (temp_ar[k]) & mod;
				}
			}
			speedtest_print_cpucycles("Gen(A) A[i]");
		}
	} else {
		for (j = 0; j < SABER_K; j++) {
			for (i = 0; i < SABER_K; i++) {
				BS2POLq(buf + (i * SABER_K + j) * one_vector, temp_ar);
				for (k = 0; k < SABER_N; k++) {
					a[i].vec[j].coeffs[k] = (temp_ar[k]) & mod;
				}
			}
			speedtest_print_cpucycles("Gen(A) A[i]");
		}
	}
}

void GenSecret_opt_speedtest(polyvec* r, const unsigned char *seed) {
	//r[SABER_K][SABER_N]
	//卧槽，要传polyvec的话，必须传地址才对，不传地址的话win运行正确，esp32运行错误
	uint16_t i;
	int32_t buf_size = SABER_MU * SABER_N * SABER_K / 8;//8*256*3/8=768
	uint8_t buf[buf_size];
	speedtest_print_cpucycles("Gen(s) prepare");
	shake128(buf, buf_size, seed, SABER_NOISESEEDBYTES);
	speedtest_print_cpucycles("Gen(s) shake128(32B->768)");
	for (i = 0; i < SABER_K; i++) {
		cbd((r->vec[i].coeffs), buf + i * SABER_MU * SABER_N / 8);
		speedtest_print_cpucycles("Gen(s) s[i]");
	}

}
//
int crypto_kem_keypair_opt_speedtest(unsigned char *pk, unsigned char *sk) {
	speedtest_print_title("crypto_kem_keypair_opt", 1);
	polyvec a[SABER_K];
	polyvec skpv; //[SABER_K][SABER_N];
	polyvec res;//uint16_t res[SABER_K][SABER_N]; res = A * skpv;
	unsigned char seed[SABER_SEEDBYTES];//(32B)
	unsigned char noiseseed[SABER_COINBYTES];//(32B)
	speedtest_reset_startcpucycles();
	randombytes(seed, SABER_SEEDBYTES);//(random 32B)
	speedtest_print_cpucycles("rand(seed)");

	shake128(seed, SABER_SEEDBYTES, seed, SABER_SEEDBYTES);// for not revealing system RNG state
	speedtest_print_cpucycles("shake128(seed)");

	randombytes(noiseseed, SABER_COINBYTES);//(random 32B)
	speedtest_print_cpucycles("random(noise)");

	GenMatrix_opt(a, seed, 0);//
	speedtest_print_cpucycles("GenMatrix(seed)->A");

	GenSecret_opt(&skpv, noiseseed);//
	speedtest_print_cpucycles("GenSecret(noise)->sk");

	MatrixVectorMul_rounding_parallel(a, &skpv, &res, SABER_Q_MASK, 0);//matrix vector multiplication and rounding
	speedtest_print_cpucycles("MulRounding(A,sk)->res");

	POLVECq2BS_opt(sk, &skpv);//
	speedtest_print_cpucycles("Poly2Bs(sk)->sk");

	POLVECp2BS_opt(pk, &res);//
	speedtest_print_cpucycles("Poly22Bs(res)->pk");

	memcpy(pk + SABER_POLYVECCOMPRESSEDBYTES, seed, SABER_SEEDBYTES);
	speedtest_print_cpucycles("Insert(seed)->pk");

	memcpy(sk + SABER_INDCPA_SECRETKEYBYTES, pk, SABER_INDCPA_PUBLICKEYBYTES);
	speedtest_print_cpucycles("Insert(pk)->sk");

	sha3_256(sk + SABER_SECRETKEYBYTES - 64, pk, SABER_INDCPA_PUBLICKEYBYTES);// Then hash(pk) is appended.
	speedtest_print_cpucycles("Insert(sha256(pk))->sk");

	randombytes(sk + SABER_SECRETKEYBYTES - SABER_KEYBYTES, SABER_KEYBYTES);// Remaining part of sk contains a pseudo-random number.
	speedtest_print_cpucycles("Insert(rand(32B))->sk");// This is output when check in crypto_kem_dec() fails.
	speedtest_print_title("crypto_kem_keypair_opt", 0);
	return (0);
}


int crypto_kem_enc_opt_speedtest(unsigned char *c, unsigned char *k,
		const unsigned char *pk) {
	speedtest_print_title("crypto_kem_enc_opt", 1);
	unsigned char kr[64];   // Will contain key, coins (k || random-noise)
	unsigned char buf[64];
	randombytes(buf, 32);
	sha3_256(buf, buf, 32); // BUF[0:31] <-- random message (will be used as the key for client) Note: hash doesnot release system RNG output
	sha3_256(buf + 32, pk, SABER_INDCPA_PUBLICKEYBYTES); // BUF[32:63] <-- Hash(public key);  Multitarget countermeasure for coins + contributory KEM
	sha3_512(kr, buf, 64);   // kr[0:63] <-- Hash(buf[0:63]);// K^ <-- kr[0:31]// noiseseed (r) <-- kr[32:63];
	//indcpa_kem_enc(buf, kr + 32, pk, c);// buf[0:31] contains message; kr[32:63] contains randomness r;
	//void indcpa_kem_enc(unsigned char *message_received, unsigned char *noiseseed,
	//	const unsigned char *pk, unsigned char *ciphertext)
	//========== INDCPA.ENC ===========
	unsigned char *message_received = buf;
	unsigned char *noiseseed = kr + 32;
	//unsigned char *ciphertext = c;
	uint32_t i, j;
	polyvec a[SABER_K];	// skpv;
	//unsigned char seed[SABER_SEEDBYTES];
	polyvec pkcl;	//uint16_t pkcl[SABER_K][SABER_N];		//public key of received by the client
	polyvec skpv1;	//uint16_t skpv1[SABER_K][SABER_N];
	uint16_t message[SABER_KEYBYTES * 8];
	polyvec res;	//uint16_t res[SABER_K][SABER_N];
	uint16_t mod_p = SABER_P - 1;
	uint16_t vprime[SABER_N];

	//unsigned char rec_c[SABER_RECONBYTES_KEM];

	//	for (i = 0; i < SABER_SEEDBYTES; i++) { // extract the seedbytes from Public Key.
	//		seed[i] = pk[ SABER_POLYVECCOMPRESSEDBYTES + i];
	//	}
	//memcpy(seed, pk + SABER_POLYVECCOMPRESSEDBYTES, SABER_SEEDBYTES);
	//GenMatrix_opt(a, seed, 1);
	speedtest_reset_startcpucycles();
	GenMatrix_opt(a, pk + SABER_POLYVECCOMPRESSEDBYTES, 1);
	speedtest_print_cpucycles("Gen(A)");

	GenSecret_opt(&skpv1, noiseseed);	//generate secret from constant-time binomial distribution
	speedtest_print_cpucycles("Gen(s)");
	MatrixVectorMul_rounding_parallel(a, &skpv1, &res, SABER_Q_MASK, 1);
	speedtest_print_cpucycles("MulRouding(A,sk1)->res");
	//POLVECp2BS_opt(ciphertext, &res);
	POLVECp2BS_opt(c, &res);
	speedtest_print_cpucycles("Polyp2Bs(res)->c");

	//------now calculate the v'
	BS2POLVECp_opt(pk, &pkcl);	//unpack the public_key, pkcl is the b in the protocol
	speedtest_print_cpucycles("Bsp2Poly(pk)->pk");
	//skpv1[i][j] = skpv1[i][j] & (mod_p);//可以省略
	VectorMul_opt(&pkcl, &skpv1, mod_p, vprime);
	speedtest_print_cpucycles("VectorMul(pk,sk1)->vprime");
	for (j = 0; j < SABER_KEYBYTES; j++) {	// unpack message_received;
		for (i = 0; i < 8; i++) {
			message[8 * j + i] = ((message_received[j] >> i) & 0x01);
		}
	}
	speedtest_print_cpucycles("unpack(message_received)->message");
	for (i = 0; i < SABER_N; i++) {
		message[i] = (message[i] << 9);	// message encoding
		vprime[i] = vprime[i] + message[i];
	}
	speedtest_print_cpucycles("vprime += (message << 9)");

	ReconDataGen(vprime, c + SABER_POLYVECCOMPRESSEDBYTES);
	speedtest_print_cpucycles("ReconDataGen(vprime)->rec_c");
	//============== INDCPA.ENC END ================

	sha3_256(kr + 32, c, SABER_BYTES_CCA_DEC);	//1088 B
	speedtest_print_cpucycles("Insert(sha256(c(1088B))->kr[32:63]");
	sha3_256(k, kr, 64);
	speedtest_print_cpucycles("Insert(sha256(kr(64B)))->k");
	speedtest_print_title("crypto_kem_enc_opt", 0);
	return (0);
}
//
int crypto_kem_dec_opt_speedtest(unsigned char *k, const unsigned char *c,
		const unsigned char *sk) {
	speedtest_print_title("crypto_kem_dec_opt", 1);
	int i, j, fail;
	unsigned char cmp[SABER_BYTES_CCA_DEC];
	unsigned char buf[64];	//buf[0:31]->message; buf[32:43]->hash(pk)
	unsigned char kr[64];// Will contain key, coins (kr=sha512(message||hash(pk)))
	//kr[0:31] is Kprime(first 32B of hash(mssage||pk) or last 32B random in sk)
	//kr[32:63] is noise(last 32B of hash(mssage||pk) for enc, hash(c) for create K)
	const unsigned char *pk = sk + SABER_INDCPA_SECRETKEYBYTES;

	//==== INDCPA.dec begin =====
	//indcpa_kem_dec(sk, c, buf);	// dec(sk, c) -> message_dec(in buf[0:31])
	//indcpa_kem_dec(*sk, *ciphertext, *message_dec)
	const unsigned char *ciphertext = c;
	polyvec sksv;//uint16_t sksv[SABER_K][SABER_N]; //secret key of the server
	polyvec pksv;//uint16_t pksv[SABER_K][SABER_N];
	//uint8_t recon_ar[SABER_RECONBYTES_KEM];
	uint16_t message_dec_unpacked[SABER_KEYBYTES * 8];// one element containes on decrypted bit;
	uint16_t mod_p = SABER_P - 1;
	uint16_t v[SABER_N];
	speedtest_reset_startcpucycles();

	BS2POLVECq_opt(sk, &sksv);//sksv is the secret-key
	speedtest_print_cpucycles("Bs2Poly(sk)->skpv");
	BS2POLVECp_opt(ciphertext, &pksv);//pksv is the ciphertext
	speedtest_print_cpucycles("Bs2Poly(ciphertext)->pkpv");
	//sksv[i][j] = sksv[i][j] & (mod_p);//省略
	VectorMul_opt(&pksv, &sksv, mod_p, v);//vector-vector scalar multiplication with mod p
	speedtest_print_cpucycles("VectorMul(pkpv,skpv)->v");

	//recon_ar = ciphertext + SABER_POLYVECCOMPRESSEDBYTES
	Recon(v, ciphertext + SABER_POLYVECCOMPRESSEDBYTES, message_dec_unpacked);
	speedtest_print_cpucycles("Recon(v,recon_data)->message_dec_unpacked");
	POL2MSG_opt(message_dec_unpacked, buf);//pack(message_dec)->buf[0:31]
	speedtest_print_cpucycles("Poly2Msg(message_dec_unpacked)->buf[0:31]");
	//=== INDCPA.dec end =====

	// hash(pk) is in sk last[-64:-32]
	memcpy(buf + 32, sk + SABER_SECRETKEYBYTES - 64, 32);//hash(pk) -> buf[32:63];
	speedtest_print_cpucycles("Copy(hash(pk)_in_sk)->buf[32:63]");
	sha3_512(kr, buf, 64);//sha256(buf[0:63])->kr[0:63]; kr=sha512(message||pk)
	speedtest_print_cpucycles("sha512(buff[0:63]:message_dec||hash(pk))->kr");

	//(buf[0:31] is message, kr+32 is noise, pk) -> cmp (cmp is the ciphertext to be compared)

	//========== INDCPA.enc begin ======
	//indcpa_kem_enc(buf, kr + 32, pk, cmp);
	//indcpa_kem_enc(*message_received, *noiseseed, *pk, *ciphertext)
	unsigned char *message_received = buf;//buf[0:31]
	unsigned char *noiseseed = kr + 32;
	//unsigned char *ciphertext = c;
	polyvec a[SABER_K];// skpv;
	//unsigned char seed[SABER_SEEDBYTES];
	polyvec pkcl;//uint16_t pkcl[SABER_K][SABER_N];		//public key of received by the client
	polyvec skpv1;//uint16_t skpv1[SABER_K][SABER_N];
	uint16_t message[SABER_KEYBYTES * 8];
	polyvec res;//uint16_t res[SABER_K][SABER_N];
	//uint16_t mod_p = SABER_P - 1;
	uint16_t vprime[SABER_N];

	//unsigned char rec_c[SABER_RECONBYTES_KEM];

	//seed = pk + SABER_POLYVECCOMPRESSEDBYTES
	speedtest_reset_startcpucycles();
	GenMatrix_opt(a, pk + SABER_POLYVECCOMPRESSEDBYTES, 1);
	speedtest_print_cpucycles("Gen_A(seed_in_pk)->A");

	GenSecret_opt(&skpv1, noiseseed);//generate secret from constant-time binomial distribution
	speedtest_print_cpucycles("Gen_s(noise=kr[32:63])->skpv1");

	MatrixVectorMul_rounding_parallel(a, &skpv1, &res, SABER_Q_MASK, 1);
	speedtest_print_cpucycles("MatrixMulRounding(A,skpv1)->res");

	POLVECp2BS_opt(cmp, &res);
	speedtest_print_cpucycles("Poly2Bs(res)->cmp");

	//------now calculate the v'
	BS2POLVECp_opt(pk, &pkcl);//unpack the public_key, pkcl is the b in the protocol
	speedtest_print_cpucycles("Bs2Poly(pk)->pkpv1(pkcl)");
	//skpv1[i][j] = skpv1[i][j] & (mod_p);//可以省略
	VectorMul_opt(&pkcl, &skpv1, mod_p, vprime);
	speedtest_print_cpucycles("VectorMul(skpv1,pkpv1)->vprime");

	for (j = 0; j < SABER_KEYBYTES; j++) {	// unpack message_received;
		for (i = 0; i < 8; i++) {
			message[8 * j + i] = ((message_received[j] >> i) & 0x01);
		}
	}
	speedtest_print_cpucycles("unpack(message_received in buf[0:31])->message");
	for (i = 0; i < SABER_N; i++) {
		message[i] = (message[i] << 9);	// message encoding
		vprime[i] = vprime[i] + message[i];
	}
	speedtest_print_cpucycles("vprime+=(message<<9)");

	ReconDataGen(vprime, cmp + SABER_POLYVECCOMPRESSEDBYTES);
	speedtest_print_cpucycles("ReconDataGen(vprime)->cmp[3*320,~]");
	//============== INDCPA.ENC END ================

	//if c==cmp then k = hash(kr_r,c) else k = hash(z,c)  (z is in sk last 32B)
	fail = verify(c, cmp, SABER_BYTES_CCA_DEC);//0 for equal strings, 1 for non-equal strings
	//fail = memcmp(c, cmp, SABER_BYTES_CCA_DEC);//a little slower than verify
	speedtest_print_cpucycles("verify(c==cmp)");

	sha3_256(kr + 32, c, SABER_BYTES_CCA_DEC);// overwrite coins in kr with h(c)
	speedtest_print_cpucycles("sha256(c)->kr[32:63]");

	//cmov(kr, sk + SABER_SECRETKEYBYTES - SABER_KEYBYTES, SABER_KEYBYTES, fail);
	if (fail == 0) {   //success
	} else {
		//if not success, random in sk last 32B -> kr[0:31]
		memcpy(kr, sk + SABER_SECRETKEYBYTES - SABER_KEYBYTES, SABER_KEYBYTES);
	}
	speedtest_print_cpucycles("check verify");

	sha3_256(k, kr, 64);    // hash concatenation of pre-k and h(c) to k
	speedtest_print_cpucycles("sha256(kr[0:63]->K(32B)");
	speedtest_print_title("crypto_kem_dec_opt", 0);
	return (0);
}

