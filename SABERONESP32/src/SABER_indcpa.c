#include <string.h>
#include <stdint.h>
#include "SABER_indcpa.h"
#include "poly.h"
#include "pack_unpack.h"
#include "poly_mul.h"
#include "recon.h"
#include "rng.h"
#include "fips202.h"
#include "kem_opt.h"
#include "SABER_params.h"
#include "polymul_hw.h"

#define TYPE_SIGLECORE  0//   KEM_UNROLLED_TYPE_SIGLECORE
#define TYPE_SPEED      1//KEM_UNROLLED_TYPE_SPEED
#define TYPE_DUALCORE   2//KEM_UNROLLED_TYPE_DUALCORE

#define TYPE_DEF  0//KEM_UNROLLED_TYPE

/*-----------------------------------------------------------------------------------
 This routine generates a=[Matrix K x K] of 256-coefficient polynomials
 -------------------------------------------------------------------------------------*/

void MatrixVectorMul(polyvec *a, uint16_t skpv[SABER_K][SABER_N], uint16_t res[SABER_K][SABER_N], uint16_t mod,
		int16_t transpose) {
	MatrixVectorMul_parallel(a, skpv, res, mod, transpose);
//	uint16_t i;
//	if (transpose == 1) {
//		for (i = 0; i < SABER_K; i++) {
//			saber_karatsuba_parallel_vector((uint16_t *) &a[0].vec[i], (uint16_t *) &a[1].vec[i],
//					(uint16_t *) &a[2].vec[i], skpv[0], skpv[1], skpv[2], res[i]);
//		}
//	} else {
//		//transpose=0
//		for (i = 0; i < SABER_K; i++) {
//			saber_karatsuba_parallel_vector((uint16_t *) &a[i].vec[0], (uint16_t *) &a[i].vec[1],
//					(uint16_t *) &a[i].vec[2], skpv[0], skpv[1], skpv[2], res[i]);
//		}
//	}

}

//VectorMul(pkcl, skpv1, mod_p, vprime);
//(pkcl: pk received by client; skpv1: sk gen by client; vprime: v' )
void VectorMul(uint16_t pkcl[SABER_K][SABER_N], uint16_t skpv[SABER_K][SABER_N],
		uint16_t mod, uint16_t res[SABER_N]) {
	saber_karatsuba_parallel_vector(pkcl[0], pkcl[1],
			pkcl[2], skpv[0], skpv[1], skpv[2], res);
}

void POL2MSG(uint16_t *message_dec_unpacked, unsigned char *message_dec) {
	int32_t i, j;
	for (j = 0; j < SABER_KEYBYTES; j++) {
		message_dec[j] = 0;
		for (i = 0; i < 8; i++)
			message_dec[j] = message_dec[j] | (message_dec_unpacked[j * 8 + i] << i);
	}
}

void GenMatrix(polyvec *a, const unsigned char *seed) {
	unsigned int one_vector = 13 * SABER_N / 8;
	unsigned int byte_bank_length = SABER_K * SABER_K * one_vector;
	unsigned char buf[byte_bank_length];

	uint16_t temp_ar[SABER_N];

	int i, j, k;
	const uint16_t mod = SABER_Q_MASK; //(SABER_Q-1);

	shake128(buf, byte_bank_length, seed, SABER_SEEDBYTES);

	for (i = 0; i < SABER_K; i++) {
		for (j = 0; j < SABER_K; j++) {
			BS2POLq(buf + (i * SABER_K + j) * one_vector, temp_ar);
			for (k = 0; k < SABER_N; k++) {
				a[i].vec[j].coeffs[k] = (temp_ar[k]) & mod;
			}
		}
	}
}

void GenMatrix_dualcore(polyvec *a, const unsigned char *seed, const int column_first) {
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

#if TYPE_DEF == TYPE_SIGLECORE

void indcpa_kem_keypair(unsigned char *pk, unsigned char *sk) {

	polyvec a[SABER_K];

	uint16_t skpv[SABER_K][SABER_N];

	unsigned char seed[SABER_SEEDBYTES]; //(32B)
	unsigned char noiseseed[SABER_COINBYTES];//(32B)
	//int32_t i, j;

	uint16_t res[SABER_K][SABER_N];

	randombytes(seed, SABER_SEEDBYTES);//(random 32B)
	shake128(seed, SABER_SEEDBYTES, seed, SABER_SEEDBYTES);// for not revealing system RNG state
	randombytes(noiseseed, SABER_COINBYTES);//(random 32B)

	GenMatrix(a, seed);//(shake128) sample matrix A

	GenSecret(skpv, noiseseed);//(shake128) generate secret from constant-time binomial distribution

	//------do the matrix vector multiplication and rounding--------

	MatrixVectorMul_rounding_parallel(a, skpv, res, SABER_Q_MASK, 0);

	//------------------unload and pack sk=3 x (256 coefficients of 14 bits)-------

	POLVECq2BS(sk, skpv);

	//------------------unload and pack pk=256 bits seed and 3 x (256 coefficients of 11 bits)-------

	POLVECp2BS(pk, res);// load the public-key coefficients

//	for (i = 0; i < SABER_SEEDBYTES; i++) { // now load the seedbytes in PK. Easy since seed bytes are kept in byte format.
//		pk[SABER_POLYVECCOMPRESSEDBYTES + i] = seed[i];
//	}
	//(opt)
	memcpy(pk + SABER_POLYVECCOMPRESSEDBYTES, seed, SABER_SEEDBYTES);

}

void indcpa_kem_enc(unsigned char *message_received, unsigned char *noiseseed,
		const unsigned char *pk, unsigned char *ciphertext) {
	uint32_t i, j, k;
	polyvec a[SABER_K];		// skpv;
	unsigned char seed[SABER_SEEDBYTES];
	uint16_t pkcl[SABER_K][SABER_N];//public key of received by the client
	uint16_t skpv1[SABER_K][SABER_N];
	uint16_t message[SABER_KEYBYTES * 8];
	uint16_t res[SABER_K][SABER_N];
	//uint16_t acc[SABER_N];
	//uint16_t mod=SABER_Q-1;
	uint16_t mod_p = SABER_P - 1;

	uint16_t vprime[SABER_N];

	//uint16_t ciphertext_temp[SABER_N];

	unsigned char rec_c[SABER_RECONBYTES_KEM];

	for (i = 0; i < SABER_SEEDBYTES; i++) { // extract the seedbytes from Public Key.
		seed[i] = pk[ SABER_POLYVECCOMPRESSEDBYTES + i];
	}

	GenMatrix(a, seed);

	GenSecret(skpv1, noiseseed); //generate secret from constant-time binomial distribution

	//-----------------matrix-vector multiplication and rounding

//	for (i = 0; i < SABER_K; i++) {
//		for (j = 0; j < SABER_N; j++) {
//			res[i][j] = 0;
//		}
//	}

//	MatrixVectorMul(a, skpv1, res, SABER_Q_MASK, 1);

	//-----now rounding

//	for (i = 0; i < SABER_K; i++) { //shift right 3 bits
//		for (j = 0; j < SABER_N; j++) {
//			res[i][j] = res[i][j] + 4;
//			res[i][j] = (res[i][j] >> 3);
//		}
//	}

	MatrixVectorMul_rounding_parallel(a, skpv1, res, SABER_Q_MASK, 1);

	POLVECp2BS(ciphertext, res);

	//***client matrix-vector multiplication ends********

	//------now calculate the v'

	//-------unpack the public_key

	BS2POLVECp(pk, pkcl);//pkcl is the b in the protocol

//	for (i = 0; i < SABER_N; i++){
//		vprime[i] = 0;
//	}

	for (i = 0; i < SABER_K; i++) {
		for (j = 0; j < SABER_N; j++) {
			skpv1[i][j] = skpv1[i][j] & (mod_p);
		}
	}

	// vector-vector scalar multiplication with mod p
	VectorMul(pkcl, skpv1, mod_p, vprime);

	// unpack message_received;
	for (j = 0; j < SABER_KEYBYTES; j++) {
		for (i = 0; i < 8; i++) {
			message[8 * j + i] = ((message_received[j] >> i) & 0x01);
		}
	}

	// message encoding
	for (i = 0; i < SABER_N; i++) {
		message[i] = (message[i] << 9);
	}

	for (k = 0; k < SABER_N; k++) {
		vprime[k] = vprime[k] + message[k];
	}

	ReconDataGen(vprime, rec_c);

//	for (j = 0; j < SABER_RECONBYTES_KEM; j++) {
//		ciphertext[SABER_POLYVECCOMPRESSEDBYTES + j] = rec_c[j];
//	}
	//(opt)
	memcpy(ciphertext + SABER_POLYVECCOMPRESSEDBYTES, rec_c, SABER_RECONBYTES_KEM);

}

void indcpa_kem_dec(const unsigned char *sk, const unsigned char *ciphertext,
		unsigned char message_dec[]) {

	uint32_t i, j;

	uint16_t sksv[SABER_K][SABER_N]; //secret key of the server

	uint16_t pksv[SABER_K][SABER_N];

	//uint16_t recon_ar[SABER_N];
	uint8_t recon_ar[SABER_RECONBYTES_KEM];
	uint16_t message_dec_unpacked[SABER_KEYBYTES * 8];// one element containes on decrypted bit;

	//uint16_t acc[SABER_N];

	uint16_t mod_p = SABER_P - 1;

	uint16_t v[SABER_N];

	BS2POLVECq(sk, sksv);//sksv is the secret-key
	BS2POLVECp(ciphertext, pksv);//pksv is the ciphertext

	// vector-vector scalar multiplication with mod p
	for (i = 0; i < SABER_N; i++) {
		v[i] = 0;
	}

	for (i = 0; i < SABER_K; i++) {
		for (j = 0; j < SABER_N; j++) {
			sksv[i][j] = sksv[i][j] & (mod_p);
		}
	}

	VectorMul(pksv, sksv, mod_p, v);

	//Extraction
//	for (i = 0; i < SABER_RECONBYTES_KEM; i++) {
//		recon_ar[i] = ciphertext[SABER_POLYVECCOMPRESSEDBYTES + i];
//	}
	//(opt)
	memcpy(recon_ar, ciphertext + SABER_POLYVECCOMPRESSEDBYTES, SABER_RECONBYTES_KEM);

	Recon(v, recon_ar, message_dec_unpacked);

	// pack decrypted message
	POL2MSG(message_dec_unpacked, message_dec);
}

#elif TYPE_DEF == TYPE_SPEED

void indcpa_kem_keypair(unsigned char *pk, unsigned char *sk) {
	speedtest_print_title("IND_CPA.Gen", 1);
	polyvec a[SABER_K];
	uint16_t skpv[SABER_K][SABER_N];
	unsigned char seed[SABER_SEEDBYTES]; //(32B)
	unsigned char noiseseed[SABER_COINBYTES];//(32B)
	uint16_t res[SABER_K][SABER_N];
	speedtest_reset_startcpucycles();

	randombytes(seed, SABER_SEEDBYTES);//(random 32B)
	speedtest_print_cpucycles("randombytes");

	shake128(seed, SABER_SEEDBYTES, seed, SABER_SEEDBYTES);// for not revealing system RNG state
	speedtest_print_cpucycles("shake128");

	randombytes(noiseseed, SABER_COINBYTES);//(random 32B)
	speedtest_print_cpucycles("randombytes");

	GenMatrix(a, seed);//(shake128) sample matrix A
	speedtest_print_cpucycles("GenMatrix(seed)->A");

	GenSecret(skpv, noiseseed);//(shake128) generate secret from constant-time binomial distribution
	speedtest_print_cpucycles("GenSecret(noise)->sk");

	MatrixVectorMul_rounding_parallel(a, skpv, res, SABER_Q_MASK, 0);//matrix vector multiplication and rounding
	speedtest_print_cpucycles("MatrixVectorMul_rounding(A,sk)->res");

	POLVECq2BS(sk, skpv);//unload and pack sk=3 x (256 coefficients of 14 bits)
	speedtest_print_cpucycles("POLVECq2BS(sk)->sk");

	//unload and pack pk=256 bits seed and 3 x (256 coefficients of 11 bits)
	POLVECp2BS(pk, res);// load the public-key coefficients
	speedtest_print_cpucycles("POLVECp2BS(res)->pk");
	//now load the seedbytes in PK. Easy since seed bytes are kept in byte format.
	memcpy(pk + SABER_POLYVECCOMPRESSEDBYTES, seed, SABER_SEEDBYTES);
	speedtest_print_cpucycles("pk.append(seed)");

	speedtest_print_title("IND_CPA.Gen", 0);

}

void indcpa_kem_enc(unsigned char *message_received, unsigned char *noiseseed,
		const unsigned char *pk, unsigned char *ciphertext) {
	speedtest_print_title("IND_CPA.Enc", 1);
	speedtest_reset_startcpucycles();
	speedtest_print_cpucycles("input(pk, noise)");
	uint32_t i, j, k;
	polyvec a[SABER_K];		// skpv;
	//unsigned char seed[SABER_SEEDBYTES];
	uint16_t pkcl[SABER_K][SABER_N];//public key of received by the client

	uint16_t skpv1[SABER_K][SABER_N];

	uint16_t message[SABER_KEYBYTES * 8];

	uint16_t res[SABER_K][SABER_N];
	//uint16_t acc[SABER_N];  
	//uint16_t mod=SABER_Q-1;
	uint16_t mod_p = SABER_P - 1;

	uint16_t vprime[SABER_N];

	//uint16_t ciphertext_temp[SABER_N];

	unsigned char rec_c[SABER_RECONBYTES_KEM];

	speedtest_reset_startcpucycles();

//	for (i = 0; i < SABER_SEEDBYTES; i++) { // extract the seedbytes from Public Key.
//		seed[i] = pk[ SABER_POLYVECCOMPRESSEDBYTES + i];
//	}//省略
	//memcpy(seed, pk + SABER_POLYVECCOMPRESSEDBYTES, SABER_SEEDBYTES);
	unsigned char *seed = &pk[SABER_POLYVECCOMPRESSEDBYTES];
	speedtest_print_cpucycles("extract seed from pk");

	GenMatrix(a, seed);
	speedtest_print_cpucycles("GenMatrix(seed)->A");

	GenSecret(skpv1, noiseseed);//generate secret from constant-time binomial distribution
	speedtest_print_cpucycles("GenSecret(noise)->sk1");

	//-----------------matrix-vector multiplication and rounding

	MatrixVectorMul_rounding_parallel(a, skpv1, res, SABER_Q_MASK, 1);
	speedtest_print_cpucycles("MatrixVectorMul_rounding(A,sk1)->res");

	POLVECp2BS(ciphertext, res);
	speedtest_print_cpucycles("POLVECp2BS(ciphertext, res)");

	//***client matrix-vector multiplication ends********
	//------now calculate the v'
	//-------unpack the public_key

	BS2POLVECp(pk, pkcl);//pkcl is the b in the protocol
	speedtest_print_cpucycles("BS2POLVECp(pk, pkcl)");

//	for (i = 0; i < SABER_K; i++) {
//		for (j = 0; j < SABER_N; j++) {
//			skpv1[i][j] = skpv1[i][j] & (mod_p);
//		}
//	}//可以省略
//	speedtest_print_cpucycles("skpv1 & mod_p");

	// vector-vector scalar multiplication with mod p
	VectorMul(pkcl, skpv1, mod_p, vprime);
	speedtest_print_cpucycles("VectorMul(pkcl, skpv1)");

	// unpack message_received;
	for (j = 0; j < SABER_KEYBYTES; j++) {
		for (i = 0; i < 8; i++) {
			message[8 * j + i] = ((message_received[j] >> i) & 0x01);
		}
	}
	speedtest_print_cpucycles("unpack message_received");

	// message encoding
	for (i = 0; i < SABER_N; i++) {
		message[i] = (message[i] << 9);
	}
	speedtest_print_cpucycles("message encoding");

	for (k = 0; k < SABER_N; k++) {
		vprime[k] = vprime[k] + message[k];
	}
	speedtest_print_cpucycles("vprime+=message");

	ReconDataGen(vprime, rec_c);
	speedtest_print_cpucycles("ReconDataGen(vprime, rec_c)");

	memcpy(ciphertext + SABER_POLYVECCOMPRESSEDBYTES, rec_c, SABER_RECONBYTES_KEM);
	speedtest_print_cpucycles("ciphertext.append(rec_c)");

	speedtest_print_title("IND_CPA.Enc", 0);
}

void indcpa_kem_dec(const unsigned char *sk, const unsigned char *ciphertext,
		unsigned char message_dec[]) {
	//uint32_t i, j;
	uint16_t sksv[SABER_K][SABER_N];//secret key of the server
	uint16_t pksv[SABER_K][SABER_N];
	//uint8_t recon_ar[SABER_RECONBYTES_KEM];
	uint16_t message_dec_unpacked[SABER_KEYBYTES * 8];// one element containes on decrypted bit;
	uint16_t mod_p = SABER_P - 1;
	uint16_t v[SABER_N];
	speedtest_print_title("IND_CPA.Dec", 1);
	speedtest_reset_startcpucycles();

	BS2POLVECq(sk, sksv);//sksv is the secret-key
	speedtest_print_cpucycles("BS2POLVECq(sk, sksv)");
	BS2POLVECp(ciphertext, pksv);//pksv is the ciphertext
	speedtest_print_cpucycles("BS2POLVECp(ciphertext, pksv)");

	// vector-vector scalar multiplication with mod p
//	for (i = 0; i < SABER_N; i++) {
//		v[i] = 0;
//	}
//	for (i = 0; i < SABER_K; i++) {
//		for (j = 0; j < SABER_N; j++) {
//			sksv[i][j] = sksv[i][j] & (mod_p);
//		}
//	}//可以省略这一步
	speedtest_print_cpucycles("sksv & mod_p");

	VectorMul(pksv, sksv, mod_p, v);
	speedtest_print_cpucycles("VectorMul(pk,sk)->v");

	//Extraction
//	for (i = 0; i < SABER_RECONBYTES_KEM; i++) {
//		recon_ar[i] = ciphertext[SABER_POLYVECCOMPRESSEDBYTES + i];
//	}
	//(opt)
	//memcpy(recon_ar, ciphertext + SABER_POLYVECCOMPRESSEDBYTES, SABER_RECONBYTES_KEM);
	uint8_t *recon_ar = &ciphertext[SABER_POLYVECCOMPRESSEDBYTES];
	speedtest_print_cpucycles("Extract recon_ar from ciphertext");

	Recon(v, recon_ar, message_dec_unpacked);
	speedtest_print_cpucycles("Recon((v, recon_ar) -> message_dec_unpacked)");

	// pack decrypted message
	POL2MSG(message_dec_unpacked, message_dec);
	speedtest_print_cpucycles("POL2MSG(message_dec_unpacked->message_dec)");
	speedtest_print_title("IND_CPA.Dec", 0);
}

#elif ( TYPE_DEF == TYPE_DUALCORE)

void indcpa_kem_keypair(unsigned char *pk, unsigned char *sk) {
	speedtest_print_title("IND_CPA.Gen", 1);
	polyvec a[SABER_K];
	uint16_t skpv[SABER_K][SABER_N];
	unsigned char seed[SABER_SEEDBYTES]; //(32B)
	unsigned char noiseseed[SABER_COINBYTES]; //(32B)
	uint16_t res[SABER_K][SABER_N];
	speedtest_reset_startcpucycles();

	randombytes(seed, SABER_SEEDBYTES); //(random 32B)
	speedtest_print_cpucycles("randombytes");

	shake128(seed, SABER_SEEDBYTES, seed, SABER_SEEDBYTES); // for not revealing system RNG state
	speedtest_print_cpucycles("shake128");

	randombytes(noiseseed, SABER_COINBYTES); //(random 32B)
	speedtest_print_cpucycles("randombytes");

	GenMatrix_dualcore(a, seed, 0);	//(shake128) sample matrix A
	speedtest_print_cpucycles("GenMatrix(seed)->A");

	GenSecret(skpv, noiseseed);	//(shake128) generate secret from constant-time binomial distribution
	speedtest_print_cpucycles("GenSecret(noise)->sk");

	MatrixVectorMul_rounding_parallel(a, skpv, res, SABER_Q_MASK, 0);	//matrix vector multiplication and rounding
	speedtest_print_cpucycles("MatrixVectorMul_rounding(A,sk)->res");

	POLVECq2BS(sk, skpv);	//unload and pack sk=3 x (256 coefficients of 14 bits)
	speedtest_print_cpucycles("POLVECq2BS(sk)->sk");

	//unload and pack pk=256 bits seed and 3 x (256 coefficients of 11 bits)
	POLVECp2BS(pk, res); // load the public-key coefficients
	speedtest_print_cpucycles("POLVECp2BS(res)->pk");
	//now load the seedbytes in PK. Easy since seed bytes are kept in byte format.
	memcpy(pk + SABER_POLYVECCOMPRESSEDBYTES, seed, SABER_SEEDBYTES);
	speedtest_print_cpucycles("pk.append(seed)");

	speedtest_print_title("IND_CPA.Gen", 0);

}

void indcpa_kem_enc(unsigned char *message_received, unsigned char *noiseseed,
		const unsigned char *pk, unsigned char *ciphertext) {
	speedtest_print_title("IND_CPA.Enc", 1);
	speedtest_reset_startcpucycles();
	speedtest_print_cpucycles("input(pk, noise)");
	uint32_t i, j, k;
	polyvec a[SABER_K];		// skpv;
	//unsigned char seed[SABER_SEEDBYTES];
	uint16_t pkcl[SABER_K][SABER_N]; 	//public key of received by the client

	uint16_t skpv1[SABER_K][SABER_N];

	uint16_t message[SABER_KEYBYTES * 8];

	uint16_t res[SABER_K][SABER_N];
	//uint16_t acc[SABER_N];
	//uint16_t mod=SABER_Q-1;
	uint16_t mod_p = SABER_P - 1;

	uint16_t vprime[SABER_N];

	//uint16_t ciphertext_temp[SABER_N];

	unsigned char rec_c[SABER_RECONBYTES_KEM];

	speedtest_reset_startcpucycles();

//	for (i = 0; i < SABER_SEEDBYTES; i++) { // extract the seedbytes from Public Key.
//		seed[i] = pk[ SABER_POLYVECCOMPRESSEDBYTES + i];
//	}//省略
	//memcpy(seed, pk + SABER_POLYVECCOMPRESSEDBYTES, SABER_SEEDBYTES);
	unsigned char *seed = &pk[SABER_POLYVECCOMPRESSEDBYTES];
	speedtest_print_cpucycles("extract seed from pk");

	GenMatrix_dualcore(a, seed, 1);
	speedtest_print_cpucycles("GenMatrix(seed)->A");

	GenSecret(skpv1, noiseseed); //generate secret from constant-time binomial distribution
	speedtest_print_cpucycles("GenSecret(noise)->sk1");

	//-----------------matrix-vector multiplication and rounding

	MatrixVectorMul_rounding_parallel(a, skpv1, res, SABER_Q_MASK, 1);
	speedtest_print_cpucycles("MatrixVectorMul_rounding(A,sk1)->res");

	POLVECp2BS(ciphertext, res);
	speedtest_print_cpucycles("POLVECp2BS(ciphertext, res)");

	//***client matrix-vector multiplication ends********
	//------now calculate the v'
	//-------unpack the public_key

	BS2POLVECp(pk, pkcl); //pkcl is the b in the protocol
	speedtest_print_cpucycles("BS2POLVECp(pk, pkcl)");

//	for (i = 0; i < SABER_K; i++) {
//		for (j = 0; j < SABER_N; j++) {
//			skpv1[i][j] = skpv1[i][j] & (mod_p);
//		}
//	}//可以省略
//	speedtest_print_cpucycles("skpv1 & mod_p");

	// vector-vector scalar multiplication with mod p
	VectorMul(pkcl, skpv1, mod_p, vprime);
	speedtest_print_cpucycles("VectorMul(pkcl, skpv1)");

	// unpack message_received;
	for (j = 0; j < SABER_KEYBYTES; j++) {
		for (i = 0; i < 8; i++) {
			message[8 * j + i] = ((message_received[j] >> i) & 0x01);
		}
	}
	speedtest_print_cpucycles("unpack message_received");

	// message encoding
	for (i = 0; i < SABER_N; i++) {
		message[i] = (message[i] << 9);
	}
	speedtest_print_cpucycles("message encoding");

	for (k = 0; k < SABER_N; k++) {
		vprime[k] = vprime[k] + message[k];
	}
	speedtest_print_cpucycles("vprime+=message");

	ReconDataGen(vprime, rec_c);
	speedtest_print_cpucycles("ReconDataGen(vprime, rec_c)");

	memcpy(ciphertext + SABER_POLYVECCOMPRESSEDBYTES, rec_c, SABER_RECONBYTES_KEM);
	speedtest_print_cpucycles("ciphertext.append(rec_c)");

	speedtest_print_title("IND_CPA.Enc", 0);
}

void indcpa_kem_dec(const unsigned char *sk, const unsigned char *ciphertext,
		unsigned char message_dec[]) {
	//uint32_t i, j;
	uint16_t sksv[SABER_K][SABER_N]; //secret key of the server
	uint16_t pksv[SABER_K][SABER_N];
	//uint8_t recon_ar[SABER_RECONBYTES_KEM];
	uint16_t message_dec_unpacked[SABER_KEYBYTES * 8];	// one element containes on decrypted bit;
	uint16_t mod_p = SABER_P - 1;
	uint16_t v[SABER_N];
	speedtest_print_title("IND_CPA.Dec", 1);
	speedtest_reset_startcpucycles();

	BS2POLVECq(sk, sksv); //sksv is the secret-key
	speedtest_print_cpucycles("BS2POLVECq(sk, sksv)");
	BS2POLVECp(ciphertext, pksv); //pksv is the ciphertext
	speedtest_print_cpucycles("BS2POLVECp(ciphertext, pksv)");

	// vector-vector scalar multiplication with mod p
//	for (i = 0; i < SABER_N; i++) {
//		v[i] = 0;
//	}
//	for (i = 0; i < SABER_K; i++) {
//		for (j = 0; j < SABER_N; j++) {
//			sksv[i][j] = sksv[i][j] & (mod_p);
//		}
//	}//可以省略这一步
	speedtest_print_cpucycles("sksv & mod_p");

	VectorMul(pksv, sksv, mod_p, v);
	speedtest_print_cpucycles("VectorMul(pk,sk)->v");

	//Extraction
//	for (i = 0; i < SABER_RECONBYTES_KEM; i++) {
//		recon_ar[i] = ciphertext[SABER_POLYVECCOMPRESSEDBYTES + i];
//	}
	//(opt)
	//memcpy(recon_ar, ciphertext + SABER_POLYVECCOMPRESSEDBYTES, SABER_RECONBYTES_KEM);
	uint8_t *recon_ar = &ciphertext[SABER_POLYVECCOMPRESSEDBYTES];
	speedtest_print_cpucycles("Extract recon_ar from ciphertext");

	Recon(v, recon_ar, message_dec_unpacked);
	speedtest_print_cpucycles("Recon((v, recon_ar) -> message_dec_unpacked)");

	// pack decrypted message
	POL2MSG(message_dec_unpacked, message_dec);
	speedtest_print_cpucycles("POL2MSG(message_dec_unpacked->message_dec)");
	speedtest_print_title("IND_CPA.Dec", 0);
}
#endif
