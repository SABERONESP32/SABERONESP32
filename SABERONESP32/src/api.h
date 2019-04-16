#ifndef API_H
#define API_H

// Available algorithms for different security levels
#define LightSaber 1
#define Saber 2
#define FireSaber 3

#define CRYPTO_ALGNAME "Saber"

// Change the algorithm name 
//#define SABER_TYPE LightSaber
#define SABER_TYPE Saber
//#define SABER_TYPE FireSaber

//  Set these three values apropriately for your algorithm
#if SABER_TYPE == LightSaber
	#define CRYPTO_SECRETKEYBYTES 1568
	#define CRYPTO_PUBLICKEYBYTES (2*320+32)
	#define CRYPTO_BYTES 32
	#define CRYPTO_CIPHERTEXTBYTES 736
	#define Saber_type 1
#elif SABER_TYPE == Saber
	#define CRYPTO_SECRETKEYBYTES 2304
	#define CRYPTO_PUBLICKEYBYTES (3*320+32)
	#define CRYPTO_BYTES 32
	#define CRYPTO_CIPHERTEXTBYTES 1088
	#define Saber_type 2
#elif SABER_TYPE == FireSaber
	#define CRYPTO_SECRETKEYBYTES 3040
	#define CRYPTO_PUBLICKEYBYTES (4*320+32)
	#define CRYPTO_BYTES 32
	#define CRYPTO_CIPHERTEXTBYTES 1472
	#define Saber_type 3
#endif


#ifdef __cplusplus
extern "C" {
#endif

//int crypto_kem_keypair(unsigned char *pk, unsigned char *sk);
//int crypto_kem_enc(unsigned char *ct, unsigned char *ss, const unsigned char *pk);
//int crypto_kem_dec(unsigned char *ss, const unsigned char *ct, const unsigned char *sk);

#ifdef __cplusplus
}
#endif

#endif /* api_h */