/*
 * bignum_lite.h
 *  精简以下 大整数加法和减法
 *  Created on: 2019年1月3日
 *      Author: Mai
 */

#ifndef BIGNUM_LITE_H_
#define BIGNUM_LITE_H_

//#include "bignum.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef   signed long t_sint;
typedef unsigned long t_uint;

#define MPI_CHK(f) if( ( ret = f ) != 0 ) goto cleanup

/*
 * Maximum size MPIs are allowed to grow to in number of limbs.
 */
#define POLARSSL_MPI_MAX_LIMBS                             10000

/*
 * Maximum window size used for modular exponentiation. Default: 6
 * Minimum value: 1. Maximum value: 6.
 *
 * Result is an array of ( 2 << POLARSSL_MPI_WINDOW_SIZE ) MPIs used
 * for the sliding window calculation. (So 64 by default)
 *
 * Reduction in size, reduces speed.
 */
#define POLARSSL_MPI_WINDOW_SIZE                           6        /**< Maximum windows size used. */

/*
 * Maximum size of MPIs allowed in bits and bytes for user-MPIs.
 * ( Default: 512 bytes => 4096 bits )
 *
 * Note: Calculations can results temporarily in larger MPIs. So the number
 * of limbs required (POLARSSL_MPI_MAX_LIMBS) is higher.
 */
#define POLARSSL_MPI_MAX_SIZE                              512//512      /**< Maximum number of bytes for usable MPIs. */
#define POLARSSL_MPI_MAX_BITS                              ( 8 * POLARSSL_MPI_MAX_SIZE )    /**< Maximum number of bits for usable MPIs. */


/**
 * \brief          MPI structure
 */
typedef struct
{
    int s;              /*!<  integer sign      */
    size_t n;           /*!<  total # of limbs  */
    t_uint *p;          /*!<  pointer to limbs  */
}
mpi;




/**
 * \brief           Initialize one MPI
 *
 * \param X         One MPI to initialize.
 */
void mpi_init( mpi *X );

/**
 * \brief          Unallocate one MPI
 *
 * \param X        One MPI to unallocate.
 */
void mpi_free( mpi *X );

/**
 * \brief          Enlarge to the specified number of limbs
 *
 * \param X        MPI to grow
 * \param nblimbs  The target number of limbs
 *
 * \return         0 if successful,
 *                 POLARSSL_ERR_MPI_MALLOC_FAILED if memory allocation failed
 */
int mpi_grow( mpi *X, size_t nblimbs );

/**
 * \brief          Baseline multiplication: X = A * B
 *
 * \param X        Destination MPI
 * \param A        Left-hand MPI
 * \param B        Right-hand MPI
 *
 * \return         0 if successful,
 *                 POLARSSL_ERR_MPI_MALLOC_FAILED if memory allocation failed
 */
int mpi_mul_mpi( mpi *X, const mpi *A, const mpi *B );

void mpi_add_abs_lite(mpi *X, const mpi *A, const mpi *B);
void mpi_sub_abs_lite(mpi *X, const mpi *A, const mpi *B);
void mpi_shift_right_1bit_lite(mpi *X) ;


#ifdef __cplusplus
}
#endif

#endif /* BIGNUM_LITE_H_ */
