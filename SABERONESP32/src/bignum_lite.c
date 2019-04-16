/*
 * bignum_lite.c
 *
 *  Created on: 2019年1月3日
 *  big integer operations based on mbedtls for fixed bit length.
 */

#include "bignum_lite.h"

#define ciL    4 // (sizeof(t_uint))         /* chars in limb  */
#define biL    (ciL << 3)  //即32             /* bits  in limb  */
#define biH    (ciL << 2)               /* half limb size */

/*
 * Initialize one MPI
 */
void mpi_init(mpi *X) {
	if (X == NULL)
		return;

	X->s = 1;
	X->n = 0;
	X->p = NULL;
}

/*
 * Unallocate one MPI
 */
void mpi_free(mpi *X) {
	if (X == NULL)
		return;

	if (X->p != NULL) {
		memset(X->p, 0, X->n * ciL);
		free(X->p);
	}

	X->s = 1;
	X->n = 0;
	X->p = NULL;
}

/*
 * Enlarge to the specified number of limbs（增长到多少个uint32）
 */
int mpi_grow(mpi *X, size_t nblimbs) {
	t_uint *p;

	if (X->n < nblimbs) {
		if ((p = (t_uint *) malloc(nblimbs * ciL)) == NULL)
			return (-1);

		memset(p, 0, nblimbs * ciL);

		if (X->p != NULL) {
			memcpy(p, X->p, X->n * ciL);
			memset(X->p, 0, X->n * ciL);
			free(X->p);
		}

		X->n = nblimbs;
		X->p = p;
	}

	return (0);
}

/*
 * Set value from integer
 */
int mpi_lset(mpi *X, t_sint z) {
	int ret;

	MPI_CHK(mpi_grow(X, 1));
	memset(X->p, 0, X->n * ciL);

	X->p[0] = (z < 0) ? -z : z;
	X->s = (z < 0) ? -1 : 1;

	cleanup:

	return (ret);
}

/*
 * Copy the contents of Y into X
 */
int mpi_copy(mpi *X, const mpi *Y) {
	int ret;
	size_t i;

	if (X == Y)
		return (0);

	for (i = Y->n - 1; i > 0; i--)
		if (Y->p[i] != 0)
			break;
	i++;

	X->s = Y->s;

	MPI_CHK(mpi_grow(X, i));

	memset(X->p, 0, X->n * ciL);
	memcpy(X->p, Y->p, i * ciL);

	cleanup:

	return (ret);
}

//X:=Y,  X,Y长度必须都是n
//void mpi_copy_lite(mpi *X, const mpi *Y) {
//	memcpy(X->p, Y->p, Y->n * ciL);
//}

//X,A,B必须是不同的指针！X,A,B必须是同样的长度！
void mpi_add_abs_lite(mpi *X, const mpi *A, const mpi *B) {
	size_t i, j;
	t_uint *o, *p, c;

	//mpi_copy(X, A);
	//mpi_copy_lite(X, A);
	memcpy(X->p, A->p, A->n * ciL);

	//X->s = 1;

//	for (j = B->n; j > 0; j--)
//		if (B->p[j - 1] != 0)
//			break;

	j = B->n;

	//MPI_CHK(mpi_grow(X, j));

	o = B->p;
	p = X->p;
	c = 0;

	for (i = 0; i < j; i++, o++, p++) {
		*p += c;
		c = (*p < c);
		*p += *o;
		c += (*p < *o);
	}

	while (c != 0) {
		if (i >= X->n) {
			//MPI_CHK(mpi_grow(X, i + 1));
			p = X->p + i;
		}
		*p += c;
		c = (*p < c);
		i++;
	}
}

/*
 * Helper for mpi substraction
 */
static inline void mpi_sub_hlp_lite(size_t n, t_uint *s, t_uint *d) {
	size_t i;
	t_uint c, z;

	for (i = c = 0; i < n; i++, s++, d++) {
		z = (*d < c);
		*d -= c;
		c = (*d < *s) + z;
		*d -= *s;
	}

	while (c != 0) {
		z = (*d < c);
		*d -= c;
		c = z;
		i++;
		d++;
	}
}

/*
 * Unsigned substraction: X = |A| - |B|  (HAC 14.9)
 */
void mpi_sub_abs_lite(mpi *X, const mpi *A, const mpi *B) {
	//mpi TB;
//	int ret;
//	size_t n;

//	mpi_init(&TB);
//	if (X == B) {
//		MPI_CHK(mpi_copy(&TB, B));
//		B = &TB;//这个的本意是传入的X和B相同了，需要保留下来B的数据
//	}
//	if (X != A)
//		MPI_CHK(mpi_copy(X, A));

	//mpi_copy_lite(X, A);
	memcpy(X->p, A->p, A->n * ciL);

	//X->s = 1;// X should always be positive as a result of unsigned substractions.
//	for (n = B->n; n > 0; n--)
//		if (B->p[n - 1] != 0)
//			break;
//	n = B->n;
	mpi_sub_hlp_lite(B->n, B->p, X->p);
//	cleanup:
//	mpi_free(&TB);
}

//右移1比特，即除以2
void mpi_shift_right_1bit_lite(mpi *X) {
	//size_t i, v0, v1;
	t_uint r0 = 0, r1;
	int i;
	for (i = X->n; i > 0; i--) {
		r1 = X->p[i - 1] << 31; // (biL - v1);
		X->p[i - 1] >>= 1; //v1;
		X->p[i - 1] |= r0;
		r0 = r1;
	}
}

#define MULADDC_INIT                    \
{                                       \
    t_uint s0, s1, b0, b1;              \
    t_uint r0, r1, rx, ry;              \
    b0 = ( b << biH ) >> biH;           \
    b1 = ( b >> biH );

#define MULADDC_CORE                    \
    s0 = ( *s << biH ) >> biH;          \
    s1 = ( *s >> biH ); s++;            \
    rx = s0 * b1; r0 = s0 * b0;         \
    ry = s1 * b0; r1 = s1 * b1;         \
    r1 += ( rx >> biH );                \
    r1 += ( ry >> biH );                \
    rx <<= biH; ry <<= biH;             \
    r0 += rx; r1 += (r0 < rx);          \
    r0 += ry; r1 += (r0 < ry);          \
    r0 +=  c; r1 += (r0 <  c);          \
    r0 += *d; r1 += (r0 < *d);          \
    c = r1; *(d++) = r0;

#define MULADDC_STOP                 \
}

static void mpi_mul_hlp(size_t i, t_uint *s, t_uint *d, t_uint b) {
	t_uint c = 0, t = 0;

#if defined(MULADDC_HUIT)
	for(; i >= 8; i -= 8 )
	{
		MULADDC_INIT
		MULADDC_HUIT
		MULADDC_STOP
	}

	for(; i > 0; i-- )
	{
		MULADDC_INIT
		MULADDC_CORE
		MULADDC_STOP
	}
#else
	for (; i >= 16; i -= 16) {
		MULADDC_INIT
			MULADDC_CORE
			MULADDC_CORE
			MULADDC_CORE
			MULADDC_CORE
			MULADDC_CORE
			MULADDC_CORE
			MULADDC_CORE
			MULADDC_CORE

			MULADDC_CORE
			MULADDC_CORE
			MULADDC_CORE
			MULADDC_CORE
			MULADDC_CORE
			MULADDC_CORE
			MULADDC_CORE
			MULADDC_CORE
		MULADDC_STOP
	}

	for (; i >= 8; i -= 8) {
		MULADDC_INIT
			MULADDC_CORE
			MULADDC_CORE
			MULADDC_CORE
			MULADDC_CORE

			MULADDC_CORE
			MULADDC_CORE
			MULADDC_CORE
			MULADDC_CORE
		MULADDC_STOP
	}

	for (; i > 0; i--) {
		MULADDC_INIT
			MULADDC_CORE
		MULADDC_STOP
	}
#endif

	t++;

	do {
		*d += c;
		c = (*d < c);
		d++;
	} while (c != 0);
}
/*
 * Baseline multiplication: X = A * B  (HAC 14.12)
 */
int mpi_mul_mpi(mpi *X, const mpi *A, const mpi *B) {
	int ret;
	size_t i, j;
	mpi TA, TB;

	mpi_init(&TA);
	mpi_init(&TB);

	if (X == A) {
		MPI_CHK(mpi_copy(&TA, A));
		A = &TA;
	}
	if (X == B) {
		MPI_CHK(mpi_copy(&TB, B));
		B = &TB;
	}

	for (i = A->n; i > 0; i--)
		if (A->p[i - 1] != 0)
			break;

	for (j = B->n; j > 0; j--)
		if (B->p[j - 1] != 0)
			break;

	MPI_CHK(mpi_grow(X, i + j));
	MPI_CHK(mpi_lset(X, 0));

	for (i++; j > 0; j--)
		mpi_mul_hlp(i - 1, A->p, X->p + j - 1, B->p[j - 1]);

	X->s = A->s * B->s;		//符号

	cleanup:

	mpi_free(&TB);
	mpi_free(&TA);

	return (ret);
}

