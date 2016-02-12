// @author Alexander Pletzer
//
// This code is provided with the hope that 
// it will be useful. It comes with no guarantee whatsoever.

#ifndef LAPACK_NAMES_H
#define LAPACK_NAMES_H

// Fortran
#define FC_GLOBAL(name,NAME) name##_

#ifdef FC_GLOBAL
#define _GESV_ FC_GLOBAL(dgesv,DGESV)
#define _GELS_ FC_GLOBAL(dgels,DGELS)
#define _GETRF_ FC_GLOBAL(dgetrf,DGETRF)
#define _GELSS_ FC_GLOBAL(dgelss,DGELSS)
#define _GELSY_ FC_GLOBAL(dgelsy,DGELSY)
#else
// default mangling, used when CMake failed to determine the mangling
#define _GESV_ dgesv_
#define _GELS_ dgels_
#define _GETRF_ dgetrf_
#define _GELSS_ dgelss_
#define _GELSY_ dgelsy_
#endif

//
// Lapack subroutines
//

#if defined(MKL_FOUND) && defined(LINK_WITH_MKL)
#include <mkl.h>
#else

// solution to a dense linear system
#ifdef __cplusplus
extern "C" {
#endif
void _GESV_(int *, int *,
            double *, int *, int *,
            double *, int *, int *);

// least square solves
void _GELS_(const char*,
            const int*, const int*, const int*,
            double*, const int*,
            double*, const int*,
            double*, const int*,
            int *);

// LU factorization
void _GETRF_(int* m, int* n, double* a,
             int* lda, int* ipv, int* info);

// minimum norm solution to a linear least-square problem
// using a singular value decomposition
void _GELSS_(int*, int*, int*, double*, int*,
             double*, int*, double*, double*,
             int*, double*, int*, int*);

// computes the minimum-norm solution to a real linear least
// squares problem with possible deficient rank matrix
void _GELSY_(int*, int*, int*, double*, int*, double*, int*, int*,
         double*, int*, double*, int*, int*);
#ifdef __cplusplus
}
#endif
#endif // mkl

#endif // LAPACk_NAMES_H
