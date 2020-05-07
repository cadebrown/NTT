/* ntt.h - Number Theoretic Transform (NTT) library header
 *
 * 
 * 
 * 
 */

#pragma once
#ifndef NTT_H__
#define NTT_H__

/* configured file */
#include <ntt-config.h>

/* std includes */
#include <stdio.h>
#include <stdlib.h>

#include <stdbool.h>
#include <stdint.h>

#include <string.h>

#include <math.h>



/* NTT types */


// ntt_plan_gemm_t - plan for GEMM (Matrix-Multiplication) based NTT
typedef struct {

    // N, the number of points in the transform
    int64_t N;

    // p, the prime number related to the transform, such that:
    //   p = Nk+1
    int64_t p;

    // N^-1 (mod p)
    int64_t N_inv;

    // the matrix for the forward NTT
    // (size NxN)
    int64_t* mNTT;

    // the matrix for the inverse NTT
    // (size NxN)
    int64_t* mINTT;


} ntt_plan_gemm_t;

// generate the empty GEMM-based plan
#define NTT_PLAN_GEMM_EMPTY ((ntt_plan_gemm_t){ .N = 0, .p = 0, .mNTT = NULL, .mINTT = NULL })

// Initialize a GEMM-based plan, with 'N' points, mod 'p'
// NOTE: p = Nk + 1, or, if p==0, then 'p' will be calculated as the smallest
//   prime of the form (Nk+1)
void ntt_plan_gemm_init(ntt_plan_gemm_t* plan, int64_t N, int64_t p);

// Do forward NTT:
// out = NTT(inp)
void ntt_plan_gemm_NTT(ntt_plan_gemm_t* plan, int64_t* inp, int64_t* out);

// Do inverse NTT (INTT):
// out = INTT(inp)
void ntt_plan_gemm_INTT(ntt_plan_gemm_t* plan, int64_t* inp, int64_t* out);


// ntt_plan_bfly_t - plan for butterfly-based NTT codes
typedef struct {

    // N, the number of points in the transform
    int64_t N;

    // p, the prime number related to the transform, such that:
    //   p = Nk+1
    int64_t p;

    // N^-1 (mod p)
    int64_t N_inv;

    // twiddle factors:
    // W = 1, w, w^2, w^3, ... w^(n-1)
    // W = 1, w^-1, w^-2, w^-3, ... w^(1-n)
    // Where 'w' is an 'N'th root of unity
    int64_t* W;
    int64_t* IW;

} ntt_plan_bfly_t;

#define NTT_PLAN_BFLY_EMPTY ((ntt_plan_bfly_t){ .N = 0, .p = 0, .W = 0 })

//
void ntt_plan_bfly_init(ntt_plan_bfly_t* plan, int64_t N, int64_t p);

// Do forward NTT:
// out = NTT(inp)
void ntt_plan_bfly_NTT(ntt_plan_bfly_t* plan, int64_t* inp, int64_t* out);

// Do inverse NTT (INTT):
// out = INTT(inp)
void ntt_plan_bfly_INTT(ntt_plan_bfly_t* plan, int64_t* inp, int64_t* out);



/* NTT NT utils */

// Compute gcd(a, b), the largest number which divides into both 'a' and 'b'
NTT_API int64_t ntt_gcd(int64_t a, int64_t b);

// Compute gcd(a, b), AND calculate numbers 'x' and 'y' such that:
//   x * a + y * b = gcd(a, b)
// NOTE: xy[0] is x, xy[1] is y
NTT_API int64_t ntt_egcd(int64_t a, int64_t b, int64_t* xy);

// Compute a^-1 (mod N), which, by definition means: (a * (a^-1)) = 1 (mod N)
// NOTE: Returns '0' if 'a' is not invertible (mod N)
NTT_API int64_t ntt_modinv(int64_t a, int64_t N);

// Compute a^b (mod N) (always returning positive)
// NOTE: If 'b<0', then modular inversing is used. The result will be '-1' if no
//   modular inverse exists
NTT_API int64_t ntt_modpow(int64_t a, int64_t b, int64_t m);

// Return whether or not 'a' is a prime number
// NOTE: Uses a deterministic variant of the miller-rabin test, and so should be fairly fast
NTT_API bool ntt_isprime(int64_t a);

// Compute Euler's Totient function (phi(n))
NTT_API int64_t ntt_tot(int64_t n);


// Compute the first root of unity, w, such that:
// w^i != 1 for i = 1 through phi(n), but that w^phi(n) = 1 (mod n)
// Or, returns '0' if there was no root of unity
NTT_API int64_t ntt_prim_root_unity(int64_t n);


// Compute an 'nth' root ofunity modulo p
NTT_API int64_t ntt_nth_root_unity(int64_t n, int64_t p);

/* Advanced utils like factoring */

// Factor a number 'n', but return unique unsorted prime factors (UUP)
// Return the number of factors, and set (*facts)[:] to the factors
// NOTE: *facts can be NULL, and will be allocated using 'realloc',
//   so only use 'free' with it afterwards
NTT_API int ntt_factor_uup(int64_t n, int64_t** facts);



/* Actual NTT algorithm */


#endif /* NTT_H__ */
