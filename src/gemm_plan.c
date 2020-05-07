/* gemm_plan.c - Matrix-Multiply based NTT plan */

#include "ntt.h"

// create plan with given size
// if p==0, calcaulte it as the smallest prime of the form (Nk+1)
void ntt_plan_gemm_init(ntt_plan_gemm_t* plan, int64_t N, int64_t p) {

    // search for appropriate 'p'
    if (p == 0) {
        p = N + 1;
        while (!ntt_isprime(p)) p += N;
    }

    // p = Nk+1
    int64_t k = (p - 1) / N;

    // store some basic stuff in the plan
    plan->N = N;
    plan->p = p;

    // calculate N^-1 (mod p)
    plan->N_inv = ntt_modinv(N, p);

    // allocate the matrices for the forward and inverse transform
    plan->mNTT = realloc(plan->mNTT, sizeof(*plan->mNTT) * N * N);
    plan->mINTT = realloc(plan->mINTT, sizeof(*plan->mINTT) * N * N);

    // calculate a primitive root of unity
    int64_t rt_p = ntt_prim_root_unity(p);

    // now, calculate N'th root of that (using a trick via 'k')
    int64_t w = ntt_modpow(rt_p, k, p);

    // and the inverse root
    int64_t w_inv = ntt_modinv(w, p);

    // now, calculate according the plan:
    // [1 1 1 ... 1]
    // [1 w w^2 w^3 ...]
    // [1 w^2 w^4 ...]
    // [1 w^3 ...]
    // And, that matrix's inverse
    // Store in row-major order (even though they are symmetric)
    int i;
    for (i = 0; i < N; ++i) {
        int j;
        for (j = 0; j < N; ++j) {
            plan->mNTT[i * N + j] = ntt_modpow(w, i * j, p);
            plan->mINTT[i * N + j] = ntt_modpow(w_inv, i * j, p);
        }
    }
}


// Do forward NTT:
// out = NTT(inp)
void ntt_plan_gemm_NTT(ntt_plan_gemm_t* plan, int64_t* inp, int64_t* out) {
    // store plan variables
    int64_t N = plan->N, p = plan->p;

    // do a matrix multiplication:
    // mNTT * inp -> out
    int64_t i;
    for (i = 0; i < N; ++i) {
        // current row to dot
        int64_t* row = &plan->mNTT[i * N];
        
        // result of the current row dotted with 'inp'
        int64_t r = 0;
        
        int64_t j;
        for (j = 0; j < N; ++j, ++row) {
            r += (*row) * inp[j];
            r %= p;
        }

        // now, set the output
        out[i] = r;
        // and adjust to ensure it is positive
        if (out[i] < 0) out[i] += p;
    }
}

// Do inverse NTT (INTT):
// out = INTT(inp)
void ntt_plan_gemm_INTT(ntt_plan_gemm_t* plan, int64_t* inp, int64_t* out) {

    // store plan variables
    int64_t N = plan->N, p = plan->p, N_inv = plan->N_inv;


    // do a matrix multiplication:
    // mINTT * inp -> out
    int64_t i;
    for (i = 0; i < N; ++i) {
        // current row to dot
        int64_t* row = &plan->mINTT[i * N];
        
        // result of the current row dotted with 'inp'
        int64_t r = 0;
        
        int64_t j;
        for (j = 0; j < N; ++j, ++row) {
            r += (*row) * inp[j];
            r %= p;
        }

        // now, set the output (and multiply by corrective factor)
        out[i] = (r * N_inv) % p;
        // and adjust to ensure it is positive
        if (out[i] < 0) out[i] += p;
    }
}




