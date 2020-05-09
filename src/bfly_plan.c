/* bfly_plan.c - butterfly-based NTT code, which is faster than GEMM-based */

#include "ntt.h"


// Get the reversed bits of 'x'



// shuffle a sequence with bit-reversal index mapping
static void shuffle_bitrev(int64_t* inp, int64_t N) {
    int64_t oN = N;

    N = 1;
    while (N < oN) {
        N *= 2;
    }

    // not power of 2
    if (N != oN) return;

    int64_t i, j = 0;
    for (i = 1; i < N; ++i) {
        int64_t b = (N >> 1);
        while (j >= b) {
            j -= b;
            b >>= 1;
        }
        j += b;
        if (j > i) {
            int64_t t = inp[i];
            inp[i] = inp[j];
            inp[j] = t;
        }
    }

}

// initialize butterfly-based plan
void ntt_plan_bfly_init(ntt_plan_bfly_t* plan, int64_t N, int64_t p) {

    // search for appropriate 'p'
    if (p == 0) {
        p = N + 1;
        while (!ntt_isprime(p)) p += N;
    }

    // save NTT data
    plan->N = N;
    plan->p = p;
    plan->N_inv = ntt_modinv(N, p);

    // allocate twiddle tables
    plan->W = realloc(plan->W, sizeof(*plan->W) * N);
    plan->IW = realloc(plan->IW, sizeof(*plan->IW) * N);

    // bit reversed
    plan->W_br = realloc(plan->W_br, sizeof(*plan->W_br) * N);

    int64_t k = (p - 1) / N;

    // calculate a primitive root of unity and it's inverse
    int64_t rt_p = ntt_prim_root_unity(p);

    // calculate roots of unity, using NT
    int64_t w = ntt_modpow(rt_p, k, p);
    int64_t w_inv = ntt_modinv(w, p);

    int64_t Wi = 1, Wi_inv = 1;

    // caculate twiddle factors w^i (mod p)
    int64_t i;
    for (i = 0; i < N; ++i) {
        plan->W[i] = Wi;
        plan->IW[i] = Wi_inv;
        Wi = ntt_modmul(Wi, w, p);
        Wi_inv = ntt_modmul(Wi_inv, w_inv, p);
    }

    //memcpy(plan->W_br, plan->W, sizeof(*plan->W) * N);
    //shuffle_bitrev(plan->W_br, N);
}

// Do forward NTT:
// out = NTT(inp)
void ntt_plan_bfly_NTT(ntt_plan_bfly_t* plan, int64_t* inp, int64_t* out) {
    // do in place on output
    memcpy(out, inp, sizeof(*inp) * plan->N);
    shuffle_bitrev(out, plan->N);

    // store plan variables as locals
    int64_t N = plan->N, p = plan->p;

    // temporary variables
    int64_t i, j, U, V;

    // current transform size (powers of 2)
    int64_t m = 1;

    while (m <= N) {
        int64_t m2 = m / 2;
        int64_t bnd = (N / m) * m;

        for (i = 0; i < m2; ++i) {
            // current root of unity
            int64_t wi = plan->W[i * N / m];

            // interio transform
            for (j = i; j < bnd + i; j += m) {
                U = out[j];
                V = (out[j + m2] * wi) % p;

                out[j] = (U + V) % p;
                out[j + m2] = (U - V) % p;
            }
        }

        // keep growing up the transform size
        m *= 2;
    }

    // adjust for modulo 'p'
    for (i = 0; i < plan->N; ++i) {
        if (out[i] < 0) out[i] += p;
    }

}

// Do inverse NTT (INTT):
// out = INTT(inp)
void ntt_plan_bfly_INTT(ntt_plan_bfly_t* plan, int64_t* inp, int64_t* out) {
    // do in place on output
    memcpy(out, inp, sizeof(*inp) * plan->N);
    shuffle_bitrev(out, plan->N);

    // store plan variables as locals
    int64_t N = plan->N, p = plan->p;

    // temp indicies/values
    int64_t i, j, U, V;

    // current transform stride
    int64_t m = 2;

    while (m <= N) {
        // temporary variables for saving
        int64_t m2 = m / 2;
        int64_t bnd = (N / m) * m;

        for (i = 0; i < m2; ++i) {
            // current root of unity
            int64_t wi = plan->IW[i * N / m];

            // calculate interior butterfly
            for (j = i; j < bnd + i; j += m) {
                U = out[j];
                V = (out[j + m2] * wi) % p;

                out[j] = (U + V) % p;
                out[j + m2] = (U - V) % p;
            }
        }

        m *= 2;
    }

    // now, multiply by corrective force and adjust modulo 'p'
    for (i = 0; i < plan->N; ++i) {
        out[i] = (out[i] * plan->N_inv) % p;
        if (out[i] < 0) out[i] += p;
    }
}

