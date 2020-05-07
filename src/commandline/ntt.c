

#include "ntt.h"

#include <time.h>
#include <sys/time.h>

// get the start time (initialize it in 'ks_init')
static struct timeval ntt_start_time = (struct timeval){ .tv_sec = 0, .tv_usec = 0 };

// return the time since it started
double ntt_time() {
    struct timeval curtime;
    gettimeofday(&curtime, NULL);
    return (curtime.tv_sec - ntt_start_time.tv_sec) + 1.0e-6 * (curtime.tv_usec - ntt_start_time.tv_usec);
}


// print an array of items
void printarr(int64_t* N, int n) {
    int i;
    for (i = 0; i < n; ++i) {
        if (i != 0) printf(" ");
        printf("%3li", N[i]);
    }
}

// print an 'nxn' matrix
void printmat(int64_t* N, int n) {
    int i, j;
    for (i = 0; i < n; ++i) {
        printf("[");
        for (j = 0; j < n; ++j) {
            if (j != 0) printf(" ");
            printf("%3li", N[i * n + j]);
        }
        printf("]\n");
    }
}


int main(int argc, char** argv) {
    gettimeofday(&ntt_start_time, NULL);
    srand(time(NULL));

    // temp. vars
    int i, j;

    // transform sizes
    int64_t N = 1024 * 1024;

    // find suitable prime
    int64_t p = N + 1;
    while (!ntt_isprime(p)) {
        p += N;
    }

    // these are good sizes

    // generate NTT plan
    ntt_plan_gemm_t plan_G = NTT_PLAN_GEMM_EMPTY;
    ntt_plan_gemm_init(&plan_G, N, p);

    ntt_plan_bfly_t plan_B = NTT_PLAN_BFLY_EMPTY;
    ntt_plan_bfly_init(&plan_B, N, p);

    int64_t* A = malloc(sizeof(*A) * N);
    int64_t* B = malloc(sizeof(*B) * N);
    int64_t* C = malloc(sizeof(*C) * N);

    int64_t* nttA = malloc(sizeof(*nttA) * N);
    int64_t* nttB = malloc(sizeof(*nttB) * N);
    int64_t* nttC = malloc(sizeof(*nttC) * N);

    for (i = 0; i < N/2; ++i) {
        A[i] = rand() % p;
    }

    for (i = 0; i < N/2; ++i) {
        B[i] = (i == 1) ? 1 : 0;
    }

    for (; i < N; ++i) A[i] = B[i] = 0;


    printf("A       : ");
    printarr(A, N);
    printf("\n");

    printf("B       : ");
    printarr(B, N);
    printf("\n");

    ntt_plan_bfly_NTT(&plan_B, A, nttA);
    ntt_plan_bfly_NTT(&plan_B, B, nttB);

    printf("ntt(A)  : ");
    printarr(nttA, N);
    printf("\n");
    printf("ntt(B)  : ");
    printarr(nttB, N);
    printf("\n");

    for (i = 0; i < N; ++i) {
        nttC[i] = (nttA[i] * nttB[i]) % plan_G.p;
    }

    printf("ntt(C)  : ");
    printarr(nttC, N);
    printf("\n");

    ntt_plan_bfly_INTT(&plan_B, nttC, C);

    printf("C       : ");
    printarr(C, N);
    printf("\n");

    /*
    for (i = 0 ; i < N; ++i) {
        if (C[i] != A[i]) {
            printf("ER\n");
        }
    }*/

}
