

#include "ntt.h"

#include <time.h>
#include <sys/time.h>

// get the start time (initialize it in 'ks_init')
static struct timeval ntt_start_time = (struct timeval){ .tv_sec = 0, .tv_usec = 0 };

// return the time since it started
static double ntt_time() {
    struct timeval curtime;
    gettimeofday(&curtime, NULL);
    return (curtime.tv_sec - ntt_start_time.tv_sec) + 1.0e-6 * (curtime.tv_usec - ntt_start_time.tv_usec);
}

// print an array of items
static void printarr(int64_t* N, int n) {
    int i;
    printf("[");
    for (i = 0; i < n; ++i) {
        if (i != 0) printf(", ");
        printf("%lli", (long long int)N[i]);
    }
    printf("]");
}

// print an 'nxn' matrix
static void printmat(int64_t* N, int n) {
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

// get a hex digit from a character
static int gethexdig(char c) {
    /**/ if (c >= '0' && c <= '9') return c - '0';
    else if (c >= 'a' && c <= 'z') return 10 + c - 'a';
    else if (c >= 'A' && c <= 'Z') return 10 + c - 'A';
    return 0;
}

// digit string used for base conversion
static const char digstr[] = "0123456789abcdefghijklmnopqrstuvwxyz";

// get character from digit
static char gethexchar(int dig) {
    return digstr[dig];
}

// read a sequence, return '0' if there was an error
static int64_t readseq(char* fname, int64_t** to) {
    FILE* fp = fopen(fname, "r");
    if (fp == NULL) return 0;

    int n = 0;
    long long int ci;
    while (fscanf(fp, "%lli", &ci) == 1) {

        // append
        n++;
        *to = realloc(*to, sizeof(**to) * n);
        (*to)[n - 1] = ci;
    }


    fclose(fp);
    return n;
}

// read hex integer (which can be from a file, or provided in literal form in 'src')
// hdpw -> the number of hex digits per word to cram into 'to'
static int64_t readhexint(char* src, int64_t** to, int hdpw) {
    char* ext = strrchr(src, '.');
    char* hexdata = NULL;
    char* tofree = NULL;
    if (ext != NULL) {
        // read file
        FILE* fp = fopen(src, "r");
        if (fp == NULL) {
            return 0;
        }

        fseek(fp, 0, SEEK_END);
        int64_t fsize = ftell(fp);
        tofree = hexdata = malloc(fsize + 1);
        fseek(fp, 0, SEEK_SET);
        
        fread(hexdata, 1, fsize, fp);
        hexdata[fsize] = '\0';

        fclose(fp);
    } else {
        // treat it as from the source
        hexdata = src;
    }

    // remove leading signifiers
    if (hexdata[0] == '0' && hexdata[1] == 'x') hexdata += 2;

    int64_t hds = strlen(hexdata);

    // calculate size requirements
    int64_t npts = hds / hdpw + 1;
    *to = realloc(*to, sizeof(**to) * npts);

    // fill in as zeros
    int64_t i = 0, j; 
    while (i < npts) (*to)[i++] = 0;

    // temporary string
    char tmps[256];

    int64_t oi;
    for (oi = 0, i = hds - 1; oi < npts && i >= 0; ++oi) {
        // calculate current output word
        int64_t word = 0;

        int tmpsi = 0;

        // try sub-digits of the word
        for (j = 0; j < hdpw; ++j) {

            // start searching for next valid char
            while (i >= 0 && (hexdata[i] == ' ' || hexdata[i] == '\t' || hexdata[i] == '\n')) i--;

            // reached the beginning of the input
            if (i < 0) break;

            // add digit to string
            tmps[tmpsi++] = hexdata[i--];

        }

        // Nul-terminate
        tmps[tmpsi] = '\0';

        // add digits to binary form
        for (j = tmpsi - 1; j >= 0; --j) {
            // append digit
            word = word * 16 + gethexdig(tmps[j]);
        }

        (*to)[oi] = word;
    }

    if (oi < 0) oi++;

    // now, read them in
    free(tofree);
    return oi;
}




int main(int argc, char** argv) {
    gettimeofday(&ntt_start_time, NULL);
    srand(time(NULL));

    if (argc == 1 || (argc > 1 && strcmp(argv[1], "help") == 0)) {
        fprintf(stderr, "Usage: ntt [cmd] args...\n");
        fprintf(stderr, " [cmd]:\n");
        fprintf(stderr, "   help:                  prints this help message\n");
        fprintf(stderr, "   ntt [file] [p=0]:      calculates the NTT of an sequence of integers from a file (optional modulus p)\n");
        fprintf(stderr, "   intt [file] [p=0]:     calculates the INTT of an sequence of integers from a file (optional modulus p)\n");
        fprintf(stderr, "   mul [A] [B]            Uses 'NTT' to calculate A*B\n");
        
        return 1;
    }


    char* cmd = argv[1];
    if (strcmp(cmd, "ntt") == 0) {
        if (argc < 3 || argc > 4) {
            fprintf(stderr, "Expected it to be 'ntt ntt [file] [p=0]'\n");
            return 1;
        }

        int64_t p = 0;

        int64_t* x = NULL;
        int64_t N = readseq(argv[2], &x);
        if (N == 0) {
            fprintf(stderr, "Could not open '%s'\n", argv[2]);
            return 1;
        }

        if ((N & (N - 1)) != 0) {
            fprintf(stderr, "Sequence length must be a power of 2! (got N=%lli)\n", (long long int)N);
            free(x);
            return 1;
        }

        if (argc >= 4) {
            long long int p_read = 0;
            sscanf(argv[3], "%lli", &p_read);
            p = p_read;
            if (!ntt_isprime(p) || p % N != 1) {
                fprintf(stderr, "Invalid choice 'p' (given %lli) for N=%lli\n", p_read, N);
                free(x);
                return 1;
            }
        } else {
            // generate it
            p = N + 1;
            while (!ntt_isprime(p)) {
                p += N;
            }
        }

        // now, calculate plan
        ntt_plan_bfly_t plan_B = NTT_PLAN_BFLY_EMPTY;
        ntt_plan_bfly_init(&plan_B, N, p);

        int64_t* ntt_x = malloc(sizeof(*ntt_x) * N);

        ntt_plan_bfly_NTT(&plan_B, x, ntt_x);

        // print it out
        int i;
        for (i = 0; i < N; ++i) {
            if (i != 0) printf(" ");
            printf("%lli", (long long int)(ntt_x[i]));
        }

        printf("\n");

        free (x);
        free (ntt_x);

    } else if (strcmp(cmd, "intt") == 0) {
        // inverse transform
        if (argc < 3 || argc > 4) {
            fprintf(stderr, "Expected it to be 'ntt intt [file] [p=0]'\n");
            return 1;
        }

        int64_t p = 0;

        int64_t* x = NULL;
        int64_t N = readseq(argv[2], &x);
        if (N == 0) {
            fprintf(stderr, "Could not open '%s'\n", argv[2]);
            return 1;
        }

        if ((N & (N - 1)) != 0) {
            fprintf(stderr, "Sequence length must be a power of 2! (got N=%lli)\n", (long long int)N);
            free(x);
            return 1;
        }

        if (argc >= 4) {
            long long int p_read = 0;
            sscanf(argv[3], "%lli", &p_read);
            p = p_read;
            if (!ntt_isprime(p) || p % N != 1) {
                fprintf(stderr, "Invalid choice 'p' (given %lli) for N=%lli\n", p_read, N);
                free(x);
                return 1;
            }
        } else {
        }

        // now, calculate plan
        ntt_plan_bfly_t plan_B = NTT_PLAN_BFLY_EMPTY;
        ntt_plan_bfly_init(&plan_B, N, p);

        int64_t* ntt_x = malloc(sizeof(*ntt_x) * N);

        ntt_plan_bfly_INTT(&plan_B, x, ntt_x);

        // print it out
        int i;
        for (i = 0; i < N; ++i) {
            if (i != 0) printf(" ");
            printf("%lli", (long long int)(ntt_x[i]));
        }

        printf("\n");

        free (x);
        free (ntt_x);


    } else if (strcmp(cmd, "mulhex") == 0) {
        // read 2 hex numbers and multiply them

        // hex digits per word
        int hdpw = 2;

        // tmp vars
        int64_t i, j;

        // our main variables, 'A' and 'B'
        int64_t* A = NULL, *B = NULL;

        int64_t nA = readhexint(argv[2], &A, hdpw);
        if (nA < 1) {
            fprintf(stderr, "Could not get hex int from '%s'\n", argv[2]);
            return 1;
        }
        int64_t nB = readhexint(argv[3], &B, hdpw);
        if (nB < 1) {
            fprintf(stderr, "Could not get hex int from '%s'\n", argv[3]);
            return 1;
        }

        // calculate the required transform size (rounded up to the next power of 2)
        int64_t N = 1 << 2;
        while (N < (nA + nB + 1)) {
            N = N * 2;
        }

        // maximum limb residual
        int64_t limb_res = 1ULL << (4 * hdpw);

        // now, pad 'A' and 'B'
        A = realloc(A, sizeof(*A) * N);
        B = realloc(B, sizeof(*B) * N);

        // pad with 0s
        for (i = nA; i < N; ++i) A[i] = 0;
        for (i = nB; i < N; ++i) B[i] = 0;

        // output variable
        int64_t* C = malloc(sizeof(*C) * N);

        // calculate the multiplication utility
        ntt_multer_t multer = NTT_MULTER_EMPTY;
        ntt_multer_init(&multer, N);


        /*
        printf("A: ");
        printarr(A, N);
        printf("\n");

        printf("B: ");
        printarr(B, N);
        printf("\n");
        */

        /*
        // print the transform size
        printf("N = %lli\n", multer.N);

        // debug the transform modulus(es)
        printf("prod(p[:]) = %lli = ", multer.prod_p);
        for (i = 0; i < multer.n_plans; ++i) {
            if (i > 0) printf(" * ");
            printf("%lli", multer.plans[i].p);
        }
        printf("\n");
        */
        /* computation */

        double st = ntt_time();

        // perform multiplication
        ntt_multer_mult(&multer, A, B, C);

        st = ntt_time() - st;
        fprintf(stderr, "time: %.3lf\n", st);

        // handle carry propogation
        int64_t carry = 0;
        for (i = 0; i < N; ++i) {
            int64_t csum = C[i] + carry;
            carry = csum / limb_res;
            csum %= limb_res;
            C[i] = csum;
        }

        /*
        printf("C: ");
        printarr(C, N);
        printf("\n");
        */

        char* Cs = malloc(hdpw * N + 1);

        int out_i = 0;

        // temporary string
        char tmps[256];

        // output hex string
        for (i = 0; i < N; ++i) {

            uint64_t ci = C[i];

            for (j = 0; j < hdpw; ++j) {
                // extract digits per word
                char cc = gethexchar(ci % 16);
                Cs[out_i++] = cc;

                ci /= 16;
            }
        }

        Cs[out_i] = '\0';

        // remove extra white space
        while (out_i > 1 && Cs[out_i - 1] == '0') {
            out_i--;
        }

        // reverse the string
        for (i = 0; 2 * i < out_i; ++i) {
            char t = Cs[i];
            Cs[i] = Cs[out_i - i - 1];
            Cs[out_i - i - 1] = t;
        }

        Cs[out_i] = '\0';
        
        printf("0x%s\n", Cs);
        free(Cs);

    } else {
        fprintf(stderr, "Error! Invalid cmd, run `ntt help` for help\n");
        return 1;
    }

/*
    // temp. vars
    int i, j;

    // transform sizes
    int64_t N = 4096;

    // find suitable prime
    int64_t p = N * 1024 + 1;
    while (!ntt_isprime(p)) {
        p += N;
    }

    // these are good sizes

    // generate NTT plan
    ntt_plan_gemm_t plan_G = NTT_PLAN_GEMM_EMPTY;
    ntt_plan_gemm_init(&plan_G, N, p);

    ntt_plan_bfly_t plan_B = NTT_PLAN_BFLY_EMPTY;
    ntt_plan_bfly_init(&plan_B, N, p);

    printf(" -*-  PLANS -*-\n");
    printf("GEMM:\n");
    printf("  N: %li\n", plan_G.N);
    printf("  p: %li\n", plan_G.p);

    printf("BFLY:\n");
    printf("  N: %li\n", plan_B.N);
    printf("  p: %li\n", plan_B.p);



    int64_t* A = malloc(sizeof(*A) * N);
    int64_t* B = malloc(sizeof(*B) * N);
    int64_t* C = malloc(sizeof(*C) * N);

    int64_t* nttA = malloc(sizeof(*nttA) * N);
    int64_t* nttB = malloc(sizeof(*nttB) * N);
    int64_t* nttC = malloc(sizeof(*nttC) * N);

    // zero them out
    for (i = 0; i < N; ++i) A[i] = B[i] = 0;

    // now, set 'A' and 'B'
    A[1] = 2;
    B[0] = 91;
    B[2] = 93;

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


    // now, perform shift

    int64_t c = 0, m = 100;
    for (i = 0; i < N; ++i) {
        int64_t sm = C[i] + c;

        c = sm / m;
        sm %= m;
        C[i] = sm;
    }

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
