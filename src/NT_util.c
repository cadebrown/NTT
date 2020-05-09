/* NT_util.c - Number Theory utilities */

#include "ntt.h"



// Returns the GCD (Greatest Common Denominator) of 'A' and 'B'
int64_t ntt_gcd(int64_t a, int64_t b) {
    if (b == 0) {
        return a;
    } else {
        return ntt_gcd(b, a % b);
    }
}

// Return the GCD (Greatest Common Denominator) of 'A' and 'B'
// AND calculate numbers 'x' and 'y' such that:
//   x * a + y * b = gcd(a, b)
// NOTE: xy[0] is x, xy[1] is y
int64_t ntt_egcd(int64_t a, int64_t b, int64_t* xy) {
	if (a == 0) {
        // base case
        xy[0] = 0;
        xy[1] = 1;
		return b;
	} else {
        // calculate a new 'xy' from reduction
        int64_t new_xy[2];
        int64_t gcd = ntt_egcd(b % a, a, new_xy);

        xy[0] = new_xy[1] - (b / a) * new_xy[0];
        xy[1] = new_xy[0];

        return gcd;
    }
}


// Return a^-1 (mod N)
// Or, return '0' if there was no inverse (i.e. it was not invertible)
int64_t ntt_modinv(int64_t a, int64_t N) {
    int64_t xy[2];
    int64_t gcd = ntt_egcd(a, N, xy);
    if (gcd != 1) {
        // not invertible!
        return 0;
    } else {
        int64_t r = xy[0] % N;
        if (r < 0) r += N;
        return r;
    }
}


// Calculate a*b (mod m)
uint64_t ntt_modmul(uint64_t a, uint64_t b, uint64_t m) {
    int64_t res = 0;
    while (a != 0) {
        if (a & 1) res = (res + b) % m;
        a >>= 1;
        b = (b << 1) % m;
    }
    return res;
}

// Calculate a^b (mod m)
int64_t ntt_modpow(int64_t a, int64_t b, int64_t m) {
    if (b < 0) {
        int64_t a_inv = ntt_modinv(a, m);
        // no modular inverse
        if (a_inv == 0) return -1;

        // return it using positive branch
        return ntt_modpow(a_inv, -b, m);
    } else {
        // convert to the right range
        a %= m;
        if (a < 0) a += m;

        // result
        int64_t res = 1;

        // now, calculate using repeated squaring
        while (b > 0) {
            /*if (m == 4309503329) {
                printf("%lli::%lli,%lli\n", res, a, b);
            }*/
            if (b % 2 == 1) {
                // multiply by active bit
                res = ntt_modmul(res, a, m);
            }

            // repeated squaring
            a = ntt_modmul(a, a, m);
            b /= 2;
        }

        // ensure result is positive
        if (res < 0) res += m;
        return res;
    }
}

/// Internal miller rabin trial test
static bool i_milrab(int64_t n, int64_t a) {


    if (n % a == 0) return false;

    // Decompose: n = 2^r * d + 1
    int64_t r = 0;
    int64_t d = n - 1;

    // take out powers of 2
    while (d % 2 == 0) {
        r++;
        d /= 2;
    }

    // calculate a^d (mod n)
    int64_t x = ntt_modpow(a, d, n);

    if (x == 1 || x == n - 1) {
        // still might be prime
        return true;
    } else {
        int64_t i;
        // complete the test 'r-1' times
        for (i = 0; i < r - 1; ++i) {
            x = ntt_modmul(x, x, n);
            // still might be prime

            if (x == n - 1) return true;
        }
        // not prime, definitely composite
        return false;
    }
}

// Return 0 if 'val' is composite, non-zero if 'val' is prime
bool ntt_isprime(int64_t val) {
    // special cases
    if (val < 2) return 0;
    if (val < 4 || val == 5) return 1;
    if (val % 2 == 0 || val % 3 == 0) return 0;

    /**/ if (val < 2047) return i_milrab(val, 2);
    else if (val < 1373653) return i_milrab(val, 2) && i_milrab(val, 3);
    else if (val < 9080191) return i_milrab(val, 31) && i_milrab(val, 73);
    else if (val < 25326001) return i_milrab(val, 2) && i_milrab(val, 3) && i_milrab(val, 5);
    else if (val < 3215031751UL) return i_milrab(val, 2) && i_milrab(val, 3) && i_milrab(val, 5) && i_milrab(val, 7);
    else if (val < 4759123141UL) return i_milrab(val, 2) && i_milrab(val, 7) && i_milrab(val, 61);
    else if (val < 1122004669633UL) return i_milrab(val, 2) && i_milrab(val, 13) && i_milrab(val, 23) && i_milrab(val, 1662803);
    else if (val < 2152302898747UL) return i_milrab(val, 2) && i_milrab(val, 3) && i_milrab(val, 5) && i_milrab(val, 7) && i_milrab(val, 11);
    else if (val < 3474749660383UL) return i_milrab(val, 2) && i_milrab(val, 3) && i_milrab(val, 5) && i_milrab(val, 7) && i_milrab(val, 11) && i_milrab(val, 13);
    else if (val < 341550071728321UL) return i_milrab(val, 2) && i_milrab(val, 3) && i_milrab(val, 5) && i_milrab(val, 7) && i_milrab(val, 11) && i_milrab(val, 13) && i_milrab(val, 17);
    else {
        // assume not prime
        return 0;
    }
}

// Compute Euler's totient function of 'n'
int64_t ntt_tot(int64_t n) {
	int64_t tot = n, i;

	for (i = 2; i * i <= n; i += 2) {
		if (n % i == 0) {
			while (n % i == 0) n /= i;
			tot -= tot / i;
		}

        // special case for first even number
        if (i == 2) i = 1;

	}
 
	if (n > 1) tot -= tot / n;
	return tot;
}


// Factor(n), yielding only unique prime factors
int64_t ntt_factor_uup(int64_t n, int64_t** facts) {
    int64_t nfacs = 0;

    // macro to add a factor
    #define ADD_FAC(_x) { \
        nfacs++; \
        *facts = realloc(*facts, (nfacs) * sizeof(**facts)); \
        (*facts)[nfacs - 1] = (_x); \
    }

    int64_t oN = n;

    if (n % 2 == 0) ADD_FAC(2);
    if (n % 3 == 0) ADD_FAC(3);

    // max to check
    int64_t mxN = n / 2;
    //int64_t mxN = (int64_t)floor(sqrt(n));

    int64_t i = 5;
    while (i <= mxN) {

        if (n % i == 0 && ntt_isprime(i)) {
            ADD_FAC(i);
            do {
                n /= i;
            } while (n % i == 0);
        }

        if (n % (i + 2) == 0 && ntt_isprime(i + 2)) {
            ADD_FAC(i + 2);
            do {
                n /= i + 2;
            } while (n % (i + 2) == 0);
        }

        i += 6;
    }

    // add remaining factor
    if (n != 1) ADD_FAC(n);

    #undef ADD_FAC
    return nfacs;
}


// Return the first primitive root of unity (mod n), or 0 if none exists
int64_t ntt_prim_root_unity(int64_t n) {

    // totient(p)
    int64_t tot = ntt_tot(n);

    // totient factors
    int64_t* tot_facts = NULL;
    int64_t tot_n_facts = ntt_factor_uup(tot, &tot_facts);
    
    // keep testing out tries
    int64_t a = 2;
    while (a <= n) {
        // whether it equalled one
        int64_t eq1 = 0, i;
        for (i = 0; i < tot_n_facts; ++i) {
            // calculate: a^(phi(n)/p_i) (mod p)
            int64_t r = ntt_modpow(a, tot / tot_facts[i], n);
            //printf("%i ^ %i = %i (mod %i)\n", a, tot / tot_facts[i], r, n);
            if (r == 1) {
                eq1 = 1;
                break;
            }
        }

        if (eq1 == 0) {
            free(tot_facts);
            return a;
        }
        a++;
    }

    // none
    free(tot_facts);
    return 0;
}

// Calculate n'th root of unity mod p
int64_t ntt_nth_root_unity(int64_t n, int64_t p) {
    int64_t i;
    for (i = 2; i < n; ++i) {
        if (ntt_modpow(i, n, p) == 1) {
            return i;
        }
    }
    // no nth root of unity found!
    return 0;

}
