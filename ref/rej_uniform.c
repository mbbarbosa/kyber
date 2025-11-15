#include <stdint.h>
#include "params.h"
#include "rej_uniform.h"
#include "poly.h"
#include "polyvec.h"
#include "symmetric.h"

#if(XOF_BLOCKBYTES % 3)
#error "Implementation of gen_matrix assumes that XOF_BLOCKBYTES is a multiple of 3"
#endif

//#define TEMPO_ALGORITHM 1
#ifndef TEMPO_ALGORITHM


/*************************************************
* Name:        rej_uniform
*
* Description: Run rejection sampling on uniform random bytes to generate
*              uniform random integers mod q
*
* Arguments:   - int16_t *r: pointer to output buffer
*              - unsigned int len: requested number of 16-bit integers (uniform mod q)
*              - const uint8_t *buf: pointer to input buffer (assumed to be uniformly random bytes)
*              - unsigned int buflen: length of input buffer in bytes
*
* Returns number of sampled 16-bit integers (at most len)
**************************************************/
static unsigned int rej_uniform(int16_t *r,
                                unsigned int len,
                                const uint8_t *buf,
                                unsigned int buflen)
{
  unsigned int ctr, pos;
  uint16_t val0, val1;

  ctr = pos = 0;
  while(ctr < len && pos + 3 <= buflen) {
    val0 = ((buf[pos+0] >> 0) | ((uint16_t)buf[pos+1] << 8)) & 0xFFF;
    val1 = ((buf[pos+1] >> 4) | ((uint16_t)buf[pos+2] << 4)) & 0xFFF;
    pos += 3;

    if(val0 < KYBER_Q)
      r[ctr++] = val0;
    if(ctr < len && val1 < KYBER_Q)
      r[ctr++] = val1;
  }

  return ctr;
}

/*************************************************
* Name:        gen_poly
*
* Description: Deterministically generate polynomial p from a seed. 
*              Performs rejection sampling on output of a XOF
*
* Arguments:   - poly *a: pointer to output poly p
*              - const xof_state *state: pointer to XOF state after absorb

**************************************************/
#define GEN_MATRIX_NBLOCKS ((12*KYBER_N/8*(1 << 12)/KYBER_Q + XOF_BLOCKBYTES)/XOF_BLOCKBYTES)

static void gen_poly(int16_t a[KYBER_N], xof_state *state)
{
  unsigned int ctr, k;
  unsigned int buflen, off;
  uint8_t buf[GEN_MATRIX_NBLOCKS*XOF_BLOCKBYTES+2];

  xof_squeezeblocks(buf, GEN_MATRIX_NBLOCKS, state);
  buflen = GEN_MATRIX_NBLOCKS*XOF_BLOCKBYTES;
  ctr = rej_uniform(a, KYBER_N, buf, buflen);

  while(ctr < KYBER_N) {
    off = buflen % 3;
    for(k = 0; k < off; k++) buf[k] = buf[buflen - off + k];
    xof_squeezeblocks(buf + off, 1, state);
    buflen = off + XOF_BLOCKBYTES;
    ctr += rej_uniform(a + ctr, KYBER_N - ctr, buf, buflen);
  }
}

#elif defined(TEMPO_ALGORITHM) && (TEMPO_ALGORITHM == 1)

/*************************************************
* Name:        rej_uniform
*
* Description: Run rejection sampling on uniform random bytes to generate
*              uniform random integers mod q
*
* Arguments:   - int16_t *r: pointer to output buffer
*              - unsigned int ctr: number of accepted 16-bit integers (uniform mod q)
*              - const uint8_t *buf: pointer to input buffer (assumed to be uniformly random bytes)
*              - unsigned int buflen: length of input buffer in bytes
*
* Returns number of sampled 16-bit integers (at most len)
**************************************************/

#define lt_1mask(x,y) (uint16_t)((((int16_t) x) - ((int16_t) y)) >> 15)      // 0xffffffff if  x < y, 0x00000000 otherwise
#define diff_1mask(x,y) (uint16_t)((0-((int16_t)((x^y) & 0x7fff))) >> 15)    // 0xffffffff if x != y, 0x00000000 otherwise
#define eq_1mask(x,y) (uint16_t)(~diff_1mask(x,y))                           // 0xffffffff if x == y, 0x00000000 otherwise

static unsigned int rej_uniform(int16_t *p, unsigned int ctr,
                                const uint8_t *buf,
                                unsigned int buflen)
{
  unsigned int pos;
  uint16_t d1, d2, acceptable_d1, acceptable_d2, match;

  pos = 0;
  while(pos + 3 <= buflen) {
    d1 = ((buf[pos+0] >> 0) | ((uint16_t)buf[pos+1] << 8)) & 0xFFF;
    d2 = ((buf[pos+1] >> 4) | ((uint16_t)buf[pos+2] << 4)) & 0xFFF;
    pos += 3;

    for (unsigned int j = 0; j < KYBER_N; j++) {
      match = eq_1mask(j,ctr);
      p[j] = (~match & p[j]) | (match & d1);
    }
    acceptable_d1 = lt_1mask(d1,KYBER_Q); 
    ctr += acceptable_d1 & 1;
    
    for (unsigned int j = 0; j < KYBER_N; j++) {
      match = eq_1mask(j,ctr);
      p[j] = (~match & p[j]) | (match & d2);
    }
    acceptable_d2 = lt_1mask(d2,KYBER_Q); 
    ctr += acceptable_d2 & 1;
  }

  return ctr;
}

/*************************************************
* Name:        gen_poly
*
* Description: Deterministically generate polynomial p from a seed. 
*              Performs rejection sampling on output of a XOF in constant
*              time by always writing to every position in the output poly
*
* Arguments:   - poly *a: pointer to output poly p
*              - const xof_state *state: pointer to XOF state after absorb

**************************************************/
#define MAX_ITER 8

static void gen_poly(int16_t a[KYBER_N], xof_state *state) {
  uint16_t ctr = 0;
  unsigned int k, buflen, off;
  uint8_t buf[XOF_BLOCKBYTES+2];

  xof_squeezeblocks(buf, 1, state);
  buflen = XOF_BLOCKBYTES;
  ctr = rej_uniform(a, ctr, buf, buflen);

  for (unsigned int iter = 0; iter < MAX_ITER; iter++) {
    off = buflen % 3;
    for(k = 0; k < off; k++) buf[k] = buf[buflen - off + k];
    xof_squeezeblocks(buf + off, 1, state);
    buflen = off + XOF_BLOCKBYTES;
    ctr += rej_uniform(a, ctr, buf, buflen);
  }
}

#elif defined(TEMPO_ALGORITHM) && (TEMPO_ALGORITHM == 2)

void gen_poly(uint16_t a[KYBER_N], const xof_state *state) {
    keccak_state ctx;
    uint8_t C[400]; // 400 bytes = 3200 bits

    shake128_init(&ctx);
    shake128_absorb(&ctx, seed, 34);
    shake128_finalize(&ctx);
    shake128_squeeze(C, sizeof(C), &ctx);

    // create OpenSSL BN context and numbers
    BN_CTX *bn_ctx = BN_CTX_new();
    BIGNUM *x_bn = BN_new();
    BIGNUM *m_bn = BN_new();
    BIGNUM *quotient = BN_new();
    BIGNUM *remainder = BN_new();

    // set m = KYBER_Q
    BN_set_word(m_bn, KYBER_Q);
    
    // set x_bn from the 400-byte buffer C 
    // (load little-endian, as other algorithms in this library)
    BN_lebin2bn(C, sizeof(C), x_bn);

    for (uint16_t j = 0; j < KYBER_N; j++) {
        // quotient = x_bn / m_bn
        // remainder = x_bn % m_bn
        BN_div(quotient, remainder, x_bn, m_bn, bn_ctx);

        // remainer fits in 16 bits
        a[j] = (uint16_t)BN_get_word(remainder);

        // x_bn = quotient for next iteration
        BN_copy(x_bn, quotient);
    }

    BN_free(x_bn);
    BN_free(m_bn);
    BN_free(quotient);
    BN_free(remainder);
    BN_CTX_free(bn_ctx);
}

#elif defined(TEMPO_ALGORITHM) && (TEMPO_ALGORITHM == 3a)

void algorithm3a(uint16_t a[KYBER_N], const uint8_t seed[34]) {
    keccak_state ctx;
    uint8_t C[24];

    shake128_init(&ctx);
    shake128_absorb(&ctx, seed, 34);
    shake128_finalize(&ctx);

    // create OpenSSL BN context and numbers
    BN_CTX *bn_ctx = BN_CTX_new();
    BIGNUM *x_bn = BN_new();        // 192-bit integer
    BIGNUM *m_bn = BN_new();        // modulus
    BIGNUM *r_bn = BN_new();        // remainder

    // set m = KYBER_Q
    BN_set_word(m_bn, KYBER_Q);

    for (int i = 0; i < KYBER_N; i++) {
        shake128_squeeze(C, 24, &ctx);

        // set x_bn from the 24-byte buffer C 
        // (load little-endian, as Algorithm 3B, to cross-check calculations)
        BN_lebin2bn(C, 24, x_bn);

        // compute remainder: r_bn = x_bn `mod` m_bn
        BN_mod(r_bn, x_bn, m_bn, bn_ctx);

        // get the remainder as a 16-bit integer
        a[i] = (uint16_t)BN_get_word(r_bn);
    }

    BN_free(x_bn);
    BN_free(m_bn);
    BN_free(r_bn);
    BN_CTX_free(bn_ctx);
}

#elif defined(TEMPO_ALGORITHM) && (TEMPO_ALGORITHM == 3b)

/*
TODO: check PQCLEAN_MLKEM512_CLEAN_montgomery_reduce(int32_t a) and
PQCLEAN_MLKEM512_CLEAN_barrett_reduce(int16_t a) implementations;
it might be possible to do this even faster
*/

// constants and pre-computations
const uint16_t m = 3329;
const uint32_t mu = 1290167;
const uint16_t pow2_mod_m[] = {
    1,    // 2^{i*32} = 2^0 mod 3329 (for x0)
    1353, // 2^{i*32} = 2^32 mod 3329 (for x1)
    2988, // 2^{i*32} = 2^64 mod 3329 (for x2)
    1358, // 2^{i*32} = 2^96 mod 3329 (for x3)
    3095, // 2^{i*32} = 2^128 mod 3329 (for x4)
    2982  // 2^{i*32} = 2^160 mod 3329 (for x5)
};

/*
Based on Barrett reduction described in HAC, Algorithm 14.42.

b = 2^16
k = 1
m = m0 \in [0,2^16) = KYBER_Q = 3329
x = (x1,x0) \in [0,2^32)
mu = floor(2^32/m) = 1290167
*/
uint16_t barrett_reduce(uint32_t x) {
    uint64_t q3;
    uint32_t r;

    q3 = ((uint64_t)x * mu) >> 32;
    r = x - (uint32_t)(q3 * (uint32_t)m);

    // if r < 0 then r = r + 2^32
    // r is declared as uint32_t, so it is *unsigned*
    // we check if r would have been negative by checking if the MSB is set
    uint32_t mask = -(r >= 0x80000000); // mask is 0xFFFFFFFF if r is negative, 0 otherwise
    r += (mask & 0x100000000ULL);

    // while r >= m do r = r - m, at most 2 loops
    // this follows from HAC, Fact 14.43
    // TODO: check if we can remove 1 loop, based on this result: https://blog.zksecurity.xyz/posts/barrett-tighter-bound/ 
    uint32_t mask1 = -(r >= m); // mask1 is 0xFFFFFFFF if r >= m, 0 otherwise
    r -= (mask1 & m);
    uint32_t mask2 = -(r >= m); // mask2 is 0xFFFFFFFF if r >= m, 0 otherwise
    r -= (mask2 & m);

    return (uint16_t)r;
}

/*
x = x5*2^160 + x4*2^128 + x3*2^92 + x2*2^64 + x1*2^32 + x0

r = x `mod` m
  = (x5*2^160 + x4*2^128 + x3*2^92 + x2*2^64 + x1*2^32 + x0) `mod` m
  = ((x5 `mod` m) * (2^160 `mod` m) + 
     (x4 `mod` m) * (2^128 `mod` m) + 
     (x3 `mod` m) * (2^92 `mod` m)  + 
     (x2 `mod` m) * (2^64 `mod` m)  + 
     (x1 `mod` m) * (2^32 `mod` m)  +
     (x0 `mod` m)) `mod` m
*/
uint16_t barrett_reduce192(uint32_t x[6]) {
    uint32_t temp = 0;

    for (int i = 5; i >= 0; i--)
        temp += (uint32_t)barrett_reduce(x[i]) * (uint32_t)pow2_mod_m[i];
    
    return barrett_reduce(temp);
}

void algorithm3b(uint16_t a[KYBER_N], const uint8_t seed[34]) {
    keccak_state ctx;
    uint8_t C[24];
    uint32_t x[6];

    shake128_init(&ctx);
    shake128_absorb(&ctx, seed, 34);
    shake128_finalize(&ctx);

    for (size_t j = 0; j < KYBER_N; j++) {
        shake128_squeeze(C, 24, &ctx);
        
        for (int i = 0; i < 6; i++) {
            x[i] = ( uint32_t)C[i * 4 + 0]        |
                   ((uint32_t)C[i * 4 + 1] << 8)  |
                   ((uint32_t)C[i * 4 + 2] << 16) |
                   ((uint32_t)C[i * 4 + 3] << 24);
        }
        a[j] = barrett_reduce192(x);
    }
}

#endif

/*************************************************
* Name:        gen_matrix
*
* Description: Deterministically generate matrix A (or the transpose of A)
*              from a seed. Entries of the matrix are polynomials that look
*              uniformly random. Performs rejection sampling on output of
*              a XOF
*
* Arguments:   - polyvec *a: pointer to ouptput matrix A
*              - const uint8_t *seed: pointer to input seed
*              - int transposed: boolean deciding whether A or A^T is generated
**************************************************/

// Not static for benchmarking
void gen_matrix(polyvec *a, const uint8_t seed[KYBER_SYMBYTES], int transposed)
{
  unsigned int i, j;
  xof_state state;

  for(i=0;i<KYBER_K;i++) {
    for(j=0;j<KYBER_K;j++) {
      if(transposed)
        xof_absorb(&state, seed, i, j);
      else
        xof_absorb(&state, seed, j, i);

      gen_poly(a[i].vec[j].coeffs,&state);
      }
    }
 }


/*************************************************
* Name:        gen_vector
*
* Description: Deterministically generate vector v from a seed. 
*              Entries of the vector are polynomials that look
*              uniformly random. Performs rejection sampling on 
*              output of a XOF
*
* Arguments:   - polyvec *a: pointer to output vector v
*              - const uint8_t *seed: pointer to input seed
**************************************************/

void gen_vector(polyvec *v, const uint8_t seed[KYBER_SYMBYTES])
{
  unsigned int i;
  xof_state state;

  for(i=0;i<KYBER_K;i++) {
    xof_absorb(&state, seed, i, 0); // take row 0
    gen_poly(v->vec[i].coeffs,&state);
  }
}

