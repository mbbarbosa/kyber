#ifndef REJ_SAMPLING_H
#define REJ_SAMPLING_H

#include <stddef.h>
#include <stdint.h>
#include "params.h"
#include "polyvec.h"

#ifndef TEMPO_VECTOR_ALG
#define gen_poly_vector gen_poly_original
#elif (defined (TEMPO_VECTOR_ALG) && (TEMPO_VECTOR_ALG == 1))
#define gen_poly_vector gen_poly_tmp1
#elif (defined (TEMPO_VECTOR_ALG) && (TEMPO_VECTOR_ALG == 2))
#define gen_poly_vector gen_poly_tmp2
#elif (defined (TEMPO_VECTOR_ALG) && (TEMPO_VECTOR_ALG == 3))
#define gen_poly_vector gen_poly_tmp3a
#elif (defined (TEMPO_VECTOR_ALG) && (TEMPO_VECTOR_ALG == 4))
#define gen_poly_vector gen_poly_tmp3b
#endif

#ifndef TEMPO_MATRIX_ALG
#define gen_poly_matrix gen_poly_original
#elif (defined (TEMPO_MATRIX_ALG) && (TEMPO_MATRIX_ALG == 1))
#define gen_poly_matrix gen_poly_tmp1
#elif (defined (TEMPO_MATRIX_ALG) && (TEMPO_MATRIX_ALG == 2))
#define gen_poly_matrix gen_poly_tmp2
#elif (defined (TEMPO_MATRIX_ALG) && (TEMPO_MATRIX_ALG == 3))
#define gen_poly_matrix gen_poly_tmp3a
#elif (defined (TEMPO_MATRIX_ALG) && (TEMPO_MATRIX_ALG == 4))
#define gen_poly_matrix gen_poly_tmp3b
#endif

#define gen_vector KYBER_NAMESPACE(gen_vector)
void gen_vector(polyvec *v, const uint8_t seed[KYBER_SYMBYTES]);

#define gen_a(A,B)  gen_matrix(A,B,0)
#define gen_at(A,B) gen_matrix(A,B,1)

#define gen_matrix KYBER_NAMESPACE(gen_matrix)
void gen_matrix(polyvec *a, const uint8_t seed[KYBER_SYMBYTES], int transposed);


#endif
