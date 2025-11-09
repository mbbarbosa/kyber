#ifndef REJ_SAMPLING_H
#define REJ_SAMPLING_H

#include <stddef.h>
#include <stdint.h>
#include "params.h"
#include "polyvec.h"


#define gen_vector KYBER_NAMESPACE(gen_vector)
void gen_vector(polyvec *v, const uint8_t seed[KYBER_SYMBYTES]);

#define gen_a(A,B)  gen_matrix(A,B,0)
#define gen_at(A,B) gen_matrix(A,B,1)

#define gen_matrix KYBER_NAMESPACE(gen_matrix)
void gen_matrix(polyvec *a, const uint8_t seed[KYBER_SYMBYTES], int transposed);


#endif
