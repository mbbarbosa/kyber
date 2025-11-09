#ifndef REJ_SAMPLING_H
#define REJ_SAMPLING_H

#include <stddef.h>
#include <stdint.h>
#include "params.h"
#include "polyvec.h"

#define rej_uniform KYBER_NAMESPACE(rej_uniform)
unsigned int rej_uniform(int16_t *r, unsigned int len, const uint8_t *buf, unsigned int buflen);

#endif
