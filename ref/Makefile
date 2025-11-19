CC ?= /usr/bin/cc
CFLAGS += -Wall -Wextra -Wpedantic -Wmissing-prototypes -Wredundant-decls \
  -Wshadow -Wpointer-arith -O3 -fomit-frame-pointer -z noexecstack 
NISTFLAGS += -Wno-unused-result -O3 -fomit-frame-pointer
RM = /bin/rm

SOURCES = kem.c indcpa.c polyvec.c poly.c ntt.c cbd.c reduce.c verify.c rej_uniform.c 
SOURCESKECCAK = $(SOURCES) fips202.c symmetric-shake.c
HEADERS = params.h kem.h indcpa.h polyvec.h poly.h ntt.h cbd.h reduce.c verify.h symmetric.h rej_uniform.h
HEADERSKECCAK = $(HEADERS) fips202.h

.PHONY: all speed shared clean

all: test speed shared nistkat

test: \
  test/test_kyber512 \
  test/test_kyber768 \
  test/test_kyber1024 \
  test/test_vectors512 \
  test/test_vectors768 \
  test/test_vectors1024 \
  test/test_kyber512_tmp1 \
  test/test_kyber768_tmp1 \
  test/test_kyber1024_tmp1 \
  test/test_vectors512_tmp1 \
  test/test_vectors768_tmp1 \
  test/test_vectors1024_tmp1 \
  test/test_kyber512_tmp2 \
  test/test_kyber768_tmp2 \
  test/test_kyber1024_tmp2 \
  test/test_kyber512_tmp3a \
  test/test_kyber768_tmp3a \
  test/test_kyber1024_tmp3a \
  test/test_kyber512_tmp3b \
  test/test_kyber768_tmp3b \
  test/test_kyber1024_tmp3b \

speed: \
  test/test_speed512 \
  test/test_speed768 \
  test/test_speed1024 \
  test/test_speed512_tmp1 \
  test/test_speed768_tmp1 \
  test/test_speed1024_tmp1 \
  test/test_speed512_tmp2 \
  test/test_speed768_tmp2 \
  test/test_speed1024_tmp2 \
  test/test_speed512_tmp3a \
  test/test_speed768_tmp3a \
  test/test_speed1024_tmp3a \
  test/test_speed512_tmp3b \
  test/test_speed768_tmp3b \
  test/test_speed1024_tmp3b \

shared: \
  lib/libpqcrystals_kyber512_ref.so \
  lib/libpqcrystals_kyber768_ref.so \
  lib/libpqcrystals_kyber1024_ref.so \
  lib/libpqcrystals_fips202_ref.so \

nistkat: \
	nistkat/PQCgenKAT_kem512 \
	nistkat/PQCgenKAT_kem768 \
	nistkat/PQCgenKAT_kem1024 \
	nistkat/PQCgenKAT_kem512_tmp1 \
	nistkat/PQCgenKAT_kem768_tmp1 \
	nistkat/PQCgenKAT_kem1024_tmp1 \
	nistkat/PQCgenKAT_kem512_tmp3a \
	nistkat/PQCgenKAT_kem768_tmp3a \
	nistkat/PQCgenKAT_kem1024_tmp3a \
	nistkat/PQCgenKAT_kem512_tmp3b \
	nistkat/PQCgenKAT_kem768_tmp3b \
	nistkat/PQCgenKAT_kem1024_tmp3b \

lib/libpqcrystals_fips202_ref.so: fips202.c fips202.h
	mkdir -p lib
	$(CC) -shared -fPIC $(CFLAGS) fips202.c -o $@

lib/libpqcrystals_kyber512_ref.so: $(SOURCES) $(HEADERS) symmetric-shake.c
	mkdir -p lib
	$(CC) -shared -fPIC $(CFLAGS) -DKYBER_K=2 $(SOURCES) symmetric-shake.c -o $@

lib/libpqcrystals_kyber768_ref.so: $(SOURCES) $(HEADERS) symmetric-shake.c
	mkdir -p lib
	$(CC) -shared -fPIC $(CFLAGS) -DKYBER_K=3 $(SOURCES) symmetric-shake.c -o $@

lib/libpqcrystals_kyber1024_ref.so: $(SOURCES) $(HEADERS) symmetric-shake.c
	mkdir -p lib
	$(CC) -shared -fPIC $(CFLAGS) -DKYBER_K=4 $(SOURCES) symmetric-shake.c -o $@

test/test_kyber512: $(SOURCESKECCAK) $(HEADERSKECCAK) test/test_kyber.c randombytes.c
	$(CC) $(CFLAGS) -DKYBER_K=2 $(SOURCESKECCAK) randombytes.c test/test_kyber.c -o $@

test/test_kyber768: $(SOURCESKECCAK) $(HEADERSKECCAK) test/test_kyber.c randombytes.c
	$(CC) $(CFLAGS) -DKYBER_K=3 $(SOURCESKECCAK) randombytes.c test/test_kyber.c -o $@

test/test_kyber1024: $(SOURCESKECCAK) $(HEADERSKECCAK) test/test_kyber.c randombytes.c
	$(CC) $(CFLAGS) -DKYBER_K=4 $(SOURCESKECCAK) randombytes.c test/test_kyber.c -o $@

test/test_vectors512: $(SOURCESKECCAK) $(HEADERSKECCAK) test/test_vectors.c
	$(CC) $(CFLAGS) -DKYBER_K=2 $(SOURCESKECCAK) test/test_vectors.c -o $@

test/test_vectors768: $(SOURCESKECCAK) $(HEADERSKECCAK) test/test_vectors.c
	$(CC) $(CFLAGS) -DKYBER_K=3 $(SOURCESKECCAK) test/test_vectors.c -o $@

test/test_vectors1024: $(SOURCESKECCAK) $(HEADERSKECCAK) test/test_vectors.c
	$(CC) $(CFLAGS) -DKYBER_K=4 $(SOURCESKECCAK) test/test_vectors.c -o $@

test/test_speed512: $(SOURCESKECCAK) $(HEADERSKECCAK) test/cpucycles.h test/cpucycles.c test/speed_print.h test/speed_print.c test/test_speed.c randombytes.c
	$(CC) $(CFLAGS) -DKYBER_K=2 $(SOURCESKECCAK) randombytes.c test/cpucycles.c test/speed_print.c test/test_speed.c -o $@

test/test_speed768: $(SOURCESKECCAK) $(HEADERSKECCAK) test/cpucycles.h test/cpucycles.c test/speed_print.h test/speed_print.c test/test_speed.c randombytes.c
	$(CC) $(CFLAGS) -DKYBER_K=3 $(SOURCESKECCAK) randombytes.c test/cpucycles.c test/speed_print.c test/test_speed.c -o $@

test/test_speed1024: $(SOURCESKECCAK) $(HEADERSKECCAK) test/cpucycles.h test/cpucycles.c test/speed_print.h test/speed_print.c test/test_speed.c randombytes.c
	$(CC) $(CFLAGS) -DKYBER_K=4 $(SOURCESKECCAK) randombytes.c test/cpucycles.c test/speed_print.c test/test_speed.c -o $@

nistkat/PQCgenKAT_kem512: $(SOURCESKECCAK) $(HEADERSKECCAK) nistkat/PQCgenKAT_kem.c nistkat/rng.c nistkat/rng.h
	$(CC) $(NISTFLAGS) -DKYBER_K=2 -o $@ $(SOURCESKECCAK) nistkat/rng.c nistkat/PQCgenKAT_kem.c $(LDFLAGS) -lcrypto

nistkat/PQCgenKAT_kem768: $(SOURCESKECCAK) $(HEADERSKECCAK) nistkat/PQCgenKAT_kem.c nistkat/rng.c nistkat/rng.h
	$(CC) $(NISTFLAGS) -DKYBER_K=3 -o $@ $(SOURCESKECCAK) nistkat/rng.c nistkat/PQCgenKAT_kem.c $(LDFLAGS) -lcrypto

nistkat/PQCgenKAT_kem1024: $(SOURCESKECCAK) $(HEADERSKECCAK) nistkat/PQCgenKAT_kem.c nistkat/rng.c nistkat/rng.h
	$(CC) $(NISTFLAGS) -DKYBER_K=4 -o $@ $(SOURCESKECCAK) nistkat/rng.c nistkat/PQCgenKAT_kem.c $(LDFLAGS) -lcrypto

# TEMPO ALGORITHM1 VARIANTS

test/test_kyber512_tmp1: $(SOURCESKECCAK) $(HEADERSKECCAK) test/test_kyber.c randombytes.c
	$(CC) $(CFLAGS) -DTEMPO_MATRIX_ALG=1 -DKYBER_K=2 $(SOURCESKECCAK) randombytes.c test/test_kyber.c -o $@

test/test_kyber768_tmp1: $(SOURCESKECCAK) $(HEADERSKECCAK) test/test_kyber.c randombytes.c
	$(CC) $(CFLAGS) -DTEMPO_MATRIX_ALG=1 -DKYBER_K=3 $(SOURCESKECCAK) randombytes.c test/test_kyber.c -o $@

test/test_kyber1024_tmp1: $(SOURCESKECCAK) $(HEADERSKECCAK) test/test_kyber.c randombytes.c
	$(CC) $(CFLAGS) -DTEMPO_MATRIX_ALG=1 -DKYBER_K=4 $(SOURCESKECCAK) randombytes.c test/test_kyber.c -o $@

test/test_vectors512_tmp1: $(SOURCESKECCAK) $(HEADERSKECCAK) test/test_vectors.c
	$(CC) $(CFLAGS) -DTEMPO_MATRIX_ALG=1 -DKYBER_K=2 $(SOURCESKECCAK) test/test_vectors.c -o $@

test/test_vectors768_tmp1: $(SOURCESKECCAK) $(HEADERSKECCAK) test/test_vectors.c
	$(CC) $(CFLAGS) -DTEMPO_MATRIX_ALG=1 -DKYBER_K=3 $(SOURCESKECCAK) test/test_vectors.c -o $@

test/test_vectors1024_tmp1: $(SOURCESKECCAK) $(HEADERSKECCAK) test/test_vectors.c
	$(CC) $(CFLAGS) -DTEMPO_MATRIX_ALG=1 -DKYBER_K=4 $(SOURCESKECCAK) test/test_vectors.c -o $@

test/test_speed512_tmp1: $(SOURCESKECCAK) $(HEADERSKECCAK) test/cpucycles.h test/cpucycles.c test/speed_print.h test/speed_print.c test/test_speed.c randombytes.c
	$(CC) $(CFLAGS) -DTEMPO_MATRIX_ALG=1 -DKYBER_K=2 $(SOURCESKECCAK) randombytes.c test/cpucycles.c test/speed_print.c test/test_speed.c -o $@

test/test_speed768_tmp1: $(SOURCESKECCAK) $(HEADERSKECCAK) test/cpucycles.h test/cpucycles.c test/speed_print.h test/speed_print.c test/test_speed.c randombytes.c
	$(CC) $(CFLAGS) -DTEMPO_MATRIX_ALG=1 -DKYBER_K=3 $(SOURCESKECCAK) randombytes.c test/cpucycles.c test/speed_print.c test/test_speed.c -o $@

test/test_speed1024_tmp1: $(SOURCESKECCAK) $(HEADERSKECCAK) test/cpucycles.h test/cpucycles.c test/speed_print.h test/speed_print.c test/test_speed.c randombytes.c
	$(CC) $(CFLAGS) -DTEMPO_MATRIX_ALG=1 -DKYBER_K=4 $(SOURCESKECCAK) randombytes.c test/cpucycles.c test/speed_print.c test/test_speed.c -o $@

nistkat/PQCgenKAT_kem512_tmp1: $(SOURCESKECCAK) $(HEADERSKECCAK) nistkat/PQCgenKAT_kem.c nistkat/rng.c nistkat/rng.h
	$(CC) $(NISTFLAGS) -DTEMPO_MATRIX_ALG=1 -DKYBER_K=2 -o $@ $(SOURCESKECCAK) nistkat/rng.c nistkat/PQCgenKAT_kem.c $(LDFLAGS) -lcrypto

nistkat/PQCgenKAT_kem768_tmp1: $(SOURCESKECCAK) $(HEADERSKECCAK) nistkat/PQCgenKAT_kem.c nistkat/rng.c nistkat/rng.h
	$(CC) $(NISTFLAGS) -DTEMPO_MATRIX_ALG=1 -DKYBER_K=3 -o $@ $(SOURCESKECCAK) nistkat/rng.c nistkat/PQCgenKAT_kem.c $(LDFLAGS) -lcrypto

nistkat/PQCgenKAT_kem1024_tmp1: $(SOURCESKECCAK) $(HEADERSKECCAK) nistkat/PQCgenKAT_kem.c nistkat/rng.c nistkat/rng.h
	$(CC) $(NISTFLAGS) -DTEMPO_MATRIX_ALG=1 -DKYBER_K=4 -o $@ $(SOURCESKECCAK) nistkat/rng.c nistkat/PQCgenKAT_kem.c $(LDFLAGS) -lcrypto

# TEMPO ALGORITHM2 VARIANTS

test/test_kyber512_tmp2: $(SOURCESKECCAK) $(HEADERSKECCAK) test/test_kyber.c randombytes.c
	$(CC) $(CFLAGS) -DTEMPO_MATRIX_ALG=2 -DKYBER_K=2 $(SOURCESKECCAK) randombytes.c test/test_kyber.c -lcrypto -o $@

test/test_kyber768_tmp2: $(SOURCESKECCAK) $(HEADERSKECCAK) test/test_kyber.c randombytes.c
	$(CC) $(CFLAGS) -DTEMPO_MATRIX_ALG=2 -DKYBER_K=3 $(SOURCESKECCAK) randombytes.c test/test_kyber.c -lcrypto -o $@

test/test_kyber1024_tmp2: $(SOURCESKECCAK) $(HEADERSKECCAK) test/test_kyber.c randombytes.c
	$(CC) $(CFLAGS) -DTEMPO_MATRIX_ALG=2 -DKYBER_K=4 $(SOURCESKECCAK) randombytes.c test/test_kyber.c -lcrypto -o $@

test/test_speed512_tmp2: $(SOURCESKECCAK) $(HEADERSKECCAK) test/cpucycles.h test/cpucycles.c test/speed_print.h test/speed_print.c test/test_speed.c randombytes.c
	$(CC) $(CFLAGS) -DTEMPO_MATRIX_ALG=2 -DKYBER_K=2 $(SOURCESKECCAK) randombytes.c test/cpucycles.c test/speed_print.c test/test_speed.c -lcrypto -o $@

test/test_speed768_tmp2: $(SOURCESKECCAK) $(HEADERSKECCAK) test/cpucycles.h test/cpucycles.c test/speed_print.h test/speed_print.c test/test_speed.c randombytes.c
	$(CC) $(CFLAGS) -DTEMPO_MATRIX_ALG=2 -DKYBER_K=3 $(SOURCESKECCAK) randombytes.c test/cpucycles.c test/speed_print.c test/test_speed.c -lcrypto -o $@

test/test_speed1024_tmp2: $(SOURCESKECCAK) $(HEADERSKECCAK) test/cpucycles.h test/cpucycles.c test/speed_print.h test/speed_print.c test/test_speed.c randombytes.c
	$(CC) $(CFLAGS) -DTEMPO_MATRIX_ALG=2 -DKYBER_K=4 $(SOURCESKECCAK) randombytes.c test/cpucycles.c test/speed_print.c test/test_speed.c -lcrypto -o $@

# TEMPO ALGORITHM 3a VARIANTS

test/test_kyber512_tmp3a: $(SOURCESKECCAK) $(HEADERSKECCAK) test/test_kyber.c randombytes.c
	$(CC) $(CFLAGS) -DTEMPO_MATRIX_ALG=3 -DKYBER_K=2 $(SOURCESKECCAK) randombytes.c test/test_kyber.c -lcrypto -o $@

test/test_kyber768_tmp3a: $(SOURCESKECCAK) $(HEADERSKECCAK) test/test_kyber.c randombytes.c
	$(CC) $(CFLAGS) -DTEMPO_MATRIX_ALG=3 -DKYBER_K=3 $(SOURCESKECCAK) randombytes.c test/test_kyber.c -lcrypto -o $@

test/test_kyber1024_tmp3a: $(SOURCESKECCAK) $(HEADERSKECCAK) test/test_kyber.c randombytes.c
	$(CC) $(CFLAGS) -DTEMPO_MATRIX_ALG=3 -DKYBER_K=4 $(SOURCESKECCAK) randombytes.c test/test_kyber.c -lcrypto -o $@

test/test_speed512_tmp3a: $(SOURCESKECCAK) $(HEADERSKECCAK) test/cpucycles.h test/cpucycles.c test/speed_print.h test/speed_print.c test/test_speed.c randombytes.c
	$(CC) $(CFLAGS) -DTEMPO_MATRIX_ALG=3 -DKYBER_K=2 $(SOURCESKECCAK) randombytes.c test/cpucycles.c test/speed_print.c test/test_speed.c -lcrypto -o $@

test/test_speed768_tmp3a: $(SOURCESKECCAK) $(HEADERSKECCAK) test/cpucycles.h test/cpucycles.c test/speed_print.h test/speed_print.c test/test_speed.c randombytes.c
	$(CC) $(CFLAGS) -DTEMPO_MATRIX_ALG=3 -DKYBER_K=3 $(SOURCESKECCAK) randombytes.c test/cpucycles.c test/speed_print.c test/test_speed.c -lcrypto -o $@

test/test_speed1024_tmp3a: $(SOURCESKECCAK) $(HEADERSKECCAK) test/cpucycles.h test/cpucycles.c test/speed_print.h test/speed_print.c test/test_speed.c randombytes.c
	$(CC) $(CFLAGS) -DTEMPO_MATRIX_ALG=3 -DKYBER_K=4 $(SOURCESKECCAK) randombytes.c test/cpucycles.c test/speed_print.c test/test_speed.c -lcrypto -o $@

nistkat/PQCgenKAT_kem512_tmp3a: $(SOURCESKECCAK) $(HEADERSKECCAK) nistkat/PQCgenKAT_kem.c nistkat/rng.c nistkat/rng.h
	$(CC) $(NISTFLAGS) -DTEMPO_MATRIX_ALG=3 -DKYBER_K=2 -o $@ $(SOURCESKECCAK) nistkat/rng.c nistkat/PQCgenKAT_kem.c $(LDFLAGS) -lcrypto

nistkat/PQCgenKAT_kem768_tmp3a: $(SOURCESKECCAK) $(HEADERSKECCAK) nistkat/PQCgenKAT_kem.c nistkat/rng.c nistkat/rng.h
	$(CC) $(NISTFLAGS) -DTEMPO_MATRIX_ALG=3 -DKYBER_K=3 -o $@ $(SOURCESKECCAK) nistkat/rng.c nistkat/PQCgenKAT_kem.c $(LDFLAGS) -lcrypto

nistkat/PQCgenKAT_kem1024_tmp3a: $(SOURCESKECCAK) $(HEADERSKECCAK) nistkat/PQCgenKAT_kem.c nistkat/rng.c nistkat/rng.h
	$(CC) $(NISTFLAGS) -DTEMPO_MATRIX_ALG=3 -DKYBER_K=4 -o $@ $(SOURCESKECCAK) nistkat/rng.c nistkat/PQCgenKAT_kem.c $(LDFLAGS) -lcrypto

# TEMPO ALGORITHM 3b VARIANTS

test/test_kyber512_tmp3b: $(SOURCESKECCAK) $(HEADERSKECCAK) test/test_kyber.c randombytes.c
	$(CC) $(CFLAGS) -DTEMPO_MATRIX_ALG=4 -DKYBER_K=2 $(SOURCESKECCAK) randombytes.c test/test_kyber.c -lcrypto -o $@

test/test_kyber768_tmp3b: $(SOURCESKECCAK) $(HEADERSKECCAK) test/test_kyber.c randombytes.c
	$(CC) $(CFLAGS) -DTEMPO_MATRIX_ALG=4 -DKYBER_K=3 $(SOURCESKECCAK) randombytes.c test/test_kyber.c -lcrypto -o $@

test/test_kyber1024_tmp3b: $(SOURCESKECCAK) $(HEADERSKECCAK) test/test_kyber.c randombytes.c
	$(CC) $(CFLAGS) -DTEMPO_MATRIX_ALG=4 -DKYBER_K=4 $(SOURCESKECCAK) randombytes.c test/test_kyber.c -lcrypto -o $@

test/test_speed512_tmp3b: $(SOURCESKECCAK) $(HEADERSKECCAK) test/cpucycles.h test/cpucycles.c test/speed_print.h test/speed_print.c test/test_speed.c randombytes.c
	$(CC) $(CFLAGS) -DTEMPO_MATRIX_ALG=4 -DKYBER_K=2 $(SOURCESKECCAK) randombytes.c test/cpucycles.c test/speed_print.c test/test_speed.c -lcrypto -o $@

test/test_speed768_tmp3b: $(SOURCESKECCAK) $(HEADERSKECCAK) test/cpucycles.h test/cpucycles.c test/speed_print.h test/speed_print.c test/test_speed.c randombytes.c
	$(CC) $(CFLAGS) -DTEMPO_MATRIX_ALG=4 -DKYBER_K=3 $(SOURCESKECCAK) randombytes.c test/cpucycles.c test/speed_print.c test/test_speed.c -lcrypto -o $@

test/test_speed1024_tmp3b: $(SOURCESKECCAK) $(HEADERSKECCAK) test/cpucycles.h test/cpucycles.c test/speed_print.h test/speed_print.c test/test_speed.c randombytes.c
	$(CC) $(CFLAGS) -DTEMPO_MATRIX_ALG=4 -DKYBER_K=4 $(SOURCESKECCAK) randombytes.c test/cpucycles.c test/speed_print.c test/test_speed.c -lcrypto -o $@

nistkat/PQCgenKAT_kem512_tmp3b: $(SOURCESKECCAK) $(HEADERSKECCAK) nistkat/PQCgenKAT_kem.c nistkat/rng.c nistkat/rng.h
	$(CC) $(NISTFLAGS) -DTEMPO_MATRIX_ALG=4 -DKYBER_K=2 -o $@ $(SOURCESKECCAK) nistkat/rng.c nistkat/PQCgenKAT_kem.c $(LDFLAGS) -lcrypto

nistkat/PQCgenKAT_kem768_tmp3b: $(SOURCESKECCAK) $(HEADERSKECCAK) nistkat/PQCgenKAT_kem.c nistkat/rng.c nistkat/rng.h
	$(CC) $(NISTFLAGS) -DTEMPO_MATRIX_ALG=4 -DKYBER_K=3 -o $@ $(SOURCESKECCAK) nistkat/rng.c nistkat/PQCgenKAT_kem.c $(LDFLAGS) -lcrypto

nistkat/PQCgenKAT_kem1024_tmp3b: $(SOURCESKECCAK) $(HEADERSKECCAK) nistkat/PQCgenKAT_kem.c nistkat/rng.c nistkat/rng.h
	$(CC) $(NISTFLAGS) -DTEMPO_MATRIX_ALG=4 -DKYBER_K=4 -o $@ $(SOURCESKECCAK) nistkat/rng.c nistkat/PQCgenKAT_kem.c $(LDFLAGS) -lcrypto

clean:
	-$(RM) -f *.gcno *.gcda *.lcov *.o *.so
	-$(RM) -f test/test_kyber512
	-$(RM) -f test/test_kyber768
	-$(RM) -f test/test_kyber1024
	-$(RM) -f test/test_vectors512
	-$(RM) -f test/test_vectors768
	-$(RM) -f test/test_vectors1024
	-$(RM) -f test/test_kyber512_tmp1
	-$(RM) -f test/test_kyber768_tmp1
	-$(RM) -f test/test_kyber1024_tmp1
	-$(RM) -f test/test_vectors512_tmp1
	-$(RM) -f test/test_vectors768_tmp1
	-$(RM) -f test/test_vectors1024_tmp1
	-$(RM) -f test/test_speed512_tmp1
	-$(RM) -f test/test_speed768_tmp1
	-$(RM) -f test/test_speed1024_tmp1
	-$(RM) -f nistkat/PQCgenKAT_kem512_tmp1
	-$(RM) -f nistkat/PQCgenKAT_kem768_tmp1
	-$(RM) -f nistkat/PQCgenKAT_kem1024_tmp1
	-$(RM) -f test/test_kyber512_tmp2
	-$(RM) -f test/test_kyber768_tmp2
	-$(RM) -f test/test_kyber1024_tmp2
	-$(RM) -f test/test_speed512_tmp2
	-$(RM) -f test/test_speed768_tmp2
	-$(RM) -f test/test_speed1024_tmp2
	-$(RM) -f test/test_kyber512_tmp3a
	-$(RM) -f test/test_kyber768_tmp3a
	-$(RM) -f test/test_kyber1024_tmp3a
	-$(RM) -f test/test_speed512_tmp3a
	-$(RM) -f test/test_speed768_tmp3a
	-$(RM) -f test/test_speed1024_tmp3a
	-$(RM) -f nistkat/PQCgenKAT_kem512_tmp3a
	-$(RM) -f nistkat/PQCgenKAT_kem768_tmp3a
	-$(RM) -f nistkat/PQCgenKAT_kem1024_tmp3a
	-$(RM) -f test/test_kyber512_tmp3a
	-$(RM) -f test/test_kyber768_tmp3a
	-$(RM) -f test/test_kyber1024_tmp3a
	-$(RM) -f test/test_speed512_tmp3b
	-$(RM) -f test/test_speed768_tmp3b
	-$(RM) -f test/test_speed1024_tmp3b
	-$(RM) -f nistkat/PQCgenKAT_kem512_tmp3b
	-$(RM) -f nistkat/PQCgenKAT_kem768_tmp3b
	-$(RM) -f nistkat/PQCgenKAT_kem1024_tmp3b
	-$(RM) -f nistkat/*.req
	-$(RM) -f nistkat/*.rsp
	-$(RM) -rf lib/

