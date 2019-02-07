#include <inttypes.h>
#include <stdlib.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
Implementation of HyperMinHash.
Each bin is a 16 bit unsigned integer, 
6 bits for the length of the trailing zeros,
10 bits for the b-bit part with b able to vary between 0 and 10.
*/

typedef struct {
    uint8_t pref_len;
    uint8_t b;
    uint8_t nmins;
    uint64_t nbuckets;
    uint64_t nbins;
    uint16_t *bins;
} hllmh;

int hllmh_init(hllmh *sk, uint8_t prefix_bucket_len, uint8_t b, uint8_t number_of_mins);
int hllmh_serialize_to(hllmh *sk, FILE *outstream);
int hllmh_get_from(FILE* instream, hllmh *sk);
int hllmh_add(hllmh *sk, uint64_t hash);
int hllmh_estimate(hllmh *sk, long *estimate); /*not implemented yet, see HLL++ for details in the future*/
int hllmh_destroy(hllmh *sk);

#ifdef __cplusplus
}
#endif
