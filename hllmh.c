#include "hllmh.h"

#define MAX_B_LEN 10
#define MAX_NMIN 16

#define LEFTMOST_MASK = 0x80 00 00 00 00 00 00 00
const uint8_t MAGIC_NUMBER = 0xAA;

//return number > 0 if cmax is still the maximum
//return 0 if query is equal to cmax
//return number < 0 if cmax must be replaced
int is_max(uint16_t cmax, uint16_t q)
{
    int b_diff = (q >> 10) - (cmax >> 10);
    if(b_diff != 0) return b_diff;
    else return (cmax & 0x03FF) - (q & 0x03FF);
}

int hllmh_init(hllmh *sk, uint8_t prefix_bucket_len, uint8_t b, uint8_t number_of_mins)
{
    sk->pref_len = prefix_bucket_len;
    sk->nbuckets = pow(2, sk->pref_len);
    if(b > MAX_B_LEN) {
        fprintf(stderr, "The maximum dimension for the b-bit part is 10 bits\n");
        return -1;
    }
    sk->b = b;
    if(number_of_mins > MAX_NMIN) {
        fprintf(stderr, "Saving more than 16 minimums is a bit excessive, don't you think?\n");
        return -2;
    }
    sk->nmins = number_of_mins;
    sk->nbins = sk->nbuckets * sk->nmins;
    sk->bins = malloc(sk->nbins, 2);
    if(sk->bins == NULL) return 0;
    for(size_t i = 0; i < sk->nbins; ++i) sk->bins[i] = 0x03FF; //everything is at max value
    return 1;
}

int hllmh_serialize_to(hllmh *sk, FILE *outstream)
{
    if(outstream != NULL)
    {
        fwrite(&MAGIC_NUMBER, sizeof(uint8_t), 1, outstream);
        fwrite(&sk->pref_len, sizeof(uint8_t), 1, outstream);
        fwrite(&sk->b, sizeof(uint8_t), 1, outstream);
        fwrite(&sk->nmins, sizeof(uint8_t), 1, outstream);
        fwrite(sk->bins, sizeof(uint16_t), sk->nbins, outstream);
        return 1;
    } else return 0;
}

int serialize_payload_to(hllmh *sk, FILE *outstream)
{
    if(outstream != NULL)
    {
        fwrite(sk->bins, sizeof(uint16_t), sk->nbins, outstream);
        return 1;
    } else return 0;
}

int hllmh_get_from(FILE* instream, hllmh *sk)
{
    size_t read;
    uint8_t magic;
    
    read = fread(&magic, sizeof(uint8_t), 1, instream);
    if(magic != MAGIC_NUMBER) {
        fprintf(stderr, "File format not compatible. Use --force-read to read it anyway following this software version format\n");
        return 0;
    }
    if(read != sizeof(uint8_t)) {
        fprintf(stderr, "Unable to read the magic number from the stream\n");
        return -3;
    }
    
    read = fread(&sk->pref_len, sizeof(uint8_t), 1, instream);
    if(read != sizeof(uint8_t)) {
        fprintf(stderr, "Unable to read the number of bits to use as bucket index from the stream\n");
        return -3;
    }
    sk->nbuckets = pow(2, sk->pref_len);
    
    read = fread(&sk->b, sizeof(uint8_t), 1, instream);
    if(read != sizeof(uint8_t)) {
        fprintf(stderr, "Unable to read the length of the b-bit part from the stream\n");
        return -3;
    }
    if(sk->b > MAX_B_LEN) {
        fprintf(stderr, "The maximum dimension for the b-bit part is 10 bits, this is not one of my files\n");
        return -1;
    }
    
    read = fread(&sk->nmins, sizeof(uint8_t), 1, instream);
    if(read != sizeof(uint8_t)) {
        fprintf(stderr, "Unable to read the number of hashes to save for each bucket from the stream\n");
        return -3;
    }
    if(&sk->nmins > MAX_NMIN) {
        fprintf(stderr, "Saving more than 16 minimums is excessive, this is not one of my files\n");
        return -2;
    }
    
    sk->nbins = sk->nbuckets * sk->nmins;
    sk->bins = malloc(sk->nbins, 2);
    read = fread(&sk->bins, sizeof(uint16_t), sk->nbins, instream);
    if(read != sk->nbins) {
        fprintf(stderr, "Unable to read the bins from the stream\n");
        return -3;
    }
}

int hllmh_get_payload_from(FILE* instream, uint8_t prefix_bucket_len, uint8_t b, uint8_t number_of_mins, hllmh *sk)
{
    size_t read;
    sk->pref_len = prefix_bucket_len;
    sk->nbuckets = pow(2, sk->pref_len);
    if(b > MAX_B_LEN) {
        fprintf(stderr, "The maximum dimension for the b-bit part is 10 bits\n");
        return -1;
    }
    sk->b = b;
    if(number_of_mins > MAX_NMIN) {
        fprintf(stderr, "Saving more than 16 minimums is a bit excessive, don't you think?\n");
        return -2;
    }
    sk->nmins = number_of_mins;
    sk->nbins = sk->nbuckets * sk->nmins;
    sk->bins = malloc(sk->nbins, 2);
    read = fread(&sk->bins, sizeof(uint16_t), sk->nbins, instream);
    if(read != sk->nbins) {
        fprintf(stderr, "Unable to read the bins from the stream\n");
        return -3;
    }
}

int hllmh_add(hllmh *sk, uint64_t hash)
{
    uint8_t shift = 64 - pref_len;
    uint64_t bucket = hash >> shift;
    uint8_t leading_zeros;
    for(leading_zeros = 0; leading_zeros < shift && !(hash & LEFTMOST_MASK); leading_zeros) hash <<= 1;
    uint16_t h_elem = (hash >> 54) | (leading_zeros << 10); 
    
    uint64_t max_idx;
    uint16_t maximum = 0xFC00; //remember that each hash is saved as [6|10] bits where the maximum for the header is 0 while the rest of the bits follw the general arithmetic rule (remember that the header stores the number of heading zeros).
    const uint64_t bucket_start = bucket * sk->nmins;
    const uint64_t bucket_end = (bucket + 1) * sk->nmins;
    for(size_t i = bucket_start; i < bucket_end; ++i) if(is_max(shift, maximum, sk->bin[i]) < 0)
    {
        max_idx = i;
        maximum = sk->bins[i];
    }
    if(is_max(maximum, h_elem) > 0) sk->bins[max_idx] = h_elem;
    return 1;
}

int hllmh_estimate(hllmh *sk, long *estimate) /*not implemented yet, see HLL++ for details in the future*/
{
    return 1;
}

int hllmh_destroy(hllmh *sk)
{
    free(sk->bins);
    sk->bins = NULL;
    sk->pref_len = sk->b = sk->nmins = sk->nbuckets = sk->nbins = 0;
}
