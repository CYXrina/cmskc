#include "fmm_sketch.h"

void dump_hllmh_vector(hllmh *vec, size_t vlen, size_t init_step, double width_ratio, FILE *outstream)
{
    fwrite(&init_step, sizeof(init_step), 1, outstream);
    fwrite(&width_ratio, sizeof(width_ratio), 1, outstream);
    fwrite(&vlen, sizeof(vlen), 1, outstream);
    for(size_t i = 0; i < vlen; ++i)
    {
        hllmh_serialize_to(&vec[i], outstream);
    }
}
