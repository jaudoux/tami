#ifndef DNA_H
#define DNA_H

#include <stdint.h>

/*
basemap[] works by storing a very small array that maps a base to
its complement, by dereferencing the array with the ASCII char's
decimal value as the index
(int) 'A' = 65;
(int) 'C' = 67;
(int) 'G' = 71;
(int) 'T' = 84;
(int) 'a' = 97;
(int) 'c' = 99;
(int) 'g' = 103;
(int) 't' = 116;
(int) 'N' = 78;
(int) 'U' = 85;
(int) 'u' = 117;
for example: basemap['A'] => basemap[65] => 'T' etc.
*/

static const int base_to_int[255] =
{
  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, /*   0 -   9 */
  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, /*  10 -  19 */
  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, /*  20 -  29 */
  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, /*  30 -  39 */
  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, /*  40 -  49 */
  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, /*  50 -  59 */
  -1, -1, -1, -1, -1,  0, -1,  1, -1, -1, /*  60 -  69 */
  -1,  2, -1, -1, -1, -1, -1, -1, -1, -1, /*  70 -  79 */
  -1, -1, -1, -1,  3,  3, -1, -1, -1, -1, /*  80 -  89 */
  -1, -1, -1, -1, -1, -1, -1,  0, -1,  1, /*  90 -  99 */
  -1, -1, -1,  2, -1, -1, -1, -1, -1, -1, /* 100 - 109 */
  -1, -1, -1, -1, -1, -1,  3,  3, -1, -1, /* 110 - 119 */
  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, /* 120 - 129 */
  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, /* 130 - 139 */
  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, /* 140 - 149 */
  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, /* 150 - 159 */
  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, /* 160 - 169 */
  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, /* 170 - 179 */
  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, /* 180 - 189 */
  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, /* 190 - 199 */
  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, /* 200 - 209 */
  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, /* 210 - 219 */
  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, /* 220 - 229 */
  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, /* 230 - 239 */
  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, /* 240 - 249 */
  -1, -1, -1, -1, -1                                /* 250 - 254 */
};

static const char basecomp[255] =
    {
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*   0 -   9 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*  10 -  19 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*  20 -  29 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*  30 -  39 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*  40 -  49 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*  50 -  59 */
        '\0', '\0', '\0', '\0', '\0',  'T', '\0',  'G', '\0', '\0', /*  60 -  69 */
        '\0',  'C', '\0', '\0', '\0', '\0', '\0', '\0',  'N', '\0', /*  70 -  79 */
        '\0', '\0', '\0', '\0',  'A',  'A', '\0', '\0', '\0', '\0', /*  80 -  89 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0',  'T', '\0',  'G', /*  90 -  99 */
        '\0', '\0', '\0',  'C', '\0', '\0', '\0', '\0', '\0', '\0', /* 100 - 109 */
        '\0', '\0', '\0', '\0', '\0', '\0',  'A',  'A', '\0', '\0', /* 110 - 119 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 120 - 129 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 130 - 139 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 140 - 149 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 150 - 159 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 160 - 169 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 170 - 179 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 180 - 189 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 190 - 199 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 200 - 209 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 210 - 219 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 220 - 229 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 230 - 239 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 240 - 249 */
        '\0', '\0', '\0', '\0', '\0'                                /* 250 - 254 */
    };

static const char NUCLEOTIDES[4] = { 'A', 'C', 'G', 'T' };

static const int NB_NUCLEOTIDES = 4;

static const uint64_t DNA_BIT_MASK[32] = {
  3UL, 12UL, 48UL, 192UL, 768UL, 3072UL, 12288UL, 49152UL, 196608UL,
  786432UL, 3145728UL, 12582912UL, 50331648UL, 201326592UL, 805306368UL,
  3221225472UL, 12884901888UL, 51539607552UL, 206158430208UL, 824633720832UL,
  3298534883328UL, 13194139533312UL, 52776558133248UL, 211106232532992UL,
  844424930131968UL, 3377699720527872UL, 13510798882111488UL,
  54043195528445952UL, 216172782113783808UL, 864691128455135232UL,
  3458764513820540928UL, 13835058055282163712UL

};

static inline uint64_t int_revcomp(uint64_t factor, uint32_t length) {
  uint64_t mask;
  if (length == 32)
    mask = ~0;
  else
    mask =  ((uint64_t) 1 << (2*length)) - 1;

  factor ^= mask;

  uint64_t mask_lsb;
  // Corresponds to the rightmost nucleotide
  mask_lsb = 3;
  uint64_t shift = 0;
  uint64_t result = 0;
  for(int j = 0; j < length; j++){
    result <<= 2;
    // get the leftmost nucleotide and put it at the end
    result |= (factor & mask_lsb) >> shift;
    mask_lsb <<= 2;
    shift += 2;
  }

  return result;
}

// static inline uint64_t dna_to_int(const char *dna, size_t dna_length, int canonical) {
//   uint64_t dna_int = 0;
//   // TODO we should check that dna_length is smaller that 32.
//   for (size_t i = 0; i< dna_length ; i++) {
//     dna_int <<= 2;
//     dna_int |= base_to_int[(int)dna[i]];
//   }
//   // If the conversion is not "strand-specific" we calculate the reverse DNA it
//   // and return the one that has the smallest value
//   if(canonical) {
//     uint64_t rev_comp = int_revcomp(dna_int,dna_length);
//     if(rev_comp < dna_int) {
//       return rev_comp;
//     }
//   }
//   return dna_int;
// }

static inline void int_to_dna(uint64_t code, size_t dna_length, char *dna) {
  uint64_t mask = 3;
  for (int i = 0; i < dna_length; i++) {
    dna[dna_length-i-1] = NUCLEOTIDES[code & mask];
    code >>=2;
  }
}

static inline uint64_t mut_int_dna(uint64_t code, size_t dna_length, int pos, char n) {
  uint64_t x = (dna_length - pos - 1) * 2;
  uint64_t up = 1;
  if(n == 'A') {
    code &= ~(up << x);
    code &= ~(up << (x+1));
  } else if(n == 'C') {
    code |= up << x;
    code &= ~(up << (x+1));
  } else if(n == 'G') {
    code &= ~(up << x);
    code |= up << (x+1);
  } else {
    code |= up << x;
    code |= up << (x+1);
  }
  return code;
}

static inline char nuc_from_int_dna(uint64_t code, size_t dna_length, int pos) {
  uint64_t nuc_bit_pos = dna_length - pos - 1;
  return NUCLEOTIDES[(code & DNA_BIT_MASK[nuc_bit_pos]) >> (nuc_bit_pos*2)];
}


static inline uint64_t canonical_kmer(const char *dna, size_t dna_length, uint64_t *forward_kmer, uint64_t *reverse_kmer) {
  // TODO we should check that dna_length is smaller that 32.
  *forward_kmer = 0;
  for (size_t i = 0; i< dna_length ; i++) {
    *forward_kmer <<= 2;
    *forward_kmer |= base_to_int[(int)dna[i]];
  }
  *reverse_kmer = int_revcomp(*forward_kmer, dna_length);
  return *forward_kmer < *reverse_kmer? *forward_kmer : *reverse_kmer;
}

static inline uint64_t next_canonical_kmer(size_t dna_length, char new_nuc, uint64_t *forward_kmer, uint64_t *reverse_kmer) {
  uint64_t x = (dna_length - 1) * 2;
  uint64_t up = 3;
  // Unset the first nucleotide
  *forward_kmer &= ~(up << x);
  *forward_kmer <<= 2;
  *forward_kmer = mut_int_dna(*forward_kmer, dna_length, dna_length - 1, new_nuc);
  *reverse_kmer >>= 2;
  *reverse_kmer = mut_int_dna(*reverse_kmer, dna_length, 0, basecomp[(int)new_nuc]);
  return *forward_kmer < *reverse_kmer? *forward_kmer : *reverse_kmer;
}

#endif
