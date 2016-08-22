#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

typedef struct {
  uint32_t k;
  uint64_t n_kmers;
  uint32_t n_ref;
  char **ref;
} tam_header_t;

typedef struct {
  int32_t ref_id;
  uint32_t pos, n_ref_kmers, n_alt_kmers;
  uint8_t ref_seq_l, alt_seq_l;
  char *ref_seq, *alt_seq;
  uint64_t *ref_kmers, *alt_kmers;
} tam_record_t;

gzFile tam_open(const char *file, const char *mode);

tam_header_t *tam_header_init();
void tam_header_destroy(tam_header_t *h);
int tam_header_write(tam_header_t *h, gzFile fp);
int tam_header_read(tam_header_t *h, gzFile fp);

tam_record_t *tam_record_init();
void tam_record_destroy(tam_record_t *r);
int tam_record_write(tam_record_t *r, gzFile fp);
int tam_record_read(tam_record_t *r, gzFile fp);
