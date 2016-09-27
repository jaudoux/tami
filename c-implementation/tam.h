#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

#define BYTE_TO_BINARY_PATTERN "%c%c%c%c%c%c%c%c"
#define BYTE_TO_BINARY(byte)  \
  (byte & 0x01 ? '1' : '0'), \
  (byte & 0x02 ? '1' : '0'), \
  (byte & 0x04 ? '1' : '0'), \
  (byte & 0x08 ? '1' : '0'), \
  (byte & 0x10 ? '1' : '0'), \
  (byte & 0x20 ? '1' : '0'), \
  (byte & 0x40 ? '1' : '0'), \
  (byte & 0x80 ? '1' : '0')

typedef struct {
  uint32_t ref_id;
  uint32_t pos, length;
  char *seq; // nb_elements => length
  unsigned char *ref_kmers_bv; // nb_elements => (length - k + 1) / 8 + 1
} tam_target_t;

typedef struct {
  uint16_t k;
  uint64_t n_kmers;
  uint32_t n_ref;
  uint32_t n_target;
  char **ref;
  tam_target_t **target;
} tam_header_t;

typedef struct {
  uint32_t target_id;
  uint32_t pos;
  uint8_t ref_seq_l, alt_seq_l;
  uint16_t n_alt_kmers;
  char *alt_seq;
  uint64_t *alt_kmers;
} tam_record_t;

gzFile tam_open(const char *file, const char *mode);

tam_header_t *tam_header_init();
void tam_header_destroy(tam_header_t *h);
int tam_header_write(tam_header_t *h, gzFile fp);
int tam_header_read(tam_header_t *h, gzFile fp);

tam_target_t *tam_target_init();
void tam_target_destroy(tam_target_t *t);
int tam_target_write(tam_target_t *t, gzFile fp);
int tam_target_read(tam_target_t *t, gzFile fp);
void tam_target_print(FILE *stream, const tam_header_t *h, const tam_target_t *t);

tam_record_t *tam_record_init();
void tam_record_destroy(tam_record_t *r);
int tam_record_write(tam_record_t *r, gzFile fp);
int tam_record_read(tam_record_t *r, gzFile fp);
void tam_record_print(FILE *stream, const tam_header_t *h, const tam_record_t *r);
