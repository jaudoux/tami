#include "kvec.h"

typedef struct {
  char *chr;
  int start, end;
  char *seq;
} interval_t;

typedef kvec_t(interval_t*) interval_array_t;

interval_t *interval_init(char *_chr, int _start, int _end);

void interval_destroy(interval_t *interval);

int cmp_interval(const void * a, const void * b);

int interval_length(interval_t *interval);

interval_array_t *interval_array_init();

void interval_array_destroy(interval_array_t *a);

int load_intervals_from_bed(const char * bed_file, interval_array_t *a);

int merge_overlapping_intervals(interval_array_t *a);

int load_intervals_seq_from_fasta(const char* fasta_file, interval_array_t *a);
