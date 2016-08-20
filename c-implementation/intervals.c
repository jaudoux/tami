#include "intervals.h"

#include <zlib.h>

#include "kstring.h"
#include "kvec.h"
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)

interval_t *interval_init(char *_chr, int _start, int _end) {
  interval_t *interval = (interval_t*)malloc(sizeof(interval_t));
  interval->chr   = _chr;
  interval->start = _start;
  interval->end   = _end;
  interval->seq   = NULL;
  return interval;
}

void interval_destroy(interval_t *interval) {
  if(interval->seq)
    free(interval->seq);
  if(interval->chr)
    free(interval->chr);
  free(interval);
}

int cmp_interval(const void * a, const void * b) {
  const interval_t *int_a = *(const interval_t **)a;
  const interval_t *int_b = *(const interval_t **)b;
  int chr_cmp = strcmp(int_a->chr,int_b->chr);
  if(chr_cmp == 0) {
    return ( int_a->start - int_b->start );
  } else {
    return chr_cmp;
  }
}

int interval_length(interval_t *interval) {
  return interval->end - interval->start + 1;
}

interval_array_t *interval_array_init() {
  interval_array_t *a = (interval_array_t*)calloc(1, sizeof(interval_array_t));
  kv_init(*a);
  return a;
}

void interval_array_destroy(interval_array_t *a) {
  for(int i = 1; i < kv_size(*a); i++)
    interval_destroy(kv_A(*a,i));
  kv_destroy(*a);
}

int load_intervals_from_bed(const char * bed_file, interval_array_t *a) {
  gzFile fp;
	kstream_t *ks;
	kstring_t *str, *chr;
  int start, end, dret;
  interval_t *interval;

  str = calloc(1, sizeof(kstring_t));
  chr = calloc(1, sizeof(kstring_t));
  fp = strcmp(bed_file, "-")? gzopen(bed_file, "r") : gzdopen(fileno(stdin), "r");
  ks = ks_init(fp);

  while (ks_getuntil(ks, 0, str, &dret) >= 0) {
    // remove chr prefix if any
    if(str->l > 3 && strncmp("chr", str->s, 3) == 0)
      kputs(&str->s[3],chr);
    else
      kputs(str->s,chr);

    if(dret != '\n') {
      if(ks_getuntil(ks, 0, str, &dret) > 0 && isdigit(str->s[0])) {
        start = atoi(str->s);
        if(dret != '\n') {
          if(ks_getuntil(ks, 0, str, &dret) > 0 && isdigit(str->s[0])) {
            end = atoi(str->s);
            // create a new interval
            interval = interval_init(ks_release(chr),start + 1, end);
            kv_push(interval_t*,*a,interval);
          }
        }
      }
    }
    // skip the rest of the line
		if (dret != '\n') while ((dret = ks_getc(ks)) > 0 && dret != '\n');
    // release the current chr
    if(chr->l > 0) chr->l = 0;
  }
  ks_destroy(ks);
  gzclose(fp);
  free(str->s); free(str);
  free(chr->s); free(chr);
  return kv_size(*a);
}

int merge_overlapping_intervals(interval_array_t *a) {
  if(kv_size(*a) == 0) return 0;
  int nb_merged = 0, j = 0;
  // Sort and remove overlapping intervals
  qsort(a->a, kv_size(*a), sizeof(a->a[0]), cmp_interval);
  for(int i = 1; i < kv_size(*a); i++) {
    if(strcmp(a->a[j]->chr,a->a[i]->chr) == 0 &&
      a->a[i]->start <= a->a[j]->end) {
      if(a->a[i]->end > a->a[j]->end)
        a->a[j]->end = a->a[i]->end;
      nb_merged++;
      interval_destroy(a->a[i]);
    } else {
      a->a[++j] = a->a[i];
    }
  }
  // Set the new size of the interval array
  a->n = j + 1;
  return nb_merged;
}

int load_intervals_seq_from_fasta(const char* fasta_file, interval_array_t *a) {
  gzFile fp;
  int l, chr_l, int_l, nb_load = 0;
  kseq_t *seq;
  interval_t *interval;
  fp = gzopen(fasta_file, "r");
  seq = kseq_init(fp);
  while ((l = kseq_read(seq)) >= 0) {
    chr_l = strlen(seq->seq.s);
    for(int i = 0; i < kv_size(*a); i++) {
      interval = kv_A(*a,i);
      if(!interval->seq && chr_l >= interval->end &&
        strcmp(interval->chr, seq->name.s) == 0) {
        int_l = interval_length(interval);
        interval->seq = malloc(int_l + 1);
        memcpy(interval->seq, &seq->seq.s[interval->start - 1], int_l);
        interval->seq[int_l] = '\0'; // FIXME is it necessary?
        nb_load++;
      }
    }
  }
  return nb_load;
}
