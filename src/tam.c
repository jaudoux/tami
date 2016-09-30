#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <zlib.h>

#include "dna.h"
#include "tam.h"

gzFile tam_open(const char *path, const char *mode) {
  /* try gzopen */
  return gzopen(path,mode);
}

tam_header_t *tam_header_init() {
  return (tam_header_t*)calloc(1, sizeof(tam_header_t));
}

void tam_header_destroy(tam_header_t *h) {
  if(h == NULL) return;
  if(h->ref) {
    for(int i = 0; i < h->n_ref; i++)
      free(h->ref[i]);
    free(h->ref);
  }
  if(h->target) {
    for(int i = 0; i < h->n_target; i++)
      tam_target_destroy(h->target[i]);
    free(h->target);
  }
  free(h);
}

int tam_header_write(tam_header_t *h, gzFile fp) {
  char magic_number[4] = "TaMF";
  gzwrite(fp, &magic_number, 4);
  gzwrite(fp, &h->k,         sizeof(h->k));
  gzwrite(fp, &h->n_kmers,   sizeof(h->n_kmers));
  gzwrite(fp, &h->n_ref,     sizeof(h->n_ref));
  gzwrite(fp, &h->n_target,  sizeof(h->n_target));
  for(int i = 0; i < h->n_ref; i++) {
    size_t chr_l = strlen(h->ref[i]);
    gzwrite(fp, &chr_l, sizeof(size_t));
    gzwrite(fp, h->ref[i], chr_l);
  }
  for(int i = 0; i < h->n_target; i++) {
    tam_target_write(h->target[i],fp);
  }
  return 1;
}

int tam_header_read(tam_header_t *h, gzFile fp) {
  char magic_number[4];
  size_t r;
  r = gzread(fp, &magic_number, 4);
  if(r == 4) {
    if(strncmp(magic_number,"TaMF", 4) != 0) {
      fprintf(stderr, "File is not in TAM file\n");
      exit(2);
    }
  } else {
    return 0;
  }
  r = gzread(fp, &h->k, sizeof(h->k));
  if(r == sizeof(h->k))       r = gzread(fp, &h->n_kmers,   sizeof(h->n_kmers));  else return 0;
  if(r == sizeof(h->n_kmers)) r = gzread(fp, &h->n_ref,     sizeof(h->n_ref));    else return 0;
  if(r == sizeof(h->n_ref))   r = gzread(fp, &h->n_target,  sizeof(h->n_target)); else return 0;
  if(r == sizeof(h->n_target)) {
    h->ref    = (char **)malloc(sizeof(char**) * h->n_ref);
    //r = fread(kmut->chr,     sizeof(char),   chr_l, fp);
    for(int i = 0; i < h->n_ref; i++) {
      size_t chr_l;
      r = gzread(fp, &chr_l, sizeof(size_t));
      if(r == sizeof(size_t)) {
        h->ref[i] = malloc(sizeof(char) * (chr_l + 1));
        h->ref[i][chr_l] = '\0';
        r = gzread(fp, h->ref[i], chr_l);
        if(r != chr_l) return 0;
      } else {
        return 0;
      }
    }
    h->target = malloc(sizeof(tam_target_t*) * h->n_target);
    for(int i = 0; i < h->n_target; i++) {
      h->target[i] = tam_target_init();
      if(!tam_target_read(h->target[i],fp)) {
        return 0;
      }
    }
  }
  return r;
}

tam_target_t *tam_target_init() {
  return (tam_target_t*)calloc(1, sizeof(tam_target_t));
}

void tam_target_destroy(tam_target_t *t) {
  if(t->seq)           free(t->seq);
  if(t->ref_kmers_bv)  free(t->ref_kmers_bv);
}

int tam_target_write(tam_target_t *t, gzFile fp) {
  // uint32_t ref_id;
  // uint32_t pos, length;
  // char *seq; // nb_elements => length
  // unsigned char *ref_kmers_bv; // nb_elements => (length - k + 1) / 8 + 1
  uint32_t ref_kmers_bv_l = t->length / 8 + 1;
  gzwrite(fp, &t->ref_id, sizeof(t->ref_id));
  gzwrite(fp, &t->pos,    sizeof(t->pos));
  gzwrite(fp, &t->length, sizeof(t->length));
  gzwrite(fp, t->seq,    t->length);
  for(int i = 0; i < ref_kmers_bv_l; i++) {
    gzwrite(fp, &t->ref_kmers_bv[i], sizeof(unsigned char));
  }
  return 1;
}

int tam_target_read(tam_target_t *t, gzFile fp) {
  size_t r;
  r = gzread(fp, &t->ref_id, sizeof(t->ref_id));
  if(r == sizeof(t->ref_id)) r = gzread(fp, &t->pos,    sizeof(t->pos));    else return 0;
  if(r == sizeof(t->pos))    r = gzread(fp, &t->length, sizeof(t->length)); else return 0;
  if(t->seq) t->seq = realloc(t->seq, t->length + 1);
  else       t->seq = malloc(t->length + 1);
  t->seq[t->length] = '\0';
  r = gzread(fp, t->seq, t->length);
  if(r == t->length) {
    uint32_t ref_kmers_bv_l = t->length / 8 + 1;
    if(t->ref_kmers_bv) t->ref_kmers_bv = realloc(t->ref_kmers_bv, ref_kmers_bv_l);
    else                t->ref_kmers_bv = malloc(ref_kmers_bv_l);
    r = gzread(fp, t->ref_kmers_bv, ref_kmers_bv_l);
    if(r != ref_kmers_bv_l) return 0;
  }
  return r;
}

void tam_target_print(FILE *stream, const tam_header_t *h, const tam_target_t *t) {
  fprintf(stream, ">%s:%d-%d\n%s\n", h->ref[t->ref_id], t->pos, t->pos + t->length - 1, t->seq);
  //uint32_t ref_kmers_bv_l = t->length / 8 + 1;
  // for(int i = 0; i < ref_kmers_bv_l; i++) {
  //   fprintf(stream, BYTE_TO_BINARY_PATTERN" ", BYTE_TO_BINARY(t->ref_kmers_bv[i]));
  // }
  // fprintf(stream, "\n");
}

tam_record_t *tam_record_init() {
  return (tam_record_t*)calloc(1, sizeof(tam_record_t));
}

void tam_record_destroy(tam_record_t *r) {
  if(r == NULL) return;
  if(r->alt_seq) free(r->alt_seq);
  if(r->alt_kmers) free(r->alt_kmers);
  free(r);
}

int tam_record_write(tam_record_t *r, gzFile fp) {
  gzwrite(fp, &r->target_id,   sizeof(r->target_id));
  gzwrite(fp, &r->pos,         sizeof(r->pos));
  gzwrite(fp, &r->ref_seq_l,   sizeof(r->ref_seq_l));
  gzwrite(fp, &r->alt_seq_l,   sizeof(r->alt_seq_l));
  gzwrite(fp, &r->n_alt_kmers, sizeof(r->n_alt_kmers));
  gzwrite(fp, r->alt_seq,      r->alt_seq_l);
  gzwrite(fp, r->alt_kmers,    sizeof(uint64_t) * r->n_alt_kmers);
  return 1;
}

int tam_record_read(tam_record_t *d, gzFile fp) {
  size_t r;
  //fwrite(&magic_number, sizeof(char), 4, fp);
  r = gzread(fp, &d->target_id, sizeof(d->target_id));
  if(r == sizeof(d->target_id)) r = gzread(fp, &d->pos,        sizeof(d->pos)); else return 0;
  if(r == sizeof(d->pos))       r = gzread(fp, &d->ref_seq_l,  sizeof(d->ref_seq_l)); else return 0;
  if(r == sizeof(d->ref_seq_l)) r = gzread(fp, &d->alt_seq_l,  sizeof(d->alt_seq_l)); else return 0;
  if(r == sizeof(d->alt_seq_l)) r = gzread(fp, &d->n_alt_kmers, sizeof(d->n_alt_kmers)); else return 0;
  if(r == sizeof(d->n_alt_kmers)) {
    // Alocate space for sequence
    if(d->alt_seq)
      d->alt_seq = realloc(d->alt_seq, d->alt_seq_l + 1);
    else
      d->alt_seq = malloc(d->alt_seq_l + 1);
    d->alt_seq[d->alt_seq_l] = '\0';
    r = gzread(fp, d->alt_seq, d->alt_seq_l);
    if(r == d->alt_seq_l) {
      // Alocate space for k-mers
      // FIXME realloc only if needed (check array size)
      if(d->alt_kmers)
        d->alt_kmers = realloc(d->alt_kmers, sizeof(uint64_t) * d->n_alt_kmers);
      else
        d->alt_kmers = malloc(sizeof(uint64_t) * d->n_alt_kmers);
      r = gzread(fp, d->alt_kmers, sizeof(uint64_t) * d->n_alt_kmers);
      if(r != sizeof(uint64_t) * d->n_alt_kmers) return 0;
    } else {
      return 0;
    }
  } else {
    return 0;
  }
  return r;
}

void tam_record_print(FILE *stream, const tam_header_t *h, const tam_record_t *r) {
  char kmer[h->k + 1];
  kmer[h->k] = '\0';
  char ref_seq[r->ref_seq_l + 1];
  strncpy(ref_seq, &h->target[r->target_id]->seq[r->pos], r->ref_seq_l);
  ref_seq[r->ref_seq_l] = '\0';
  fprintf(stream, "%s\t%d\t%s\t%s", h->ref[h->target[r->target_id]->ref_id], h->target[r->target_id]->pos + r->pos + 1, ref_seq, r->alt_seq);
  // for(int i = 0; i < r->n_ref_kmers; i++) {
  //   int_to_dna(r->ref_kmers[i],h->k,kmer);
  //   if(i > 0) fprintf(stream, ",%s", kmer);
  //   else fprintf(stream, "\t%s", kmer);
  // }
  for(int i = 0; i < r->n_alt_kmers; i++) {
    int_to_dna(r->alt_kmers[i],h->k,kmer);
    if(i > 0) fprintf(stream, ",%s", kmer);
    else fprintf(stream, "\t%s", kmer);
  }
  fprintf(stream, "\n");
}
