#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <zlib.h>

#include "tam.h"

gzFile *tam_open(const char *path, const char *mode) {
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
  free(h);
}

int tam_header_write(tam_header_t *h, gzFile *fp) {
  char magic_number[4] = "TaMF";
  gzwrite(fp, &magic_number, 4);
  gzwrite(fp, &h->k,         sizeof(h->k));
  gzwrite(fp, &h->n_kmers,   sizeof(h->n_kmers));
  gzwrite(fp, &h->n_ref,     sizeof(h->n_ref));
  for(int i = 0; i < h->n_ref; i++) {
    size_t chr_l = strlen(h->ref[i]);
    gzwrite(fp, &chr_l, sizeof(size_t));
    gzwrite(fp, h->ref[i], chr_l);
  }
  return 1;
}

int tam_header_read(tam_header_t *h, gzFile *fp) {
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
  if(r == sizeof(h->k)) r = gzread(fp, &h->n_kmers, sizeof(h->n_kmers)); else return 0;
  if(r == sizeof(h->n_kmers)) r = gzread(fp, &h->n_ref,   sizeof(h->n_ref)); else return 0;
  if(r == sizeof(h->n_ref)) {
    h->ref = (char **)malloc(sizeof(char**) * h->n_ref);
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
  }
  return r;
}

tam_record_t *tam_record_init() {
  return (tam_record_t*)calloc(1, sizeof(tam_record_t));
}

void tam_record_destroy(tam_record_t *r) {
  if(r == NULL) return;
  if(r->ref_seq) free(r->ref_seq);
  if(r->alt_seq) free(r->alt_seq);
  if(r->ref_kmers) free(r->ref_kmers);
  if(r->alt_kmers) free(r->alt_kmers);
  free(r);
}

int tam_record_write(tam_record_t *r, gzFile *fp) {
  //fwrite(&magic_number, sizeof(char), 4, fp);
  gzwrite(fp, &r->ref_id,      sizeof(r->ref_id));
  gzwrite(fp, &r->pos,         sizeof(r->pos));
  gzwrite(fp, &r->ref_seq_l,   sizeof(r->ref_seq_l));
  gzwrite(fp, &r->alt_seq_l,   sizeof(r->alt_seq_l));
  gzwrite(fp, r->ref_seq,      r->ref_seq_l);
  gzwrite(fp, r->alt_seq,      r->alt_seq_l);
  gzwrite(fp, &r->n_ref_kmers, sizeof(r->n_ref_kmers));
  gzwrite(fp, &r->n_alt_kmers, sizeof(r->n_alt_kmers));
  gzwrite(fp, r->ref_kmers,    sizeof(uint64_t) * r->n_ref_kmers);
  gzwrite(fp, r->alt_kmers,    sizeof(uint64_t) * r->n_alt_kmers);
  return 1;
}

int tam_record_read(tam_record_t *d, gzFile *fp) {
  size_t r;
  //fwrite(&magic_number, sizeof(char), 4, fp);
  r = gzread(fp, &d->ref_id,      sizeof(d->ref_id));
  if(r == sizeof(d->ref_id))    r = gzread(fp, &d->pos,        sizeof(d->pos)); else return 0;
  if(r == sizeof(d->pos))       r = gzread(fp, &d->ref_seq_l,  sizeof(d->ref_seq_l)); else return 0;
  if(r == sizeof(d->ref_seq_l)) r = gzread(fp, &d->alt_seq_l,  sizeof(d->alt_seq_l)); else return 0;
  if(r == sizeof(d->alt_seq_l)) {
    // Alocate space for sequence
    if(d->ref_seq)
      d->ref_seq = realloc(d->ref_seq, d->ref_seq_l + 1);
    else
      d->ref_seq = malloc(d->ref_seq_l + 1);
    if(d->alt_seq)
      d->alt_seq = realloc(d->alt_seq, d->alt_seq_l + 1);
    else
      d->alt_seq = malloc(d->alt_seq_l + 1);
    d->ref_seq[d->ref_seq_l] = '\0';
    d->alt_seq[d->alt_seq_l] = '\0';
    r = gzread(fp, d->ref_seq, d->ref_seq_l);
    if(r == d->ref_seq_l) r = gzread(fp, d->alt_seq, d->alt_seq_l); else return 0;
    if(r == d->alt_seq_l) r = gzread(fp, &d->n_ref_kmers, sizeof(d->n_ref_kmers)); else return 0;
    if(r == sizeof(d->n_ref_kmers)) r = gzread(fp, &d->n_alt_kmers, sizeof(d->n_alt_kmers)); else return 0;
    if(r == sizeof(d->n_alt_kmers)) {
      // Alocate space for k-mers
      // FIXME realloc only if needed (check array size)
      if(d->ref_kmers)
        d->ref_kmers = realloc(d->ref_kmers, sizeof(uint64_t) * d->n_ref_kmers);
      else
        d->ref_kmers = malloc(sizeof(uint64_t) * d->n_ref_kmers);
      if(d->alt_kmers)
        d->alt_kmers = realloc(d->alt_kmers, sizeof(uint64_t) * d->n_alt_kmers);
      else
        d->alt_kmers = malloc(sizeof(uint64_t) * d->n_alt_kmers);
      r = gzread(fp, d->ref_kmers, sizeof(uint64_t) * d->n_ref_kmers);
      if(r == sizeof(uint64_t) * d->n_ref_kmers)
        r = gzread(fp, d->alt_kmers, sizeof(uint64_t) * d->n_alt_kmers); else return 0;
      if(r != sizeof(uint64_t) * d->n_alt_kmers) return 0;
    } else {
      return 0;
    }
  } else {
    return 0;
  }
  return r;
}
