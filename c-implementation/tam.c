#include <string.h>
#include <zlib.h>

#include "tam.h"

FILE *tam_open(const char *path, const char *mode) {
  /* try gzopen */
  gzFile * gzf = gzopen(path,mode);

  /* open file pointer */
  return funopen(gzf,
                 (int(*)(void*,char*,int))gzread,
                 (int(*)(void*,const char*,int))gzwrite,
                 (fpos_t(*)(void*,fpos_t,int))gzseek,
                 (int(*)(void*))gzclose);
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

int tam_header_write(tam_header_t *h, FILE *fp) {
  char magic_number[4] = "TaMF";
  fwrite(&magic_number, sizeof(char), 4, fp);
  fwrite(&h->k,       sizeof(h->k),       1, fp);
  fwrite(&h->n_kmers, sizeof(h->n_kmers), 1, fp);
  fwrite(&h->n_ref,   sizeof(h->n_ref),   1, fp);
  for(int i = 0; i < h->n_ref; i++) {
    size_t chr_l = strlen(h->ref[i]);
    fwrite(&chr_l,      sizeof(size_t),     1, fp);
    fwrite(h->ref[i],  sizeof(char),   chr_l, fp);
  }
  return 1;
}

int tam_header_read(tam_header_t *h, FILE *fp) {
  char magic_number[4];
  size_t r;
  r = fread(&magic_number, sizeof(char), 4, fp);
  if(r == 4) {
    if(strncmp(magic_number,"TaMF", 4) != 0) {
      fprintf(stderr, "File is not in TAM file\n");
      exit(2);
    }
  } else {
    return 0;
  }
  r = fread(&h->k, sizeof(h->k), 1, fp);
  if(r == 1) r = fread(&h->n_kmers, sizeof(h->n_kmers), 1, fp); else return 0;
  if(r == 1) r = fread(&h->n_ref,   sizeof(h->n_ref),   1, fp); else return 0;
  if(r == 1) {
    h->ref = (char **)malloc(sizeof(char**) * h->n_ref);
    //r = fread(kmut->chr,     sizeof(char),   chr_l, fp);
    for(int i = 0; i < h->n_ref; i++) {
      size_t chr_l;
      r = fread(&chr_l, sizeof(size_t), 1, fp);
      if(r == 1) {
        h->ref[i] = malloc(sizeof(char) * (chr_l + 1));
        h->ref[i][chr_l] = '\0';
        r = fread(h->ref[i], sizeof(char), chr_l, fp);
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

int tam_record_write(tam_record_t *r, FILE *fp) {
  //fwrite(&magic_number, sizeof(char), 4, fp);
  fwrite(&r->ref_id,      sizeof(r->ref_id),                 1, fp);
  fwrite(&r->pos,         sizeof(r->pos),                    1, fp);
  fwrite(&r->ref_seq_l,   sizeof(r->ref_seq_l),              1, fp);
  fwrite(&r->alt_seq_l,   sizeof(r->alt_seq_l),              1, fp);
  fwrite(r->ref_seq,     sizeof(char),           r->ref_seq_l, fp);
  fwrite(r->alt_seq,     sizeof(char),           r->alt_seq_l, fp);
  fwrite(&r->n_ref_kmers, sizeof(r->n_ref_kmers),            1, fp);
  fwrite(&r->n_alt_kmers, sizeof(r->n_alt_kmers),            1, fp);
  fwrite(r->ref_kmers,   sizeof(uint64_t),     r->n_ref_kmers, fp);
  fwrite(r->alt_kmers,   sizeof(uint64_t),     r->n_alt_kmers, fp);
  return 1;
}

int tam_record_read(tam_record_t *d, FILE *fp) {
  size_t r;
  //fwrite(&magic_number, sizeof(char), 4, fp);
  r = fread(&d->ref_id,      sizeof(d->ref_id),     1, fp);
  if(r == 1) r = fread(&d->pos,         sizeof(d->pos),        1, fp); else return 0;
  if(r == 1) r = fread(&d->ref_seq_l,   sizeof(d->ref_seq_l),  1, fp); else return 0;
  if(r == 1) r = fread(&d->alt_seq_l,   sizeof(d->alt_seq_l),  1, fp); else return 0;
  if(r == 1) {
    // Alocate space for sequence
    if(d->ref_seq) d->ref_seq = realloc(d->ref_seq, sizeof(char) * (d->ref_seq_l + 1));
    else d->ref_seq = malloc(sizeof(char) * (d->ref_seq_l + 1));
    if(d->alt_seq) d->alt_seq = realloc(d->alt_seq, sizeof(char) * (d->alt_seq_l + 1));
    else d->alt_seq = malloc(sizeof(char) * (d->alt_seq_l + 1));
    d->ref_seq[d->ref_seq_l] = '\0';
    d->alt_seq[d->alt_seq_l] = '\0';
    r = fread(d->ref_seq, sizeof(char), d->ref_seq_l, fp);
    if(r == d->ref_seq_l) r = fread(d->alt_seq, sizeof(char), d->alt_seq_l, fp); else return 0;
    if(r == d->alt_seq_l) r = fread(&d->n_ref_kmers, sizeof(d->n_ref_kmers), 1, fp); else return 0;
    if(r == 1) r = fread(&d->n_alt_kmers, sizeof(d->n_alt_kmers), 1, fp); else return 0;
    if(r == 1) {
      // Alocate space for k-mers
      // FIXME realloc only if needed (check array size)
      if(d->ref_kmers) d->ref_kmers = realloc(d->ref_kmers, sizeof(uint64_t) * d->n_ref_kmers);
      else d->ref_kmers = malloc(sizeof(uint64_t) * d->n_ref_kmers);
      if(d->alt_kmers) d->alt_kmers = realloc(d->alt_kmers, sizeof(uint64_t) * d->n_alt_kmers);
      else d->alt_kmers = malloc(sizeof(uint64_t) * d->n_alt_kmers);
      r = fread(d->ref_kmers, sizeof(uint64_t), d->n_ref_kmers, fp);
      if(r == d->n_ref_kmers) r = fread(d->alt_kmers, sizeof(uint64_t), d->n_alt_kmers, fp); else return 0;
      if(r == d->n_alt_kmers) return 1; else return 0;
    } else {
      return 0;
    }
  } else {
    return 0;
  }

}
