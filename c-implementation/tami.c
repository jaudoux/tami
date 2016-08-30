#include <stdio.h>  //fprintf
#include <stdlib.h> //free
#include <zlib.h>
#include <inttypes.h>
#include <math.h> // pow()

#include "kstring.h"
#include "ksort.h"
#include "khash.h"
#include "kvec.h"
#include "kthread.h"
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)
KHASH_MAP_INIT_INT64(kmers, uint32_t)
KHASH_MAP_INIT_STR(references, int32_t)

typedef khash_t(kmers) kmers_hash_t;

#include "api.h"
#include "dna.h"
#include "tam.h"
#include "intervals.h"

#define TAMI_VERSION "0.2.0"
#define NB_THREAD_API 10
#define DEFAULT_OUTPUT_NAME "kmers.tam"
#define DEFAULT_K_LENGTH 32

#define REFERENCE_KMER_CHECKED 0
#define REFERENCE_KMER 1
#define MUTATED_KMER 2
#define DOUBLE_MUTATED_KMER 3

long long combi(int n,int k) {
  long long ans = 1;
  k = k > n - k? n - k : k;
  int j = 1;
  for(; j <= k; j++, n--) {
    if(n % j == 0) {
      ans *= n / j;
    } else
    if(ans % j == 0) {
      ans = ans / j * n;
    } else {
      ans = (ans * n) / j;
    }
  }
  return ans;
}

unsigned long factorial(unsigned long f) {
  if (f == 0)
    return 1;
  return(f * factorial(f - 1));
}


char *get_tmp_filename(const char * f) {
  char * tmp = malloc(strlen(f) + strlen(".tmp") + 1);
  tmp[0] = '\0';
  strcat(tmp,f);
  strcat(tmp,".tmp");
  return tmp;
}

static void get_sequence_iter(void *_int_a, long i, int tid) {
  interval_t ** int_a   = (interval_t **)_int_a;
  interval_t *interval = int_a[i];
  interval->seq = get_sequence(interval->chr, interval->start, interval->end);
}

int remove_reference_kmers(char *reference_fasta, kmers_hash_t *h, int k_length) {
  int l, nb_removed_kmers = 0;
  gzFile fp;
  khiter_t k;
  kseq_t *seq;
  uint64_t kmer_int, forward_kmer_int = 0, reverse_kmer_int = 0;

  fp  = gzopen(reference_fasta, "r");
  seq = kseq_init(fp);
  while ((l = kseq_read(seq)) >= 0) {
    if(seq->seq.l >= k_length) {
      //canonical_kmer(const char *dna, size_t dna_length, uint64_t &forward_kmer, uint64_t &reverse_kmer)
      for(int i = 0; i < seq->seq.l - k_length + 1; i++) {
        if(i > 0) {
          kmer_int = next_canonical_kmer(k_length, seq->seq.s[i + k_length - 1], &forward_kmer_int, &reverse_kmer_int);
        } else {
          kmer_int = canonical_kmer(seq->seq.s, k_length, &forward_kmer_int, &reverse_kmer_int);
        }
        k = kh_get(kmers, h, kmer_int);
        if(k != kh_end(h)) {

          if(kh_value(h, k) == REFERENCE_KMER) {
            // If a reference k-mer is seen more than once, it will be deleted
            kh_value(h, k) = REFERENCE_KMER_CHECKED;
          } else {
            kh_del(kmers, h, k);
            nb_removed_kmers++;
          }
        }
      }
    }
  }

  // Reset reference k-mer checked
  for(k = kh_begin(h); k != kh_end(h); ++k) {
    if (!kh_exist(h, k)) continue;
    if(kh_value(h, k) == REFERENCE_KMER_CHECKED)
      kh_value(h, k) = REFERENCE_KMER;
  }

  return nb_removed_kmers;
}

int load_kmers(char *tam_path, kmers_hash_t *h, int set_to_zero) {
  int ret;
  khiter_t k;
  gzFile tam_file = tam_open(tam_path, "rb");
  tam_header_t *tam_header = tam_header_init();
  tam_record_t *tam_record = tam_record_init();
  if(tam_header_read(tam_header,tam_file)) {
    kh_resize(kmers, h, kh_size(h) + tam_header->n_kmers);
    while(tam_record_read(tam_record,tam_file)) {
      for(int i = 0; i < tam_record->n_ref_kmers; i++) {
        k = kh_put(kmers, h, tam_record->ref_kmers[i], &ret);
        if(set_to_zero)
          kh_value(h, k) = 0;
        else
          kh_value(h, k) = REFERENCE_KMER;
      }
      for(int i = 0; i < tam_record->n_alt_kmers; i++) {
        k = kh_put(kmers, h, tam_record->alt_kmers[i], &ret);
        if(set_to_zero)
          kh_value(h, k) = 0;
        else
          kh_value(h, k) = MUTATED_KMER;
      }
    }
  }
  tam_header_destroy(tam_header);
  tam_record_destroy(tam_record);

  return (int)kh_size(h);
}

int update_tam_file(const char *tam_path, const char *new_tam_path, kmers_hash_t *h) {
  khiter_t k;
  gzFile tam_file     = tam_open(tam_path, "rb");
  gzFile new_tam_file = tam_open(new_tam_path, "wb");
  tam_header_t *tam_header = tam_header_init();
  tam_record_t *tam_record = tam_record_init();
  int nb_record_removed = 0;

  if(tam_header_read(tam_header,tam_file)) {
    tam_header->n_kmers = kh_size(h);
    tam_header_write(tam_header, new_tam_file);
    while(tam_record_read(tam_record,tam_file)) {
      // Verify ref_kmers
      int j = 0;
      for(int i = 0; i < tam_record->n_ref_kmers; i++) {
        k = kh_get(kmers, h, tam_record->ref_kmers[i]);
        if(k != kh_end(h) && kh_value(h, k) == REFERENCE_KMER)
          tam_record->ref_kmers[j++] = tam_record->ref_kmers[i];
      }
      tam_record->n_ref_kmers = j;

      // Verify alt_kmers
      j = 0;
      for(int i = 0; i < tam_record->n_alt_kmers; i++) {
        k = kh_get(kmers, h, tam_record->alt_kmers[i]);
        if(k != kh_end(h) && kh_value(h, k) != REFERENCE_KMER)
          tam_record->alt_kmers[j++] = tam_record->alt_kmers[i];
      }
      tam_record->n_alt_kmers = j;

      if(tam_record->n_ref_kmers == 0 || tam_record->n_alt_kmers == 0)
        nb_record_removed++;
      else
        tam_record_write(tam_record, new_tam_file);
    }
  }

  gzclose(tam_file);
  gzclose(new_tam_file);
  tam_header_destroy(tam_header);
  tam_record_destroy(tam_record);
  return nb_record_removed;
}

int tami_vcfbuild(int argc, char *argv[]) {
  char *vcf_file = NULL, *reference_fasta = NULL, *output_file = DEFAULT_OUTPUT_NAME;
  int k_length = DEFAULT_K_LENGTH;

  int c;
  while ((c = getopt(argc, argv, "k:r:o:")) >= 0) {
		switch (c) {
      case 'k': k_length = atoi(optarg);  break;
      case 'r': reference_fasta = optarg; break;
      case 'o': output_file = optarg;     break;
		}
	}

  if (optind == argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   tami vcf-build [options] <in.vcf>\n\n");
		fprintf(stderr, "Options: -k INT    length of k-mers (max_value: 32) [%d]\n", k_length);
    fprintf(stderr, "         -r STR    reference FASTA (otherwise sequences are retrieved from rest.ensembl.org)\n");
    fprintf(stderr, "         -o STR    output file (default : 'kmers.tam')\n");
		fprintf(stderr, "\n");
		return 1;
	}

  // verify Options
  if(k_length > 32 || k_length < 1) {
    fprintf(stderr, "Invalid value (%d) for k_length (-k option).\n", k_length);
    return 1;
  }

  vcf_file = argv[optind++];

  int dret, ret, ref_id, pos;
  gzFile fp;
  khiter_t k, k2;
	kstream_t *ks;
  kstring_t *str, *chr, *ref_seq;
  tam_record_t *tam_record;// = tam_record_init();
  tam_header_t *tam_header;
  khash_t(references) *h_r = kh_init(references);
  khash_t(kmers) *h_k = kh_init(kmers);
  kvec_t(char*) reference_names;
  kvec_t(tam_record_t*) tam_records;
  char *kmer = malloc(k_length);
  kseq_t *seq;
  int l, chr_l;

  kv_init(reference_names);
  kv_init(tam_records);
  str = calloc(1, sizeof(kstring_t));
  chr = calloc(1, sizeof(kstring_t));
  ref_seq = calloc(1, sizeof(kstring_t));

  fp = gzopen(vcf_file, "r");
  if(!fp) { fprintf(stderr, "Failed to open %s\n", vcf_file); exit(EXIT_FAILURE); }
  ks = ks_init(fp);

  fprintf(stderr, "Reading VCF file...\n");
  while (ks_getuntil(ks, 0, str, &dret) >= 0) {

    if(str->l > 3 && strncmp("chr", str->s, 3) == 0)
      kputs(&str->s[3],chr);
    else
      kputs(str->s,chr);

    k = kh_get(references, h_r, chr->s);
    if(k == kh_end(h_r)) {
      char *chr_cpy = malloc(strlen(chr->s) + 1);
      strcpy(chr_cpy,chr->s);
      kv_push(char*,reference_names,chr_cpy);
      k = kh_put(references, h_r, chr_cpy, &ret);
      kh_value(h_r, k) = kv_size(reference_names) - 1;
      ref_id = kv_size(reference_names) - 1;
    } else {
      ref_id = kh_value(h_r, k);
    }

    if(dret != '\n' && ks_getuntil(ks, 0, str, &dret) > 0 && isdigit(str->s[0])) {
      pos = atoi(str->s);
      if(dret != '\n' && ks_getuntil(ks, 0, str, &dret) > 0) {
        // TODO now STR holds the variant id
        if(dret != '\n' && ks_getuntil(ks, 0, str, &dret) > 0) {
          if(str->l <= k_length) {
            kputs(str->s,ref_seq);
            //strcpy(tam_record->ref_seq, str->s);
            if(dret != '\n' && ks_getuntil(ks, 0, str, &dret) > 0) {
              /* get the first token */
              char * alt = strtok(str->s, ",");
              /* walk through other tokens */
              while( alt != NULL )
              {
                if(strlen(alt) <= k_length) {
                  tam_record = tam_record_init();
                  tam_record->ref_id = ref_id;
                  tam_record->pos    = pos - 1; // VCF are 1-based
                  tam_record->ref_seq = malloc(ref_seq->l + 1);
                  tam_record->alt_seq = malloc(strlen(alt) + 1);
                  strcpy(tam_record->ref_seq, ref_seq->s);
                  strcpy(tam_record->alt_seq, alt);
                  kv_push(tam_record_t*,tam_records,tam_record);
                }
                alt = strtok(NULL, ",");
              }
            }
          }
        }
      }
    }
    // skip the rest of the line
		if (dret != '\n') while ((dret = ks_getc(ks)) > 0 && dret != '\n');
    // release the current chr
    if(chr->l > 0)     chr->l = 0;
    if(ref_seq->l > 0) ref_seq->l = 0;
  }

  fprintf(stderr, "Number of mutations: %zu\n", kv_size(tam_records));

  gzFile tam_file = tam_open(output_file, "wb");

  // Set the TAM header
  tam_header = tam_header_init();
  tam_header->k = k_length;
  tam_header->n_kmers = kv_size(tam_records)*2;
  tam_header->n_ref = kv_size(reference_names);
  tam_header->ref = (char**)malloc(sizeof(char**) * tam_header->n_ref);
  for(int i = 0; i < kv_size(reference_names); i++) {
    tam_header->ref[i] = malloc(strlen(kv_A(reference_names,i)) + 1);
    strcpy(tam_header->ref[i], kv_A(reference_names,i));
  }
  tam_header_write(tam_header, tam_file);


  fprintf(stderr, "Load sequences from %s...\n", reference_fasta);
  fp = gzopen(reference_fasta, "r");
  if(!fp) { fprintf(stderr, "Failed to open %s\n", reference_fasta); exit(EXIT_FAILURE); }
  seq = kseq_init(fp);
  while ((l = kseq_read(seq)) >= 0) {
    k = kh_get(references, h_r, seq->name.s);
    if(k != kh_end(h_r)) {
      ref_id = kh_value(h_r, k);
      chr_l = strlen(seq->seq.s);
      for(int i = 0; i < kv_size(tam_records); i++) {
        tam_record = kv_A(tam_records, i);
        if(tam_record && tam_record->ref_id == ref_id && (tam_record->pos - k_length) < chr_l) {
          int ref_l = strlen(tam_record->ref_seq), alt_l = strlen(tam_record->alt_seq);
          uint64_t _for,_rev;

          // Allocate one ref and alt kmers
          tam_record->ref_kmers   = malloc(sizeof(uint64_t));
          tam_record->alt_kmers   = malloc(sizeof(uint64_t));
          tam_record->n_ref_kmers = 1;
          tam_record->n_alt_kmers = 1;

          // First set the ref_kmer
          strncpy(kmer, &seq->seq.s[tam_record->pos - (k_length - strlen(tam_record->ref_seq)) / 2], k_length);
          tam_record->ref_kmers[0] = canonical_kmer(kmer, k_length, &_for, &_rev);

          // The set the alt_kmer
          int before_length = (k_length - alt_l) / 2;
          if(before_length > 0)
            strncpy(kmer, &seq->seq.s[tam_record->pos - before_length], before_length);
          strncpy(&kmer[before_length], tam_record->alt_seq, alt_l);
          before_length += alt_l;
          strncpy(&kmer[before_length], &seq->seq.s[tam_record->pos + ref_l], k_length - before_length);
          tam_record->alt_kmers[0] = canonical_kmer(kmer, k_length, &_for, &_rev);

          k   = kh_get(kmers, h_k, tam_record->ref_kmers[0]);
          k2  = kh_get(kmers, h_k, tam_record->alt_kmers[0]);

          // No matter what, if the alt_kmers already exists, we skip this variant.
          if(k2 != kh_end(h_k)) {
            // If this alt k_mer correspond to two different variants, we
            // removed it, in order to remove the other variant.
            if(kh_value(h_k, k2) != REFERENCE_KMER) {
              kh_del(kmers, h_k, k2);
            }
          } else {
            if(k == kh_end(h_k)) {
              k = kh_put(kmers, h_k, tam_record->ref_kmers[0], &ret);
            }
            kh_value(h_k, k) = REFERENCE_KMER;
            k2 = kh_put(kmers, h_k, tam_record->alt_kmers[0], &ret);
            kh_value(h_k, k2) = MUTATED_KMER;
            tam_record_write(tam_record, tam_file);
          }
          tam_record_destroy(tam_record); // release the memory for this record
          // fprintf(stderr, "%s\t%d\t%s\t%s\n", reference_names.a[tam_record->ref_id], tam_record->pos, tam_record->ref_seq, tam_record->alt_seq);
        }
      }
    }
  }
  gzclose(tam_file);
  kv_destroy(tam_records);
  kv_destroy(reference_names);
  tam_header_destroy(tam_header);

  // fprintf(stderr, "Create mutated k-mer hash...\n");
  // load_kmers(output_file, h_k, 0);
  fprintf(stderr, "Remove mutated k-mers that are also located on the reference (this may take a while)...\n");
  int nb_removed_kmers  = remove_reference_kmers(reference_fasta, h_k, k_length);
  fprintf(stderr, "%d mutated kmers have been removed\n", nb_removed_kmers);
  char *tmp_output_file = get_tmp_filename(output_file);
  fprintf(stderr, "Updating the TAM file...\n");
  int nb_record_removed = update_tam_file(output_file,tmp_output_file,h_k);
  fprintf(stderr, "%d records have been removed\n", nb_record_removed);
  rename(tmp_output_file, output_file);

  return 0;
}

int tami_build(int argc, char *argv[]) {
  char *bed_file = NULL, *reference_fasta = NULL, *output_file = DEFAULT_OUTPUT_NAME;
  int use_derived_kmers = 0;
  int handle_insertions = 0;
  int handle_deletions  = 0;
  int indel_padding     = 4;
  int k_length          = DEFAULT_K_LENGTH;
  int k_middle          = k_length / 2;
  int debug = 0;

  int c;
  while ((c = getopt(argc, argv, "disk:r:o:")) >= 0) {
		switch (c) {
      case 'k': k_length = atoi(optarg); break;
      case 'r': reference_fasta = optarg; break;
      case 's': use_derived_kmers = 1; break;
      case 'o': output_file = optarg; break;
      case 'i': handle_deletions = 1; handle_insertions = 1; break;
      case 'd': debug = 1; break;
		}
	}

  if (optind == argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   tami build [options] <in.bed>\n\n");
		fprintf(stderr, "Options: -k INT    length of k-mers (max_value: 32) [%d]\n", k_length);
    fprintf(stderr, "         -r STR    reference FASTA (otherwise sequences are retrieved from rest.ensembl.org)\n");
    fprintf(stderr, "         -i        support 1-nt indels\n");
    fprintf(stderr, "         -s        sensitive mode (able 2 mutation per-kmer)\n");
    fprintf(stderr, "                   Warning : important memory overload and possible loss of accuracy,\n");
    fprintf(stderr, "                             this option is only adviced for amplicon sequencing.\n");
    fprintf(stderr, "         -o STR    output file (default : 'kmers.tam')\n");
    fprintf(stderr, "         -d        debug mode\n");
		fprintf(stderr, "\n");
		return 1;
	}

  bed_file = argv[optind++];

  // verify Options
  if(k_length > 32 || k_length < 1) {
    fprintf(stderr, "Invalid value (%d) for k_length (-k option).\n", k_length);
    return 1;
  }

  if(debug) {
    fprintf(stderr, "Bed file is: %s\n",bed_file);
  }

  interval_t *interval;
  tam_header_t *tam_header = tam_header_init();
  tam_record_t *tam_record = tam_record_init();
  interval_array_t *interval_array = interval_array_init();
  kvec_t(char*) reference_names;
  khiter_t k, k2;
  khash_t(kmers) *h_k = kh_init(kmers);
  khash_t(references) *h_r = kh_init(references);
  char *tmp_output_file;
  int ret;
  char *kmer = malloc(k_length);

  kv_init(reference_names);

  tmp_output_file = get_tmp_filename(output_file);

  /*****************************************************************************
  *                   READ THE BED FILE AND LOAD INTERVALS
  *****************************************************************************/
  fprintf(stderr, "Reading BED file...\n");
  load_intervals_from_bed(bed_file, interval_array);

  for(int i = 0; i < kv_size(*interval_array); i++) {
    interval = kv_A(*interval_array,i);
    if(kh_get(references, h_r, interval->chr) == kh_end(h_r)) {
      char *chr_cpy = malloc(strlen(interval->chr) + 1);
      strcpy(chr_cpy,interval->chr);
      kv_push(char*,reference_names,chr_cpy);
      k = kh_put(references, h_r, chr_cpy, &ret);
      kh_value(h_r, k) = kv_size(reference_names) - 1;
    }
  }

  merge_overlapping_intervals(interval_array);

  /*****************************************************************************
  *                  LOAD SEQUENCES FROM FASTA OR ENSEMBL API
  *                          (if not FASTA provided)
  *****************************************************************************/

  if(reference_fasta) {
    fprintf(stderr, "Load sequences from %s...\n", reference_fasta);
    load_intervals_seq_from_fasta(reference_fasta, interval_array);
    //kt_for(NB_THREAD_API, get_sequence_iter, interval_array.a, kv_size(interval_array));
  } else {
    fprintf(stderr, "Load sequences from Ensembl REST API...\n");
    kt_for(NB_THREAD_API, get_sequence_iter, interval_array->a, kv_size(*interval_array));
  }

  /*****************************************************************************
  *                      CREATE THE MUTATED K-MER HASH
  *****************************************************************************/

  fprintf(stderr, "Create mutated k-mer hash...\n");

  int ref_kmer_pos, prev_ref_kmer_pos, mut_pos, pos, ref_id;
  char ref_nuc;
  uint64_t ref_kmer, mut_kmer, reverse_mut_kmer, forward_ref_kmer = 0, reverse_ref_kmer = 0;

  /* open file pointer */
  gzFile tam_file = tam_open(output_file, "wb");

  // Write kmer_file header
  tam_header->k = k_length;
  tam_header->n_kmers = 0;
  tam_header->n_ref = kh_size(h_r);
  tam_header->ref = (char**)malloc(sizeof(char**) * tam_header->n_ref);
  for(int i = 0; i < kv_size(reference_names); i++) {
    tam_header->ref[i] = malloc(strlen(kv_A(reference_names,i)) + 1);
    strcpy(tam_header->ref[i], kv_A(reference_names,i));
  }
  tam_header_write(tam_header, tam_file);

  if(!tam_record->ref_kmers) tam_record->ref_kmers   = malloc(sizeof(uint64_t));
  if(!tam_record->alt_kmers) tam_record->alt_kmers   = malloc(sizeof(uint64_t));
  if(!tam_record->ref_seq)   tam_record->ref_seq     = malloc(255);
  if(!tam_record->alt_seq)   tam_record->alt_seq     = malloc(255);
  tam_record->n_ref_kmers = 1;
  tam_record->n_alt_kmers = 1;
  // tam_record->ref_seq_l   = 1;
  // tam_record->alt_seq_l   = 1;

  // Loop over the interval array and print sequences
  for(int i = 0; i < kv_size(*interval_array); i++) {
    interval = kv_A(*interval_array,i);

    if(!interval->seq) {
      fprintf(stderr, "No sequence found for interval %s:%d-%d\n", interval->chr, interval->start, interval->end);
      continue;
    }

    if(debug)
      fprintf(stderr, ">%s:%d..%d\n%s\n", interval->chr, interval->start, interval->end, interval->seq);

    size_t seq_length = strlen(interval->seq);
    if(seq_length < k_length)
      fprintf(stderr, "Sequence length is smaller than k-mer length, skipping.\n");

    ref_kmer = canonical_kmer(interval->seq, k_length, &forward_ref_kmer, &reverse_ref_kmer);
    prev_ref_kmer_pos = 0;
    ref_id = kh_value(h_r, kh_get(references, h_r, interval->chr));

    for (int n = 0; n < seq_length; n++) {

      // Case of the left extremity where the mutated position is not the center of the k-mer
      if(n < k_middle) {
        ref_kmer_pos  = 0;
        mut_pos       = n;
      // Case of the right extremity where the mutated position is not the center of the k-mer
      } else if(n >= seq_length - k_middle) {
        ref_kmer_pos  = seq_length - k_length;
        mut_pos       = k_length + n - seq_length;
      // Default case, the mutated position is the center of the k-mer
      } else {
        ref_kmer_pos  = n - k_middle;
        mut_pos       = k_middle;
      }

      pos = interval->start + ref_kmer_pos + mut_pos;
      ref_nuc = interval->seq[mut_pos + ref_kmer_pos];

      if(ref_kmer_pos > prev_ref_kmer_pos) {
        ref_kmer = next_canonical_kmer(k_length, interval->seq[ref_kmer_pos + k_length - 1], &forward_ref_kmer, &reverse_ref_kmer);
        prev_ref_kmer_pos = ref_kmer_pos;
      }

      tam_record->ref_id        = ref_id;
      tam_record->ref_kmers[0]  = ref_kmer;

      // Handle 1-nt deletion
      // FIXME we should verify that a stretch of nuclotide does not stretch to
      // the end of the k-mer. Same for insertions
      // We should also set the ref_seq and alt_seq to include the whole stretch
      // as in free-bayes
      if(handle_deletions && mut_pos >= indel_padding && mut_pos < (k_length - indel_padding - 1) &&
         interval->seq[mut_pos + ref_kmer_pos] != interval->seq[mut_pos + ref_kmer_pos - 1]) {
        memcpy(kmer,           &interval->seq[ref_kmer_pos],               mut_pos);
        memcpy(&kmer[mut_pos], &interval->seq[ref_kmer_pos + mut_pos + 1], k_length - mut_pos);
        uint64_t _for,_rev;
        mut_kmer = canonical_kmer(kmer, k_length, &_for, &_rev);
        tam_record->pos           = pos - 1;
        if(kh_get(kmers, h_k, mut_kmer) == kh_end(h_k)) {
          k = kh_put(kmers, h_k, mut_kmer, &ret);
          kh_value(h_k, k) = MUTATED_KMER;
          tam_record->alt_kmers[0]  = mut_kmer;
          memcpy(tam_record->ref_seq, &interval->seq[ref_kmer_pos + mut_pos - 1], 2);
          tam_record->ref_seq[2] = '\0';
          memcpy(tam_record->alt_seq, tam_record->ref_seq, 1);
          tam_record->alt_seq[1] = '\0';
          tam_record_write(tam_record, tam_file);
        }
      }

      // Handle 1-nt insertion
      if(handle_insertions && mut_pos >= indel_padding && mut_pos < (k_length - indel_padding - 1)) {
        memcpy(kmer,           &interval->seq[ref_kmer_pos],               mut_pos);
        memcpy(&kmer[mut_pos + 1], &interval->seq[ref_kmer_pos + mut_pos], k_length - mut_pos - 1);
        memcpy(tam_record->ref_seq, &kmer[mut_pos - 1], 1);
        uint64_t _for,_rev;
        tam_record->ref_seq[1] = '\0';
        tam_record->pos           = pos - 1;
        for(int p = 0; p < NB_NUCLEOTIDES; p++) {
          if(NUCLEOTIDES[p] != kmer[mut_pos - 1]) {
            kmer[mut_pos] = NUCLEOTIDES[p];
            // FIXME We could optimize this by using mut_in_dna
            mut_kmer = canonical_kmer(kmer, k_length, &_for, &_rev);
            if(kh_get(kmers, h_k, mut_kmer) == kh_end(h_k)) {
              k = kh_put(kmers, h_k, mut_kmer, &ret);
              kh_value(h_k, k) = MUTATED_KMER;
              tam_record->alt_kmers[0]  = mut_kmer;
              memcpy(tam_record->alt_seq, &kmer[mut_pos - 1], 2);
              tam_record->ref_seq[2] = '\0';
              tam_record_write(tam_record, tam_file);
            }
          }
        }
      }

      // Handle 1-nt substitution
      strncpy(tam_record->ref_seq, &ref_nuc, 1);
      tam_record->ref_seq[1] = '\0';
      tam_record->pos           = pos;

      for(int p = 0; p < NB_NUCLEOTIDES; p++) {
        if(NUCLEOTIDES[p] != ref_nuc) {
          mut_kmer = mut_int_dna(forward_ref_kmer, k_length, mut_pos, NUCLEOTIDES[p]);
          reverse_mut_kmer = int_revcomp(mut_kmer, k_length);
          if(reverse_mut_kmer < mut_kmer)
            mut_kmer = reverse_mut_kmer;

          if(kh_get(kmers, h_k, mut_kmer) == kh_end(h_k)) {
            k = kh_put(kmers, h_k, mut_kmer, &ret);
          	kh_value(h_k, k) = MUTATED_KMER;
            strncpy(tam_record->alt_seq, &NUCLEOTIDES[p], 1);
            tam_record->alt_seq[1] = '\0';
            tam_record->alt_kmers[0]  = mut_kmer;
            tam_record_write(tam_record, tam_file);
          }
        } else {
          if(kh_get(kmers, h_k, ref_kmer) == kh_end(h_k)) {
            k = kh_put(kmers, h_k, ref_kmer, &ret);
          	kh_value(h_k, k) = REFERENCE_KMER;
          }
        }
      }
    }
    //interval_destroy(interval);
  }
  gzclose(tam_file);


  /*****************************************************************************
  *                          ADD 2-nt MUTATED K-MER
  *                             (sensitive mode)
  *****************************************************************************/

  if(use_derived_kmers) {
    // We update the size of the hash to fit the derived kmers. Theoritically, it should
    // kh_size * k * 3, but we in practice it will be closer to a factor 2.
    kh_resize(kmers, h_k, kh_size(h_k)*k_length*2);
    fprintf(stderr, "Create 2-nuc mutated k-mer (-s option)...\n");
    // Add all derived k-mer with one more mutation
    tam_file = tam_open(output_file, "rb");
    gzFile tam_file_alt = tam_open(tmp_output_file, "wb");

    if(tam_header_read(tam_header,tam_file)) {

      tam_header->n_kmers = kh_size(h_k);
      tam_header_write(tam_header, tam_file_alt);

      while(tam_record_read(tam_record,tam_file)) {
        if(tam_record->n_alt_kmers == 1) {
          tam_record->alt_kmers = realloc(tam_record->alt_kmers, sizeof(uint64_t) * k_length * (NB_NUCLEOTIDES - 1) + 1);

          for (int n = 0; n < k_length; n++) {
            ref_nuc = nuc_from_int_dna(tam_record->alt_kmers[0], k_length, n);
            for(int p = 0; p < NB_NUCLEOTIDES; p++) {
              if(ref_nuc != NUCLEOTIDES[p]) {
                mut_kmer = mut_int_dna(tam_record->alt_kmers[0], k_length, n,  NUCLEOTIDES[p]);
                if(kh_get(kmers, h_k, mut_kmer) == kh_end(h_k)) {
                  k2 = kh_put(kmers, h_k, mut_kmer, &ret);
                  kh_value(h_k, k2) = DOUBLE_MUTATED_KMER;
                  tam_record->alt_kmers[tam_record->n_alt_kmers++] = mut_kmer;
                }
              }
            }
          }
          tam_record_write(tam_record,tam_file_alt);
        }
      }
    }
    gzclose(tam_file);
    gzclose(tam_file_alt);
    rename(tmp_output_file, output_file);
  }

  /*****************************************************************************
  *             REMOVE MUTATED K-MER ALSO FOUND ON THE REFERENCE
  *****************************************************************************/
  if(reference_fasta) {
    fprintf(stderr, "Remove mutated k-mers that are also located on the reference (this may take a while)...\n");
    int nb_removed_kmers = remove_reference_kmers(reference_fasta, h_k, k_length);
    fprintf(stderr, "%d mutated kmers have been removed\n", nb_removed_kmers);

    fprintf(stderr, "Updating the TAM file...\n");

    int nb_removed_records = update_tam_file(output_file, tmp_output_file, h_k);
    fprintf(stderr, "%d records have been removed\n", nb_removed_records);
    rename(tmp_output_file, output_file);
  }

  kv_destroy(*interval_array);
  kv_destroy(reference_names);
  tam_header_destroy(tam_header);
  tam_record_destroy(tam_record);
  //interval_array_destroy(interval_array); // FIXME THERE IS A PROBLEM WITH A DOUBLE FREE
	return 0;
}

int tami_scan(int argc, char *argv[]) {
  char *fastq_file, *tam_path;
  int debug = 0;
  float min_alternate_fraction = 0.2;
  int min_alternate_count = 3;
  int min_coverage = 10;
  int max_coverage = -1;

  int c;
  while ((c = getopt(argc, argv, "dF:C:m:M:")) >= 0) {
    switch (c) {
      case 'F': min_alternate_fraction = atof(optarg); break;
      case 'C': min_alternate_count = atoi(optarg); break;
      case 'm': min_coverage = atoi(optarg); break;
      case 'M': max_coverage = atoi(optarg); break;
      case 'd': debug = 1; break;
    }
  }

  if (optind + 1 >= argc) {
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage:   tami build [options] <in.tam> <in.fq>\n\n");
    fprintf(stderr, "Options: -F FLOAT  min alternate allele frequency [%.2f]\n", min_alternate_fraction);
    fprintf(stderr, "         -C INT    min alternate allele count (min_value: 1)[%d]\n", min_alternate_count);
    fprintf(stderr, "         -m INT    min coverage [%d]\n", min_coverage);
    fprintf(stderr, "         -M INT    max coverage (min_value: 1)[%d]\n", max_coverage);
    fprintf(stderr, "         -d        debug mode\n");
    fprintf(stderr, "\n");
    return 1;
  }

  tam_path = argv[optind++];

  // verify Options
  if(min_coverage < 1) {
    fprintf(stderr, "Invalid value (%d) for min_coverage (-m option).\n", min_coverage);
    return 1;
  }
  if(min_alternate_count < 1) {
    fprintf(stderr, "Invalid value (%d) for min_alternate_count (-C option).\n", min_coverage);
    return 1;
  }

  if(debug) {
    fprintf(stderr, "min_alternate_count (C): %d\n",min_alternate_count);
    fprintf(stderr, "min_alternate_fraction (F): %f\n",min_alternate_fraction);
  }

  tam_header_t *tam_header = tam_header_init();
  tam_record_t *tam_record = tam_record_init();
  gzFile  tam_file = tam_open(tam_path, "rb");
  khiter_t k, k2;
  khash_t(kmers) *h_k = kh_init(kmers);
  char *kmer;
  int k_length;
  gzFile fp;

  tam_header_read(tam_header,tam_file);
  k_length = tam_header->k;
  kmer = malloc(k_length + 1);

  fprintf(stderr, "Loading kmers from TAM file into memory...\n");
  int nb_kmers = load_kmers(tam_path, h_k, 1);
  fprintf(stderr, "%d k-mers loaded into memory\n", nb_kmers);
  gzclose(tam_file);

  /*****************************************************************************
  *                             SCAN FASTQ FILES
  *****************************************************************************/

  // Set all counters to 0
  // FIXME This should not be necessary given the load_kmers parameters
  for(k = kh_begin(h_k); k != kh_end(h_k); ++k) {
    if (!kh_exist(h_k, k)) continue;
    kh_value(h_k, k) = 0;
  }

  while(optind < argc) {
    fastq_file = argv[optind++];
    fprintf(stderr, "Reading %s FASTQ file...\n", fastq_file);

    kseq_t *seq;
    int l;
    uint64_t kmer_int, forward_kmer_int = 0, reverse_kmer_int = 0;
    int nb_reads = 0;
    fp = gzopen(fastq_file, "r");
    seq = kseq_init(fp);
    while ((l = kseq_read(seq)) >= 0) {
      nb_reads++;
      if(seq->seq.l >= k_length) {
        for(int i = 0; i < seq->seq.l - k_length + 1; i++) {
          if(i > 0) {
            kmer_int = next_canonical_kmer(k_length, seq->seq.s[i + k_length - 1], &forward_kmer_int, &reverse_kmer_int);
          } else {
            kmer_int = canonical_kmer(seq->seq.s, k_length, &forward_kmer_int, &reverse_kmer_int);
          }
          k = kh_get(kmers, h_k, kmer_int);
          if(k != kh_end(h_k)) {
            kh_value(h_k, k) = kh_value(h_k, k) + 1;
          }
        }
      }
      if (nb_reads % 50000 == 0) {
        fprintf(stderr, "*");
      }
    }
    fprintf(stderr, "\n%d reads parsed.\n",nb_reads);
    kseq_destroy(seq);
    gzclose(fp);
  }

  /*****************************************************************************
  *                             UPDATE COUNTS
  *****************************************************************************/

  fprintf(stderr, "Update counts...\n");
  tam_file = tam_open(tam_path, "rb");
  tam_header_read(tam_header,tam_file);
  while(tam_record_read(tam_record,tam_file)) {
    for(int i = 0; i < tam_record->n_ref_kmers; i++) {
      k = kh_get(kmers, h_k, tam_record->ref_kmers[i]);
      if(k != kh_end(h_k)) {
        for(int j = 0; j < tam_record->n_alt_kmers; j++) {
          k2 = kh_get(kmers, h_k, tam_record->alt_kmers[j]);
          if(k2 != kh_end(h_k))
            kh_value(h_k, k) = kh_value(h_k, k) + kh_value(h_k, k2);
        }
      }
    }
  }
  gzclose(tam_file);

  // Compute mean coverage
  // int x = 0;
  // int y = 0;
  // int N = 0;
  //
  // tam_file = tam_open(tam_path, "rb");
  // tam_header_read(tam_header,tam_file);
  // while(tam_record_read(tam_record,tam_file)) {
  //   for(int i = 0; i < tam_record->n_ref_kmers; i++) {
  //     if(tam_record->ref_kmers[i] > 0) {
  //       N++;
  //       break;
  //     }
  //   }
  // }
  // gzclose(tam_file);
  //
  // tam_file = tam_open(tam_path, "rb");
  // tam_header_read(tam_header,tam_file);
  // while(tam_record_read(tam_record,tam_file)) {
  //   int DP = 0;
  //   for(int i = 0; i < tam_record->n_ref_kmers; i++) {
  //     k = kh_get(kmers, h_k, tam_record->ref_kmers[i]);
  //     if(k != kh_end(h_k)) {
  //       DP += kh_val(h_k, k);
  //     }
  //   }
  //   if(DP > 0) {
  //     x += DP / N;
  //     int b = DP % N;
  //     if (y >= N - b) {
  //       x++;
  //       y -= N - b;
  //     } else {
  //       y += b;
  //     }
  //   }
  // }
  // gzclose(tam_file);
  //
  // double mean_coverage =  x + y / N;
  // fprintf(stderr, "Mean coverage: %.2f\n", mean_coverage);
  // double negative_mean_coverage_exp = exp(-1 * mean_coverage);

  /*****************************************************************************
  *                             PRINT OUTPUT
  *****************************************************************************/

  /* FILTER AND PRINT OUTPUT VCF */
  fprintf(stdout, "##fileformat=VCFv4.1\n");
  fprintf(stdout, "##source=TaMI v%s\n", TAMI_VERSION);
  fprintf(stdout, "##commandline=");
  for (int i = 0; i < argc; i++)
    fprintf(stdout, " %s", argv[i]);
  fprintf(stdout, "\n##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total read depth at the locus\">\n");
  fprintf(stdout, "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Estimated allele frequency in the range (0,1]\">\n");
  fprintf(stdout, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
  fprintf(stdout, "##FORMAT=<ID=RO,Number=1,Type=Integer,Description=\"Reference allele observation count\">\n");
  fprintf(stdout, "##FORMAT=<ID=AO,Number=A,Type=Integer,Description=\"Alternate allele observation count\">\n");
  fprintf(stdout, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tNA00001\n");

  tam_file = tam_open(tam_path, "rb");
  tam_header_read(tam_header,tam_file);
  while(tam_record_read(tam_record,tam_file)) {
    int AC = 0;
    int DP = 0;
    int max_alt_kmer_value = 0;
    uint64_t max_alt_kmer = 0;
    for(int i = 0; i < tam_record->n_ref_kmers; i++) {
      k = kh_get(kmers, h_k, tam_record->ref_kmers[i]);
      if(k != kh_end(h_k)) {
        DP += kh_val(h_k, k);
      }
    }
    for(int i = 0; i < tam_record->n_alt_kmers; i++) {
      k   = kh_get(kmers, h_k, tam_record->alt_kmers[i]);
      if(k != kh_end(h_k)) {
        AC += kh_val(h_k, k);
        if(kh_val(h_k, k) > max_alt_kmer_value) {
          max_alt_kmer_value = kh_val(h_k, k);
          max_alt_kmer = tam_record->alt_kmers[i];
        }
      }
    }

    if(AC < min_alternate_count)
      continue;
    if(DP < min_coverage)
      continue;
    float AF = (float) AC / DP;
    if(AF < min_alternate_fraction)
      continue;
    if(max_coverage > 0 && DP > max_coverage)
      continue;

    // Compute genotypes see LAVA paper (Shajii & al 2016)
    char GT[4] = "0/0";
    double p_error = 0.01;
    long long combi_a = combi(DP, DP - AC);
    long long combi_b = combi(DP, AC);
    double p_C_g0 = (double) combi_a * pow(1 - p_error, DP - AC) * pow(p_error, AC);
    double p_C_g1 = (double) combi_a * 1/pow(2, DP);
    double p_C_g2 = (double) combi_b * pow(p_error, DP - AC) * pow(1 - p_error, AC);
    double p_C    = p_C_g0 + p_C_g1 + p_C_g2;
    double p_g0   = pow((double) (DP - AC) / DP, 2);
    double p_g2   = pow((double) AC / DP, 2);
    double p_g1   = 1 - p_g0 - p_g2;
    double p_g0_C = p_g0 * p_C_g0 / p_C;
    double p_g1_C = p_g1 * p_C_g1 / p_C;
    double p_g2_C = p_g2 * p_C_g2 / p_C;
    // double scaling_term = pow(mean_coverage, DP) / factorial(DP) * negative_mean_coverage_exp;
    // double score, p;

    if(p_g1_C > p_g0_C && p_g1_C > p_g2_C) {
      strcpy(GT,"0/1");
      // p = p_g1_C;
    } else if(p_g2_C > p_g0_C && p_g2_C > p_g1_C) {
      strcpy(GT,"1/1");
      // p = p_g2_C;
    } else {
      // p = p_g0_C;
    }

    // score = p * scaling_term;
    // FIXME RO is not computed right, because of the DP cumultate count for ref and other alterate alleles...
    int_to_dna(max_alt_kmer,k_length,kmer);
    fprintf(stdout, "%s\t%d\t%s\t%s\t%s\t.\tPASS\tDP=%d;AF=%.2f\tGT:RO:AO\t%s:%d:%d\n", tam_header->ref[tam_record->ref_id], tam_record->pos + 1, kmer, tam_record->ref_seq, tam_record->alt_seq, DP, AF, GT, DP - AC, AC);

  }
  gzclose(tam_file);

  tam_header_destroy(tam_header);
  tam_record_destroy(tam_record);
  return 0;
}

int tami_view(int argc, char *argv[]) {
  if (argc < 2) {
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage:   tami view <in.tam>\n");
    fprintf(stderr, "\n");
    return 1;
  }

  char *tam_path = argv[1];
  tam_header_t *tam_header = tam_header_init();
  tam_record_t *tam_record = tam_record_init();
  char *kmer;
  gzFile  tam_file = tam_open(tam_path, "rb");

  tam_header_read(tam_header,tam_file);
  kmer = malloc(tam_header->k + 1);
  kmer[tam_header->k] = '\0';

  while(tam_record_read(tam_record,tam_file)) {
    tam_record_print(stdout, tam_header, tam_record);
  }
  gzclose(tam_file);

  tam_header_destroy(tam_header);
  tam_record_destroy(tam_record);
  return 0;
}


/* main function */
static int usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage:   tami <command> <arguments>\n");
	fprintf(stderr, "Version: %s\n\n", TAMI_VERSION);
	fprintf(stderr, "Command: build      Create a mutated k-mer lib from genomic intervals\n");
  fprintf(stderr, "         vcf-build  Create a mutated k-mer lib from known mutation\n");
	fprintf(stderr, "         scan       Scan a FASTQ files agains a TAM file\n");
  fprintf(stderr, "         view       Display text version of TAM files\n");
	fprintf(stderr, "\n");
	return 1;
}

int main(int argc, char *argv[])
{
	if (argc == 1) return usage();
	if (strcmp(argv[1], "build") == 0) tami_build(argc-1, argv+1);
  else if (strcmp(argv[1], "vcf-build") == 0) tami_vcfbuild(argc-1, argv+1);
	else if (strcmp(argv[1], "scan") == 0) tami_scan(argc-1, argv+1);
  else if (strcmp(argv[1], "view") == 0) tami_view(argc-1, argv+1);
	else {
		fprintf(stderr, "[main] unrecognized command '%s'. Abort!\n", argv[1]);
		return 1;
	}
	return 0;
}
