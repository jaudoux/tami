#include <stdio.h>  //fprintf
#include <stdlib.h> //free
#include <zlib.h>
#include <inttypes.h>
#include <math.h> // pow()
#include <stdint.h>

#include "kstring.h"
#include "ksort.h"
#include "khash.h"
#include "kvec.h"
#include "kthread.h"
#include "kseq.h"

#include "api.h"
#include "dna.h"
#include "tam.h"
#include "intervals.h"

#define TAMI_VERSION "0.3.0"
#define NB_THREAD_API 10
#define DEFAULT_OUTPUT_NAME "kmers.tam"
#define DEFAULT_K_LENGTH 32

#define REFERENCE_KMER_CHECKED 0
#define REFERENCE_KMER 1
#define MUTATED_KMER 2
#define DOUBLE_MUTATED_KMER 3

KSEQ_INIT(gzFile, gzread)
KHASH_MAP_INIT_INT64(kmers, uint32_t)
KHASH_MAP_INIT_STR(references, int32_t)

typedef khash_t(kmers) kmers_hash_t;

typedef struct {
  uint32_t count : 31;
  uint32_t is_reference_kmer : 1;
  uint32_t target_id;
  uint32_t pos;
} kmer_count_t;

KHASH_MAP_INIT_INT64(kmers_count, kmer_count_t*)
typedef khash_t(kmers_count) kmers_count_hash_t;

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
    if(!kh_exist(h, k)) continue;
    if(kh_value(h, k) == REFERENCE_KMER_CHECKED)
      kh_value(h, k) = REFERENCE_KMER;
  }

  return nb_removed_kmers;
}

int load_kmers(char *tam_path, kmers_count_hash_t *kch) {
  int ret;
  khiter_t k;
  gzFile tam_file = tam_open(tam_path, "rb");
  tam_header_t *h = tam_header_init();
  tam_record_t *r = tam_record_init();
  uint32_t mut_id = 0;
  if(tam_header_read(h,tam_file)) {
    kh_resize(kmers_count, kch, kh_size(kch) + h->n_kmers);
    // Load ref_kmers
    for(int i = 0; i < h->n_target; i++) {
      tam_target_t *t = h->target[i];
      for(int j = 0; j < t->length - h->k + 1; j++) {
        if((t->ref_kmers_bv[j / 8] >> j % 8) & 1) {
          uint64_t f_k, r_k;
          uint64_t kmer_int = canonical_kmer(&t->seq[j],h->k,&f_k,&r_k);
          k = kh_get(kmers_count, kch, kmer_int);
          if(k == kh_end(kch)) {
            kmer_count_t *kc = (kmer_count_t*)calloc(1, sizeof(kmer_count_t));
            kc->is_reference_kmer = 1;
            kc->target_id         = i;
            kc->pos               = j;
            k = kh_put(kmers_count, kch, kmer_int, &ret);
            kh_value(kch, k) = kc;
          }
        }
      }
    }

    // Load mutated k-mers
    while(tam_record_read(r,tam_file)) {
      // First check if the mutated k-mers overlaps unique reference k-mers
      int found_ref_kmer = 0, i = 0, j = 0;
      tam_target_t *t = h->target[r->target_id];
      while(!found_ref_kmer && i < h->k) {
        j = r->pos - i;
        // If the k-mer ref bit is up, we have at least one overlapping ref-kmer
        if(j >= 0 && (t->ref_kmers_bv[j / 8] >> j % 8) & 1) {
          found_ref_kmer = 1;
          // FIXME We should do something better than this, because we could know the exact
          // Position of the k_mer in the genome
          kmer_count_t *kc = (kmer_count_t*)calloc(1, sizeof(kmer_count_t));
          kc->is_reference_kmer = 0;
          kc->target_id         = r->target_id;
          kc->pos               = r->pos;
          for(int i = 0; i < r->n_alt_kmers; i++) {
            k = kh_put(kmers_count, kch, r->alt_kmers[i], &ret);
            kh_value(kch, k) = kc;
          }
          mut_id++;
        }
        i++;
      }
    }
  }
  tam_header_destroy(h);
  tam_record_destroy(r);
  return (int)kh_size(kch);
}

int update_tam_file(const char *tam_path, const char *new_tam_path, kmers_hash_t *hk) {
  khiter_t k;
  gzFile f        = tam_open(tam_path, "rb");
  gzFile f_n      = tam_open(new_tam_path, "wb");
  tam_header_t *h = tam_header_init();
  tam_record_t *r = tam_record_init();
  int nb_record_removed = 0;

  if(tam_header_read(h,f)) {
    h->n_kmers = kh_size(hk);
    // First check reference k-mers
    for(int i = 0; i < h->n_target; i++) {
      tam_target_t *t = h->target[i];

      uint64_t ref_kmer, f_k, r_k;
      ref_kmer = canonical_kmer(t->seq, h->k, &f_k, &r_k);
      for (int n = 0; n < t->length - h->k + 1; n++) {
        if(n > 0) {
          ref_kmer = next_canonical_kmer(h->k, t->seq[n + h->k - 1], &f_k, &r_k);
        }
        k = kh_get(kmers, hk, ref_kmer);
        if(k != kh_end(hk) && kh_value(hk, k) == REFERENCE_KMER) {
          t->ref_kmers_bv[n / 8] |= 1 << n % 8;
        } else {
          t->ref_kmers_bv[n / 8] &= ~(1 << n % 8);
        }
      }
    }

    tam_header_write(h, f_n);

    // Now verify alt kmers
    while(tam_record_read(r,f)) {
      int j = 0;
      for(int i = 0; i < r->n_alt_kmers; i++) {
        k = kh_get(kmers, hk, r->alt_kmers[i]);
        if(k != kh_end(hk) && kh_value(hk, k) != REFERENCE_KMER)
          r->alt_kmers[j++] = r->alt_kmers[i];
      }
      r->n_alt_kmers = j;

      // TODO we should also remove the record if it has no reference k-mers
      if(r->n_alt_kmers == 0)
        nb_record_removed++;
      else
        tam_record_write(r, f_n);
    }
  }

  gzclose(f);
  gzclose(f_n);
  tam_header_destroy(h);
  tam_record_destroy(r);
  return nb_record_removed;
}

int tami_build(int argc, char *argv[]) {
  char *bed_file = NULL, *reference_fasta = NULL, *output_file = DEFAULT_OUTPUT_NAME;
  // int use_derived_kmers = 0;
  int handle_insertions = 0;
  int handle_deletions  = 0;
  int indel_padding     = 4;
  int kmer_sampling     = 3;
  int k_length          = DEFAULT_K_LENGTH;
  int target_padding    = k_length;
  int debug = 0;

  int c;
  while ((c = getopt(argc, argv, "disk:r:o:n:")) >= 0) {
		switch (c) {
      case 'k': k_length = atoi(optarg); break;
      case 'r': reference_fasta = optarg; break;
      // case 's': use_derived_kmers = 1; break;
      case 'o': output_file = optarg; break;
      case 'n': kmer_sampling = atoi(optarg); break;
      case 'i': handle_deletions = 1; handle_insertions = 1; break;
      case 'd': debug = 1; break;
		}
	}

  if (optind == argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   tami build [options] <in.bed>\n\n");
		fprintf(stderr, "Options: -k INT    length of k-mers (max_value: 32) [%d]\n", k_length);
    fprintf(stderr, "         -r STR    reference FASTA (otherwise sequences are retrieved from rest.ensembl.org)\n");
    fprintf(stderr, "         -n INT    number of overlapping k-mer [%d]\n", kmer_sampling);
    fprintf(stderr, "         -p INT    target padding [%d]\n", target_padding);
    fprintf(stderr, "         -i        support 1-nt indels\n");
    // fprintf(stderr, "         -s        sensitive mode (able 2 mutation per-kmer)\n");
    // fprintf(stderr, "                   Warning : important memory overload and possible loss of accuracy,\n");
    // fprintf(stderr, "                             this option is only adviced for amplicon sequencing.\n");
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

  int k_middle          = k_length / 2;
  int k_jump            = k_length / (kmer_sampling - 1);

  interval_t *interval;
  tam_header_t *tam_header = tam_header_init();
  tam_record_t *tam_record = tam_record_init();
  interval_array_t *interval_array = interval_array_init();
  kvec_t(char*) reference_names;
  khiter_t k;//, k2;
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

  int nb_merged = merge_overlapping_intervals(interval_array);
  fprintf(stderr, "%d overlapping intervals have been merged\n", nb_merged);

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
  int ref_kmer_pos, mut_pos;
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
  kv_destroy(reference_names);

  tam_header->n_target = 0;
  tam_header->target   = malloc(kv_size(*interval_array) * sizeof(tam_target_t*));
  for(int i = 0; i < kv_size(*interval_array); i++) {
    interval = kv_A(*interval_array,i);
    if(!interval->seq) {
      fprintf(stderr, "No sequence found for interval %s:%d-%d\n", interval->chr, interval->start, interval->end);
      continue;
    }
    if(interval_length(interval) < k_length) {
      fprintf(stderr, "Sequence length is smaller than k-mer length, skipping.\n");
      continue;
    }
    tam_target_t *t = tam_target_init();
    t->ref_id = kh_value(h_r, kh_get(references, h_r, interval->chr));
    t->pos    = interval->start;
    t->length = interval_length(interval);
    t->seq    = interval->seq;
    t->ref_kmers_bv = calloc(t->length/8+1,8);
    tam_header->target[tam_header->n_target++] = t;
    interval->seq = NULL; // release seq string given to the target
  }
  interval_array_destroy(interval_array);

  tam_header_write(tam_header, tam_file);

  tam_record->alt_kmers   = malloc(sizeof(uint64_t)*kmer_sampling);
  tam_record->alt_seq     = malloc(k_length);

  for(int i = 0; i < tam_header->n_target; i++) {
    tam_target_t *t = tam_header->target[i];
    // Compute all ref kmers
    uint64_t *ref_kmers = malloc(sizeof(uint64_t) * (t->length - k_length + 1));
    uint64_t *forward_ref_kmers = malloc(sizeof(uint64_t) * (t->length - k_length + 1));
    ref_kmer = canonical_kmer(t->seq, k_length, &forward_ref_kmer, &reverse_ref_kmer);
    ref_kmers[0] = ref_kmer;
    forward_ref_kmers[0] = forward_ref_kmer;
    for (int n = 0; n < t->length - k_length + 1; n++) {
      if(n > 0) {
        ref_kmers[n] = next_canonical_kmer(k_length, t->seq[n + k_length - 1], &forward_ref_kmer, &reverse_ref_kmer);
        forward_ref_kmers[n] = forward_ref_kmer;
      }
      if(kh_get(kmers, h_k, ref_kmers[n]) == kh_end(h_k)) {
        k = kh_put(kmers, h_k, ref_kmers[n], &ret);
        kh_value(h_k, k) = REFERENCE_KMER;
      }
    }
    tam_record->target_id = i;

    // Loop over all positions in the sequences
    for (int n = 0; n < t->length; n++) {

      // Set ref_kmer_pos and mut_pos for indels
      if(n >= k_middle) {
        if (n < t->length - k_middle) {
          ref_kmer_pos = n - k_middle;
          mut_pos      = k_middle;
        } else {
          ref_kmer_pos = t->length - k_length;
          mut_pos      = k_length + n - t->length;
        }
      } else {
        ref_kmer_pos = 0;
        mut_pos      = n;
      }

      tam_record->n_alt_kmers   = 1;

      // Handle 1-nt deletion
      // FIXME we should verify that a stretch of nuclotide does not stretch to
      // the end of the k-mer. Same for insertions
      // We should also set the ref_seq and alt_seq to include the whole stretch
      // as in free-bayes
      if(handle_deletions && n >= indel_padding && n < (t->length - indel_padding - 1) &&
      t->seq[n] != t->seq[n - 1]) {
        memcpy(kmer,           &t->seq[ref_kmer_pos],               mut_pos);
        memcpy(&kmer[mut_pos], &t->seq[ref_kmer_pos + mut_pos + 1], k_length - mut_pos);
        uint64_t _for,_rev;
        mut_kmer = canonical_kmer(kmer, k_length, &_for, &_rev);
        tam_record->pos           = n - 1;
        tam_record->ref_seq_l     = 2;
        tam_record->alt_seq_l     = 1;
        if(kh_get(kmers, h_k, mut_kmer) == kh_end(h_k)) {
          k = kh_put(kmers, h_k, mut_kmer, &ret);
          kh_value(h_k, k) = MUTATED_KMER;
          tam_record->alt_kmers[0]  = mut_kmer;
          tam_record->alt_seq[0] = t->seq[ref_kmer_pos + mut_pos];
          tam_record->alt_seq[1] = '\0';
          tam_record_write(tam_record, tam_file);
        }
      }

      // Handle 1-nt insertion
      if(handle_insertions && mut_pos >= indel_padding && n < (t->length - indel_padding - 1)) {
        memcpy(kmer,           &t->seq[ref_kmer_pos],               mut_pos);
        memcpy(&kmer[mut_pos + 1], &t->seq[ref_kmer_pos + mut_pos], k_length - mut_pos - 1);
        uint64_t _for,_rev;
        tam_record->pos           = n - 1;
        tam_record->ref_seq_l     = 1;
        tam_record->alt_seq_l     = 2;
        for(int p = 0; p < NB_NUCLEOTIDES; p++) {
          if(NUCLEOTIDES[p] != kmer[mut_pos - 1]) {
            kmer[mut_pos] = NUCLEOTIDES[p];
            // FIXME We could optimize this by using mut_in_dna
            mut_kmer = canonical_kmer(kmer, k_length, &_for, &_rev);
            if(kh_get(kmers, h_k, mut_kmer) == kh_end(h_k)) {
              k = kh_put(kmers, h_k, mut_kmer, &ret);
              kh_value(h_k, k) = MUTATED_KMER;
              tam_record->alt_kmers[0]  = mut_kmer;
              strncpy(tam_record->alt_seq, &kmer[mut_pos - 1], 2);
              tam_record_write(tam_record, tam_file);
            }
          }
        }
      }

      ref_nuc = t->seq[n];

      // Handle 1-nt substitutions
      tam_record->pos       = n;
      tam_record->ref_seq_l = 1;
      tam_record->alt_seq_l = 1;

      for(int p = 0; p < NB_NUCLEOTIDES; p++) {
        if(NUCLEOTIDES[p] != ref_nuc) {
          strncpy(tam_record->alt_seq, &NUCLEOTIDES[p], 1);
          tam_record->alt_seq[1] = '\0';
          tam_record->n_alt_kmers = 0;

          // Loop over all overlapping k-mers given the sampling rate
          for (int j = 0; j <= k_length; j += k_jump) {

            // Set the mut_pos to the last nucleotide of the k_mer if j is larger
            // than k_length. This happens for the last k-mer when k is even
            if(j == k_length) {
              ref_kmer_pos = n - k_length + 1;
              mut_pos      = k_length - 1;
            } else {
              ref_kmer_pos  = n - j;
              mut_pos       = j;
            }

            // The k-mer is not completely included in the reference sequence
            if(ref_kmer_pos >= 0 && ref_kmer_pos + k_length <= t->length) {
              ref_kmer                    = ref_kmers[ref_kmer_pos];
              forward_ref_kmer            = forward_ref_kmers[ref_kmer_pos];

              // Mut the ref k-mer
              mut_kmer = mut_int_dna(forward_ref_kmer, k_length, mut_pos, NUCLEOTIDES[p]);
              reverse_mut_kmer = int_revcomp(mut_kmer, k_length);
              if(reverse_mut_kmer < mut_kmer)
                mut_kmer = reverse_mut_kmer;

              // Add the mut k-mer if it is not already in the hash
              k = kh_get(kmers, h_k, mut_kmer);
              if(k == kh_end(h_k)) {
                k = kh_put(kmers, h_k, mut_kmer, &ret);
                kh_value(h_k, k) = MUTATED_KMER;
                tam_record->alt_kmers[tam_record->n_alt_kmers++]  = mut_kmer;
              } else {
                kh_del(kmers, h_k, k);
              }
            }
          }
          // Print the record in the tam file if it has at least one ref and
          // on alt k-mer
          if(tam_record->n_alt_kmers > 0)
            tam_record_write(tam_record, tam_file);
        }
      }
    }
    free(ref_kmers);
    free(forward_ref_kmers);
  }
  gzclose(tam_file);


  // /*****************************************************************************
  // *                          ADD 2-nt MUTATED K-MER
  // *                             (sensitive mode)
  // *****************************************************************************/
  //
  // if(use_derived_kmers) {
  //   // We update the size of the hash to fit the derived kmers. Theoritically, it should
  //   // kh_size * k * 3, but we in practice it will be closer to a factor 2.
  //   kh_resize(kmers, h_k, kh_size(h_k)*k_length*2);
  //   fprintf(stderr, "Create 2-nuc mutated k-mer (-s option)...\n");
  //   // Add all derived k-mer with one more mutation
  //   tam_file = tam_open(output_file, "rb");
  //   gzFile tam_file_alt = tam_open(tmp_output_file, "wb");
  //
  //   if(tam_header_read(tam_header,tam_file)) {
  //
  //     tam_header->n_kmers = kh_size(h_k);
  //     tam_header_write(tam_header, tam_file_alt);
  //
  //     while(tam_record_read(tam_record,tam_file)) {
  //       if(tam_record->n_alt_kmers == 1) {
  //         tam_record->alt_kmers = realloc(tam_record->alt_kmers, sizeof(uint64_t) * k_length * (NB_NUCLEOTIDES - 1) + 1);
  //
  //         for (int n = 0; n < k_length; n++) {
  //           ref_nuc = nuc_from_int_dna(tam_record->alt_kmers[0], k_length, n);
  //           for(int p = 0; p < NB_NUCLEOTIDES; p++) {
  //             if(ref_nuc != NUCLEOTIDES[p]) {
  //               mut_kmer = mut_int_dna(tam_record->alt_kmers[0], k_length, n,  NUCLEOTIDES[p]);
  //               if(kh_get(kmers, h_k, mut_kmer) == kh_end(h_k)) {
  //                 k2 = kh_put(kmers, h_k, mut_kmer, &ret);
  //                 kh_value(h_k, k2) = DOUBLE_MUTATED_KMER;
  //                 tam_record->alt_kmers[tam_record->n_alt_kmers++] = mut_kmer;
  //               }
  //             }
  //           }
  //         }
  //         tam_record_write(tam_record,tam_file_alt);
  //       }
  //     }
  //   }
  //   gzclose(tam_file);
  //   gzclose(tam_file_alt);
  //   rename(tmp_output_file, output_file);
  // }
  //

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

  tam_header_destroy(tam_header);
  tam_record_destroy(tam_record);
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
  khiter_t k;
  khash_t(kmers_count) *h_k = kh_init(kmers_count);
  char *kmer;
  int k_length;
  gzFile fp;

  tam_header_read(tam_header,tam_file);
  k_length = tam_header->k;
  kmer = malloc(k_length + 1);
  gzclose(tam_file);

  fprintf(stderr, "Loading kmers from TAM file into memory...\n");
  int nb_kmers = load_kmers(tam_path, h_k);
  fprintf(stderr, "%d k-mers loaded into memory\n", nb_kmers);

  /*****************************************************************************
  *                             SCAN FASTQ FILES
  *****************************************************************************/

  while(optind < argc) {
    fastq_file = argv[optind++];
    fprintf(stderr, "Reading %s FASTQ file...\n", fastq_file);

    kseq_t *seq;
    int l;
    uint64_t kmer_int, f_k = 0, r_k = 0;
    int nb_reads = 0;
    fp = gzopen(fastq_file, "r");
    seq = kseq_init(fp);
    while ((l = kseq_read(seq)) >= 0) {
      uint32_t prev_target_id = UINT32_MAX, prev_pos = UINT32_MAX;
      uint32_t prev_is_reference_kmer = 0;
      nb_reads++;
      if(seq->seq.l >= k_length) {
        for(int i = 0; i < seq->seq.l - k_length + 1; i++) {
          if(i > 0) {
            kmer_int = next_canonical_kmer(k_length, seq->seq.s[i + k_length - 1], &f_k, &r_k);
          } else {
            kmer_int = canonical_kmer(seq->seq.s, k_length, &f_k, &r_k);
          }
          k = kh_get(kmers_count, h_k, kmer_int);
          if(k != kh_end(h_k)) {
            kmer_count_t *kc = (kmer_count_t*)kh_value(h_k, k);
            // TODO we could add some "pseudo-mapping" procedure to only count mutated k-mers
            // if a reference_kmer for mapping to the same loci is found in the read
            if(kc->is_reference_kmer) {
              // Only update the k_mer if the previous matched one was further that k
              if(!prev_is_reference_kmer || (prev_is_reference_kmer && (prev_target_id != kc->target_id || abs(kc->pos - prev_pos) >= k_length))) {
                kc->count++;
                prev_is_reference_kmer = 1;
                prev_target_id         = kc->target_id;
                prev_pos               = kc->pos;
              }
              // Only update the k-mer if the previous matched one was not the same mutation
            } else if(prev_target_id != kc->target_id || kc->pos != prev_pos) {
              kc->count++;
              prev_is_reference_kmer = 0;
              prev_target_id         = kc->target_id;
              prev_pos               = kc->pos;
            }
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

  // Compute coverage for each target
  uint32_t **target_coverage = calloc(tam_header->n_target, sizeof(uint32_t*));
  for(int i = 0; i < tam_header->n_target; i++) {
    tam_target_t *t = tam_header->target[i];
    target_coverage[i] = calloc(t->length,sizeof(uint32_t));
    uint64_t f_k, r_k;
    for(int j = 0; j < t->length - tam_header->k + 1; j++) {
      if((t->ref_kmers_bv[j / 8] >> j % 8) & 1) {
        uint64_t kmer_int = canonical_kmer(&t->seq[j], tam_header->k, &f_k, &r_k);
        k = kh_get(kmers_count, h_k, kmer_int);
        if(k != kh_end(h_k)) {
          kmer_count_t *kc = (kmer_count_t*)kh_value(h_k, k);
          for(int n = 0; n < tam_header->k; n++) {
            target_coverage[i][j+n] += kc->count;
          }
        }
      }
    }
  }

  // TODO add mutated k-mer to coverage !!
  tam_file = tam_open(tam_path, "rb");
  tam_header_read(tam_header,tam_file);
  while(tam_record_read(tam_record,tam_file)) {
    k = kh_get(kmers_count, h_k, tam_record->alt_kmers[0]);
    if(k != kh_end(h_k)) {
      kmer_count_t *kc = (kmer_count_t*)kh_value(h_k, k);
      target_coverage[tam_record->target_id][tam_record->pos] += kc->count;
    }
  }
  gzclose(tam_file);

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
  fprintf(stdout, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\n");

  tam_file = tam_open(tam_path, "rb");
  tam_header_read(tam_header,tam_file);
  while(tam_record_read(tam_record,tam_file)) {
    int DP = target_coverage[tam_record->target_id][tam_record->pos];
    int AC = 0;
    k = kh_get(kmers_count, h_k, tam_record->alt_kmers[0]);
    if(k != kh_end(h_k)) {
      kmer_count_t *kc = (kmer_count_t*)kh_value(h_k, k);
      AC = kc->count;
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
    // double p_error = 0.01;
    // long long combi_a = combi(DP, DP - AC);
    // long long combi_b = combi(DP, AC);
    //
    // double p_C_g0 = (double) combi_a * pow(1 - p_error, DP - AC) * pow(p_error, AC);
    // double p_C_g1 = (double) combi_a * 1/pow(2, DP);
    // double p_C_g2 = (double) combi_b * pow(p_error, DP - AC) * pow(1 - p_error, AC);

    double p_g0   = pow((double) (DP - AC) / DP, 2);
    double p_g2   = pow((double) AC / DP, 2);
    double p_g1   = 1 - p_g0 - p_g2;

    // double p_C    = p_C_g0 * p_g0 + p_C_g1 * p_g1 + p_C_g2 * p_g2;

    // double p_g0_C = p_g0 * p_C_g0 / p_C;
    // double p_g1_C = p_g1 * p_C_g1 / p_C;
    // double p_g2_C = p_g2 * p_C_g2 / p_C;
    // double scaling_term = pow(mean_coverage, DP) / factorial(DP) * negative_mean_coverage_exp;
    // double score, p;

    //fprintf(stderr, "p_C: %f, p-G0: %f, p-G1: %f, p-G2: %f\n", p_C, p_g0_C, p_g1_C, p_g2_C);

    // if(p_g1_C > p_g0_C && p_g1_C > p_g2_C) {
    //   strcpy(GT,"0/1");
    //   // p = p_g1_C;
    // } else if(p_g2_C > p_g0_C && p_g2_C > p_g1_C) {
    //   strcpy(GT,"1/1");
    //   // p = p_g2_C;
    // } else {
    //   // p = p_g0_C;
    // }

    if(p_g1 > p_g0 && p_g1 > p_g2) {
      strcpy(GT,"0/1");
    } else if(p_g2 > p_g0 && p_g2 > p_g1) {
      strcpy(GT,"1/1");
    }

    // score = p * scaling_term;
    // FIXME RO is not computed right, because of the DP cumultate count for ref and other alterate alleles...
    int_to_dna(tam_record->alt_kmers[tam_record->n_alt_kmers/2],k_length,kmer);
    char ref_seq[tam_record->ref_seq_l + 1];
    strncpy(ref_seq, &tam_header->target[tam_record->target_id]->seq[tam_record->pos], tam_record->ref_seq_l);
    ref_seq[tam_record->ref_seq_l] = '\0';
    fprintf(stdout, "%s\t%d\t%s\t%s\t%s\t.\tPASS\tDP=%d;AF=%.2f\tGT:RO:AO\t%s:%d:%d\n", tam_header->ref[tam_header->target[tam_record->target_id]->ref_id], tam_header->target[tam_record->target_id]->pos + tam_record->pos + 1, kmer, ref_seq, tam_record->alt_seq, DP, AF, GT, DP - AC, AC);

  }
  gzclose(tam_file);

  tam_header_destroy(tam_header);
  tam_record_destroy(tam_record);
  for(int i = 0; i < tam_header->n_target; i++)
    free(target_coverage[i]);
  free(target_coverage);
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

  fprintf(stdout, "K = %d\n", tam_header->k);
  for(int i = 0; i < tam_header->n_target; i++) {
    tam_target_print(stdout,tam_header,tam_header->target[i]);
  }

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
	else if (strcmp(argv[1], "scan") == 0) tami_scan(argc-1, argv+1);
  else if (strcmp(argv[1], "view") == 0) tami_view(argc-1, argv+1);
	else {
		fprintf(stderr, "[main] unrecognized command '%s'. Abort!\n", argv[1]);
		return 1;
	}
	return 0;
}
