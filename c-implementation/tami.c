#include <stdio.h>  //fprintf
#include <stdlib.h> //free
#include <zlib.h>
#include <inttypes.h>

#include "kstring.h"
#include "ksort.h"
#include "khash.h"
#include "kvec.h"
#include "kthread.h"
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)
KHASH_MAP_INIT_INT64(kmers, uint32_t)
KHASH_MAP_INIT_STR(references, int32_t)

#include "api.h"
#include "dna.h"
#include "tam.h"
#include "intervals.h"

#define TAMI_VERSION "0.1.0"
#define NB_THREAD_API 10

static void get_sequence_iter(void *_int_a, long i, int tid) {
  interval_t ** int_a   = (interval_t **)_int_a;
  interval_t *interval = int_a[i];
  interval->seq = get_sequence(interval->chr, interval->start, interval->end);
}

int tami_build(int argc, char *argv[]) {
  char *bed_file = NULL, *reference_fasta = NULL, *output_file = "kmers.tam";
  int use_derived_kmers = 0;
  int k_length      = 30;
  int k_middle      = k_length / 2;
  int debug = 0;

  int c;
  while ((c = getopt(argc, argv, "dsk:r:o:")) >= 0) {
		switch (c) {
      case 'k': k_length = atoi(optarg); break;
      case 'r': reference_fasta = optarg; break;
      case 's': use_derived_kmers = 1; break;
      case 'o': output_file = optarg; break;
      case 'd': debug = 1; break;
		}
	}

  if (optind == argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   tami build [options] <in.bed>\n\n");
		fprintf(stderr, "Options: -k INT    length of k-mers (max_value: 32) [%d]\n", k_length);
    fprintf(stderr, "         -r STR    reference FASTA (otherwise sequences are retrieved from rest.ensembl.org)\n");
    fprintf(stderr, "         -s        sensitive mode (able 2 mutation per-kmer)\n");
    fprintf(stderr, "                   Warning : important memory overload and possible loss of accuracy,\n");
    fprintf(stderr, "                             this option is only adviced for amplicon sequencing.\n");
    fprintf(stderr, "         -o STR    output file (default : 'kmers.tam')\n");
    fprintf(stderr, "         -v        print tami version number\n");
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

  gzFile fp;
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

  kv_init(reference_names);

  tmp_output_file= malloc(strlen(output_file) + strlen(".tmp") + 1);
  tmp_output_file[0] = '\0';
  strcat(tmp_output_file,output_file);
  strcat(tmp_output_file,".tmp");

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
  *                      REATE THE MUTATED K-MER HASH
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
  if(!tam_record->ref_seq)   tam_record->ref_seq     = malloc(sizeof(char));
  if(!tam_record->alt_seq)   tam_record->alt_seq     = malloc(sizeof(char));
  tam_record->n_ref_kmers = 1;
  tam_record->n_alt_kmers = 1;
  tam_record->ref_seq_l   = 1;
  tam_record->alt_seq_l   = 1;

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

      //fprintf(stderr, "seq_length: %d, ref_kmer_pos: %d, mut_pos: %d\n", seq_length, ref_kmer_pos, mut_pos);
      pos = interval->start + ref_kmer_pos + mut_pos;
      ref_nuc = interval->seq[mut_pos + ref_kmer_pos];

      if(ref_kmer_pos > prev_ref_kmer_pos) {
        ref_kmer = next_canonical_kmer(k_length, interval->seq[ref_kmer_pos + k_length - 1], &forward_ref_kmer, &reverse_ref_kmer);
        prev_ref_kmer_pos = ref_kmer_pos;
      }

      tam_record->ref_id        = ref_id;
      tam_record->pos           = pos;
      tam_record->ref_seq[0]    = ref_nuc;
      tam_record->ref_kmers[0]  = ref_kmer;

      for(int p = 0; p < NB_NUCLEOTIDES; p++) {
        if(NUCLEOTIDES[p] != ref_nuc) {
          mut_kmer = mut_int_dna(forward_ref_kmer, k_length, mut_pos, NUCLEOTIDES[p]);
          reverse_mut_kmer = int_revcomp(mut_kmer, k_length);
          if(reverse_mut_kmer < mut_kmer)
            mut_kmer = reverse_mut_kmer;

          if(kh_get(kmers, h_k, mut_kmer) == kh_end(h_k)) {
            k = kh_put(kmers, h_k, mut_kmer, &ret);
          	kh_value(h_k, k) = 0;
            tam_record->alt_seq[0]    = NUCLEOTIDES[p];
            tam_record->alt_kmers[0]  = mut_kmer;
            tam_record_write(tam_record, tam_file);
          }
        } else {
          if(kh_get(kmers, h_k, ref_kmer) == kh_end(h_k)) {
            k = kh_put(kmers, h_k, ref_kmer, &ret);
          	kh_value(h_k, k) = 2;
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
          tam_record->alt_kmers = realloc(tam_record->alt_kmers, sizeof(uint64_t) * tam_header->k * (NB_NUCLEOTIDES - 1) + 1);

          for (int n = 0; n < tam_header->k; n++) {
            ref_nuc = nuc_from_int_dna(tam_record->alt_kmers[0], tam_header->k, n);
            for(int p = 0; p < NB_NUCLEOTIDES; p++) {
              if(ref_nuc != NUCLEOTIDES[p]) {
                mut_kmer = mut_int_dna(tam_record->alt_kmers[0], tam_header->k, n,  NUCLEOTIDES[p]);
                if(kh_get(kmers, h_k, mut_kmer) == kh_end(h_k)) {
                  k2 = kh_put(kmers, h_k, mut_kmer, &ret);
                  kh_value(h_k, k2) = 1;
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
    int l, nb_removed_kmers = 0;
    uint64_t kmer_int, forward_kmer_int = 0, reverse_kmer_int = 0;
    fp = gzopen(reference_fasta, "r");
    kseq_t *seq = kseq_init(fp);
    while ((l = kseq_read(seq)) >= 0) {
      if(debug)
        fprintf(stderr, "Checking chr %s\n", seq->name.s);
      if(seq->seq.l >= k_length) {
        //canonical_kmer(const char *dna, size_t dna_length, uint64_t &forward_kmer, uint64_t &reverse_kmer)
        for(int i = 0; i < seq->seq.l - k_length + 1; i++) {
          // We could do something faster here by upding the previous kmer_int without
          // One new nucleotides. This should be k-time faster (in theory)
          if(i > 0) {
            kmer_int = next_canonical_kmer(k_length, seq->seq.s[i + k_length - 1], &forward_kmer_int, &reverse_kmer_int);
          } else {
            kmer_int = canonical_kmer(seq->seq.s, k_length, &forward_kmer_int, &reverse_kmer_int);
          }
          //kmer_int = dna_to_int(&seq->seq.s[i], k_length, 1);
          k = kh_get(kmers, h_k, kmer_int);
          if(k != kh_end(h_k) && kh_value(h_k, k) != 2) {
            kh_del(kmers, h_k, k);
            nb_removed_kmers++;
          }
        }
      }
    }
    fprintf(stderr, "%d mutated kmers have been removed\n", nb_removed_kmers);

    fprintf(stderr, "Updating the TAM file...\n");
    // Write the final tam file
    tam_file = tam_open(output_file, "rb");
    gzFile tam_file_clean = tam_open(tmp_output_file, "wb");

    if(tam_header_read(tam_header,tam_file)) {
      tam_header->n_kmers = kh_size(h_k);
      tam_header_write(tam_header, tam_file_clean);
      while(tam_record_read(tam_record,tam_file)) {
        // Verify alt_kmers
        int j = 0;
        for(int i = 0; i < tam_record->n_alt_kmers; i++) {
          k = kh_get(kmers, h_k, tam_record->alt_kmers[i]);
          if(k != kh_end(h_k))
            tam_record->alt_kmers[j++] = tam_record->alt_kmers[i];
        }
        tam_record->n_alt_kmers = j;
        if(j == 0) continue;
        tam_record_write(tam_record, tam_file_clean);
      }
    }
    gzclose(tam_file);
    gzclose(tam_file_clean);
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
  int ret, k_length;
  gzFile fp;

  tam_header_read(tam_header,tam_file);
  kmer = malloc(sizeof(char) * (tam_header->k + 1));
  k_length = tam_header->k;
  kh_resize(kmers, h_k, tam_header->n_kmers);
  while(tam_record_read(tam_record,tam_file)) {
    //fprintf(stderr, "TOTO\n");
    for(int i = 0; i < tam_record->n_ref_kmers; i++) {
      k = kh_put(kmers, h_k, tam_record->ref_kmers[i], &ret);
      kh_value(h_k, k) = 0;
    }
    for(int i = 0; i < tam_record->n_alt_kmers; i++) {
      k = kh_put(kmers, h_k, tam_record->alt_kmers[i], &ret);
      kh_value(h_k, k) = 0;
    }
  }
  fprintf(stderr, "%d k-mers loaded into memory\n", (int)kh_size(h_k));
  gzclose(tam_file);

  /*****************************************************************************
  *                             SCAN FASTQ FILES
  *****************************************************************************/

  // Set all counters to 0
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

  /*****************************************************************************
  *                             PRINT OUTPUT
  *****************************************************************************/

  /* FILTER AND PRINT OUTPUT VCF */
  fprintf(stdout, "##fileformat=VCFv4.1\n");
  fprintf(stdout, "##source=TaMI v%s\n", TAMI_VERSION);
  fprintf(stdout, "##commandline=");
  for (int i = 0; i < argc; i++)
    fprintf(stdout, " %s", argv[i]);
  fprintf(stdout, "\n##INFO=<ID=DP,Number=1,Type=Integer,");
  fprintf(stdout, "Description=\"Total read depth at the locus\">\n");
  fprintf(stdout, "##INFO=<ID=AC,Number=A,Type=Integer,");
  fprintf(stdout, "Description=\"Total number of alternate alleles in called genotypes\">\n");
  fprintf(stdout, "##INFO=<ID=AF,Number=A,Type=Float,");
  fprintf(stdout, "Description=\"Estimated allele frequency in the range (0,1]\">\n");
  fprintf(stdout, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");

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
    int_to_dna(max_alt_kmer,k_length,kmer);
    fprintf(stdout, "%s\t%d\t%s\t%s\t%s\t.\t.\tDP=%d;AF=%.2f;AC=%d\n", tam_header->ref[tam_record->ref_id], tam_record->pos, kmer, tam_record->ref_seq, tam_record->alt_seq, DP, AF, AC);

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
	fprintf(stderr, "Command: build     Create a mutated k-mer lib (TAM file)\n");
	fprintf(stderr, "         scan      Scan a FASTQ files agains a TAM file\n");
	fprintf(stderr, "\n");
	return 1;
}

int main(int argc, char *argv[])
{
	if (argc == 1) return usage();
	if (strcmp(argv[1], "build") == 0) tami_build(argc-1, argv+1);
	else if (strcmp(argv[1], "scan") == 0) tami_scan(argc-1, argv+1);
	else {
		fprintf(stderr, "[main] unrecognized command '%s'. Abort!\n", argv[1]);
		return 1;
	}
	return 0;
}
