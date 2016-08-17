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

#include "api.h"
#include "dna.h"

#define TAMI_VERSION "0.0.3"
#define NB_THREAD_API 10

typedef struct {
  char *chr;
  int start, end;
  char *seq;
} interval_t;

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

typedef struct {
  uint64_t mutk, refk;
  char *chr;
  int pos;
  char mutn, refn;
} kmer_mut_t;

int write_kmer_mut(kmer_mut_t *kmut, FILE *fp) {
  size_t chr_l = strlen(kmut->chr);
  fwrite(&kmut->mutk,   sizeof(kmut->mutk), 1, fp);
  fwrite(&kmut->refk,   sizeof(kmut->refk), 1, fp);
  fwrite(&chr_l,        sizeof(size_t),     1, fp);
  fwrite(kmut->chr,     sizeof(char),   chr_l, fp);
  fwrite(&kmut->pos,    sizeof(kmut->pos),  1, fp);
  fwrite(&kmut->mutn,   sizeof(kmut->mutn), 1, fp);
  return fwrite(&kmut->refn,   sizeof(kmut->refn), 1, fp);
}

int read_kmer_mut(kmer_mut_t *kmut, FILE *fp) {
  //fread(&kmut->mute_kmer,    sizeof(uint64_t),     1, kmer_file);
  size_t chr_l, r;
  r = fread(&kmut->mutk,   sizeof(kmut->mutk), 1, fp);
  if(r == 1) r = fread(&kmut->refk,   sizeof(kmut->refk), 1, fp);
  else return 0;
  if(r == 1) r = fread(&chr_l,        sizeof(size_t),     1, fp);
  else return 0;
  if(r == 1) {
    if(kmut->chr)
      kmut->chr = realloc(kmut->chr, sizeof(char) * (chr_l + 1));
    else
      kmut->chr = malloc(sizeof(char) * (chr_l + 1));
    kmut->chr[chr_l] = '\0';
    r = fread(kmut->chr,     sizeof(char),   chr_l, fp);
  } else {
    return 0;
  }
  if(r == chr_l) r = fread(&kmut->pos,    sizeof(kmut->pos),  1, fp);
  else return 0;
  if(r == 1) r = fread(&kmut->mutn,   sizeof(kmut->mutn), 1, fp);
  else return 0;
  if(r == 1) r = fread(&kmut->refn,   sizeof(kmut->refn), 1, fp);
  else return 0;
  if(r == 1) return 1;
  else return 0;
}

typedef struct {
  uint64_t altk, refk;
} kmer_alt_t;

int write_kmer_alt(kmer_alt_t *kalt, FILE *fp) {
  fwrite(&kalt->altk,   sizeof(kalt->altk), 1, fp);
  return fwrite(&kalt->refk,   sizeof(kalt->refk), 1, fp);
}

int read_kmer_alt(kmer_alt_t *kalt, FILE *fp) {
  size_t r;
  r = fread(&kalt->altk,   sizeof(kalt->altk), 1, fp);
  if(r == 1) r = fread(&kalt->refk,   sizeof(kalt->refk), 1, fp);
  else return 0;
  if(r == 1) return 1;
  else return 0;
  //return fread(&kalt->refk,   sizeof(kalt->refk), 1, fp);
}

static void get_sequence_iter(void *_int_a, long i, int tid) {
  interval_t ** int_a   = (interval_t **)_int_a;
  interval_t *interval = int_a[i];
  interval->seq = get_sequence(interval->chr, interval->start, interval->end);
}

int main(int argc, char *argv[]) {

  // Init the array that stores intervals
  kvec_t(interval_t*) interval_array;
  kv_init(interval_array);

  // Init the hash wher k-mer counts are stored
  khiter_t k, k2;
  khash_t(kmers) *h = kh_init(kmers);


  char *bed_file = NULL, *reference_fasta = NULL, *fastq_file;
  int use_derived_kmers = 0;
  float min_alternate_fraction = 0.2;
  int min_alternate_count = 3;
  int min_coverage = 10;
  int max_coverage = -1;
  int k_length      = 30;
  int k_middle      = k_length / 2;
  int debug = 0;
  int version = 0;


  int c;
  while ((c = getopt(argc, argv, "dsvF:C:m:M:d:k:r:")) >= 0) {
		switch (c) {
			case 'F': min_alternate_fraction = atof(optarg); break;
      case 'C': min_alternate_count = atoi(optarg); break;
      case 'm': min_coverage = atoi(optarg); break;
      case 'M': max_coverage = atoi(optarg); break;
      case 'k': k_length = atoi(optarg); break;
      case 's': use_derived_kmers = 1; break;
      case 'r': reference_fasta = optarg; break;
      case 'v': version = 1; break;
      case 'd': debug = 1; break;
		}
	}

  if(version) {
    fprintf(stderr, "tami v%s\n",TAMI_VERSION);
    return 1;
  } else if (optind + 1 >= argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   tami [options] <in.bed> <in.fq>\n\n");
		fprintf(stderr, "Options: -F FLOAT  min alternate allele frequency [%.2f]\n", min_alternate_fraction);
    fprintf(stderr, "         -C INT    min alternate allele count (min_value: 1)[%d]\n", min_alternate_count);
    fprintf(stderr, "         -m INT    min coverage [%d]\n", min_coverage);
    fprintf(stderr, "         -M INT    max coverage (min_value: 1)[%d]\n", max_coverage);
    fprintf(stderr, "         -k INT    length of k-mers (max_value: 32) [%d]\n", k_length);
    fprintf(stderr, "         -r STR    reference FASTA (otherwise sequences are retrieved from rest.ensembl.org)\n");
    fprintf(stderr, "         -s        sensitive mode (able 2 mutation per-kmer)\n");
    fprintf(stderr, "                   Warning : important memory overload and possible loss of accuracy,\n");
    fprintf(stderr, "                             this option is only adviced for amplicon sequencing.\n");
    fprintf(stderr, "         -v        print tami version number\n");
    fprintf(stderr, "         -d        debug mode\n");
		fprintf(stderr, "\n");
		return 1;
	}

  bed_file = argv[optind++];

  // verify Options
  if(min_coverage < 1) {
    fprintf(stderr, "Invalid value (%d) for min_coverage (-m option).\n", min_coverage);
    return 1;
  }
  if(min_alternate_count < 1) {
    fprintf(stderr, "Invalid value (%d) for min_alternate_count (-C option).\n", min_coverage);
    return 1;
  }
  if(k_length > 32 || k_length < 1) {
    fprintf(stderr, "Invalid value (%d) for k_length (-k option).\n", k_length);
    return 1;
  }

  if(debug) {
    fprintf(stderr, "Bed file is: %s\n",bed_file);
    fprintf(stderr, "min_alternate_count (C): %d\n",min_alternate_count);
    fprintf(stderr, "min_alternate_fraction (F): %f\n",min_alternate_fraction);
  }

  /* READ THE BED FILE and load intervals*/
  fprintf(stderr, "Reading BED file...\n");

  gzFile fp;
	kstream_t *ks;
  int dret;
	kstring_t *str;
  interval_t *interval;

  kstring_t *chr;
  int start, end;

  str = calloc(1, sizeof(kstring_t));
  chr = calloc(1, sizeof(kstring_t));
  fp = strcmp(bed_file, "-")? gzopen(bed_file, "r") : gzdopen(fileno(stdin), "r");
  ks = ks_init(fp);

  while (ks_getuntil(ks, 0, str, &dret) >= 0) {
    // remove chr prefix if any
    if(str->l > 3 && strncmp("chr", str->s, 3) == 0) {
      kputs(&str->s[3],chr);
    } else {
      kputs(str->s,chr);
    }

    if(dret != '\n') {
      if(ks_getuntil(ks, 0, str, &dret) > 0 && isdigit(str->s[0])) {
        start = atoi(str->s);
        if(dret != '\n') {
          if(ks_getuntil(ks, 0, str, &dret) > 0 && isdigit(str->s[0])) {
            end = atoi(str->s);
            // create a new interval
            interval = interval_init(ks_release(chr),start + 1, end);
            kv_push(interval_t*,interval_array,interval);
            if(debug)
              fprintf(stderr, "Read interval : %s:%d..%d\n", interval->chr,interval->start,interval->end);
          }
        }
      }
    }
    // skip the rest of the line
		if (dret != '\n') while ((dret = ks_getc(ks)) > 0 && dret != '\n');
    // release the current chr
    if(chr->l > 0) { chr->l = 0; };
  }
  ks_destroy(ks);
  gzclose(fp);
  free(str->s); free(str);


  /* SORT AND REMOVE OVERLAPPING INTERVALS */

  // Sort and remove overlapping intervals
  //fprintf(stderr, "Number of intervals: %d\n",kv_size(interval_array));
  qsort(interval_array.a, kv_size(interval_array), sizeof(interval_array.a[0]), cmp_interval);
  int j = 0;
  for(int i = 1; i < kv_size(interval_array); i++) {
    //interval_t *curr_interval = kv_A(interval_array,i);
    if(strcmp(interval_array.a[j]->chr,interval_array.a[i]->chr) == 0 &&
      interval_array.a[i]->start <= interval_array.a[j]->end) {
      if(interval_array.a[i]->end > interval_array.a[j]->end) {
        interval_array.a[j]->end = interval_array.a[i]->end;
      }
    } else {
      interval_array.a[++j] = interval_array.a[i];
    }
  }
  // Set the new size of the interval array
  interval_array.n = j + 1;
  //ks_mergesort(interval, kv_size(interval_array), *interval_array.a, 0);

  /* LOAD SEQUENCES FROM FASTA OR ENSEMBL API (if not FASTA provided) */
  if(reference_fasta) {
    fprintf(stderr, "Load sequences from %s...\n", reference_fasta);
    int l;
    fp = gzopen(reference_fasta, "r");
    kseq_t *seq = kseq_init(fp);
    while ((l = kseq_read(seq)) >= 0) {
      //k = kh_get(chr_intervals, h, seq->name.s);
      for(int i = 0; i < kv_size(interval_array); i++) {
        interval = kv_A(interval_array,i);
        if(!interval->seq && strcmp(interval->chr, seq->name.s) == 0) {
          int int_length = interval_length(interval);
          char *interval_seq  = (char *)malloc(int_length * sizeof(char) + 1);
          memcpy(interval_seq, &seq->seq.s[interval->start - 1], int_length);
          interval_seq[int_length] = '\0';
          interval->seq = interval_seq;
        }
      }
    }
    //kt_for(NB_THREAD_API, get_sequence_iter, interval_array.a, kv_size(interval_array));
  } else {
    fprintf(stderr, "Load sequences from Ensembl REST API...\n");
    kt_for(NB_THREAD_API, get_sequence_iter, interval_array.a, kv_size(interval_array));
  }

  /* CREATE THE MUTATED K-MER HASH */
  fprintf(stderr, "Create mutated k-mer hash...\n");
  int ret, is_missing;

  int ref_kmer_pos, mut_pos, pos;
  char ref_nuc;
  uint64_t ref_kmer, mut_kmer;
  char kmer[k_length], kmer2[k_length];;
  FILE *kmer_file = fopen("kmers.txt", "wb+"); // TODO handle error of opening

  // Loop over the interval array and print sequences
  for(int i = 0; i < kv_size(interval_array); i++) {
    interval = kv_A(interval_array,i);

    if(!interval->seq) {
      fprintf(stderr, "No sequence found for interval %s:%d-%d\n", interval->chr, interval->start, interval->end);
      continue;
    }

    //char *seq = get_sequence(interval->chr, interval->start, interval->end);
    char *seq = interval->seq;
    if(debug)
      fprintf(stderr, ">%s:%d..%d\n%s\n", interval->chr, interval->start, interval->end,seq);

    size_t seq_length = strlen(seq);

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
      ref_nuc = seq[mut_pos + ref_kmer_pos];
      ref_kmer = dna_to_int(&seq[ref_kmer_pos],k_length,1);
      memmove(kmer, &seq[ref_kmer_pos], sizeof kmer);

      // if(pos == 27581443) {
      //   fprintf(stderr, "ref_seq  : %s\n", seq);
      //   fprintf(stderr, "ref_kmer : %s\n", kmer);
      //   fprintf(stderr, "ref_nucl : %c\n", ref_nuc);
      // }


      for(int p = 0; p < NB_NUCLEOTIDES; p++) {
        if(NUCLEOTIDES[p] != ref_nuc) {
          kmer[mut_pos] = NUCLEOTIDES[p];
          // if(pos == 27581443) {
          //   fprintf(stderr, "mut_kmer %c => %c (%d): %s\n", ref_nuc, NUCLEOTIDES[p], mut_pos, kmer);
          // }
          //fprintf(stderr, "mut_kmer :       %s\n", kmer);
          mut_kmer = dna_to_int(kmer,k_length,1);

          // if(pos == 27581443) {
          //   int_to_dna(mut_kmer,k_length,kmer2);
          //   fprintf(stderr, "canonical_kmer : %s\n", kmer2);
          // }

          if(kh_get(kmers, h, mut_kmer) == kh_end(h)) {
            k = kh_put(kmers, h, mut_kmer, &ret);
          	kh_value(h, k) = 0;


            kmer_mut_t kmer_mut_struct = { mut_kmer, ref_kmer, interval->chr, pos, NUCLEOTIDES[p], ref_nuc };
            write_kmer_mut(&kmer_mut_struct, kmer_file);
            //fprintf(kmer_file, "%" PRIu64 "\t%" PRIu64 "\t%s\t%d\t%c\t%c\n",mut_kmer,ref_kmer,interval->chr,(interval->start+ref_kmer_pos+mut_pos),ref_nuc,NUCLEOTIDES[p]);
          }
        } else {
          if(kh_get(kmers, h, ref_kmer) == kh_end(h)) {
            k = kh_put(kmers, h, ref_kmer, &ret);
          	kh_value(h, k) = 1;
          }
        }
      }
    }
    interval_destroy(interval);
  }
  fclose(kmer_file);

  if(use_derived_kmers) {
    // We update the size of the hash to fit the derived kmers. Theoritically, it should
    // kh_size * k * 3, but we in practice it will be closer to a factor 2.
    kh_resize(kmers, h, kh_size(h)*k_length*2);
    fprintf(stderr, "Create 2-nuc mutated k-mer (-s option)...\n");
    // Add all derived k-mer with one more mutation
    FILE *kmer_alt_file = fopen("kmers_alt.txt", "wb+"); // TODO handle error of opening
    for(k = kh_begin(h); k != kh_end(h); ++k) {
      if (!kh_exist(h,k)) continue;
      if(kh_val(h,k) == 0) {
        ref_kmer = (uint64_t)kh_key(h,k);
        for (int n = 0; n < k_length; n++) {
          ref_nuc = nuc_from_int_dna(ref_kmer,k_length,n);
          for(int p = 0; p < NB_NUCLEOTIDES; p++) {
            if(ref_nuc != NUCLEOTIDES[p]) {
              mut_kmer = mut_int_dna(ref_kmer, k_length, n,  NUCLEOTIDES[p]);
              if(kh_get(kmers, h, mut_kmer) == kh_end(h)) {
                k2 = kh_put(kmers, h, mut_kmer, &ret);
                kh_value(h, k2) = 1;
                kmer_alt_t kmer_alt = { mut_kmer, ref_kmer };
                write_kmer_alt(&kmer_alt, kmer_alt_file);
              }
            }
          }
        }
      }
    }
    fclose(kmer_alt_file);
  }

  for(k = kh_begin(h); k != kh_end(h); ++k) {
    if (!kh_exist(h,k)) continue;
    if(kh_val(h,k) == 1) {
      kh_value(h, k) = 0;
    }
  }

  /* READ THE FASTQ FILES */
  while(optind < argc) {
    fastq_file = argv[optind++];
    fprintf(stderr, "Reading %s FASTQ file...\n", fastq_file);

    kseq_t *seq;
    int l;
    uint64_t kmer_int;
    int nb_reads = 0;
    fp = gzopen(fastq_file, "r");
    seq = kseq_init(fp);
    while ((l = kseq_read(seq)) >= 0) {
      nb_reads++;
      if(strlen(seq->seq.s) >= k_length) {
        for(int i = 0; i < strlen(seq->seq.s) - k_length + 1; i++) {
          kmer_int = dna_to_int(&seq->seq.s[i], k_length, 1);
          k = kh_get(kmers, h, kmer_int);
          if(k != kh_end(h)) {
            kh_value(h, k) = kh_value(h, k) + 1;
          }
        }
      }
      if (nb_reads % 50000 == 0){
        fprintf(stderr, "*");
      }
    }
    fprintf(stderr, "\n%d reads parsed.\n",nb_reads);
    kseq_destroy(seq);
    gzclose(fp);
  }

  if(use_derived_kmers) {
    fprintf(stderr, "Update counts with alternate kmers...\n");
    // Update the mutated k-mer count with derived k-mer counts
    FILE *kmer_alt_file= fopen("kmers_alt.txt", "r");
    kmer_alt_t kmer_alt_struct;
    while(!feof(kmer_alt_file)) {
      if(read_kmer_alt(&kmer_alt_struct, kmer_alt_file)) {
        k   = kh_get(kmers, h, kmer_alt_struct.refk);
        k2  = kh_get(kmers, h, kmer_alt_struct.altk);
        if(k != kh_end(h) && k2 != kh_end(h)) {
          kh_value(h, k) = kh_value(h, k) + kh_value(h, k2);
        }
        //fprintf(stderr, "%s:%d\t%c => %c\n", kmer_mut_struct.chr,kmer_mut_struct.pos,kmer_mut_struct.refn,kmer_mut_struct.mutn);
      }
    }
  }

  // Update reference k-mer count with mutated k-mer counts
  kmer_file = fopen("kmers.txt", "r");
  kmer_mut_t kmer_mut_struct = { 0, 0, NULL, 0, '\0', '\0' };
  fprintf(stderr, "Update counts for reference kmers...\n");
  while(!feof(kmer_file)) {
    if(read_kmer_mut(&kmer_mut_struct, kmer_file)) {
      k   = kh_get(kmers, h, kmer_mut_struct.refk);
      k2  = kh_get(kmers, h, kmer_mut_struct.mutk);
      if(k != kh_end(h) && k2 != kh_end(h)) {
        kh_value(h, k) = kh_value(h, k) + kh_value(h, k2);
      }
      //fprintf(stderr, "%s:%d\t%c => %c\n", kmer_mut_struct.chr,kmer_mut_struct.pos,kmer_mut_struct.refn,kmer_mut_struct.mutn);
    }
  }
  fclose(kmer_file);

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

  kmer_file = fopen("kmers.txt", "r");
  while(!feof(kmer_file)) {
    if(read_kmer_mut(&kmer_mut_struct, kmer_file)) {
      k   = kh_get(kmers, h, kmer_mut_struct.refk);
      k2  = kh_get(kmers, h, kmer_mut_struct.mutk);
      if(k != kh_end(h) && k2 != kh_end(h)) {
        int AC = kh_val(h,k2);
        int DP = kh_val(h,k);
        // if(kmer_mut_struct.pos == 27581443) {
        //   fprintf(stderr, "%s\t%d\t%s\t%c\t%c\tDP=%d;AF=%.2f;AC=%d\n", kmer_mut_struct.chr, kmer_mut_struct.pos, kmer, kmer_mut_struct.refn, kmer_mut_struct.mutn,DP,AF,AC);
        // }
        if(AC < min_alternate_count)
          continue;
        if(DP < min_coverage)
          continue;
        float AF = (float) AC / DP;
        if(AF < min_alternate_fraction)
          continue;
        if(max_coverage > 0 && DP > max_coverage)
          continue;
        int_to_dna(kmer_mut_struct.mutk,k_length,kmer);
        fprintf(stdout, "%s\t%d\t%s\t%c\t%c\t.\t.\tDP=%d;AF=%.2f;AC=%d\n", kmer_mut_struct.chr, kmer_mut_struct.pos, kmer, kmer_mut_struct.refn, kmer_mut_struct.mutn,DP,AF,AC);
      }
      //fprintf(stderr, "%s:%d\t%c => %c\n", kmer_mut_struct.chr,kmer_mut_struct.pos,kmer_mut_struct.refn,kmer_mut_struct.mutn);
    }
  }
  fclose(kmer_file);

  free(kmer_mut_struct.chr);
  kv_destroy(interval_array);
	return 0;
}
