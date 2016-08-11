#include <stdlib.h>
#include <stdio.h>

#include "minunit.h"
#include "string.h"
#include "../api.h"
#include "../dna.h"

MU_TEST(test_get_sequence) {
  char* seq = get_sequence("X",1000000,1000100);
  mu_check(strcmp(seq,"CTGTAGAAACATTAGCCTGGCTAACAAGGTGAAACCCCATCTCTACTAACAATACAAAATATTGGTTGGGCGTGGTGGCGGGTGCTTGTAATCCCAGCTAC") == 0);
  free(seq);
}

MU_TEST(test_dna_to_int) {
  char DNA[] = "ATGCGTAGCC";
  char DNA_BIS[10];
  int_to_dna(dna_to_int(DNA,strlen(DNA),0),10,&DNA_BIS);
  mu_check(strcmp(DNA,DNA_BIS) == 0);
}

MU_TEST_SUITE(test_suite) {
    MU_RUN_TEST(test_get_sequence);
    MU_RUN_TEST(test_dna_to_int);
}

int main(int argc, char *argv[]) {
    MU_RUN_SUITE(test_suite);
    MU_REPORT();
    return 0;
}
