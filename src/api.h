#ifndef API_H
#define API_H

#define ENSEMBL_API_URL "https://rest.ensembl.org/"

char *api_call(char* query);
char *get_sequence(char *chr, unsigned start, unsigned end);

#endif
