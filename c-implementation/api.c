#include <string.h> // memset()
#include <unistd.h> // getopt()
#include <stdlib.h>

#include "kurl.h"
#include "kstring.h"
#include "api.h"

char *api_call(char* query) {
  kurl_t *f;
	int l, l_buf = 0x10000;
	off_t rest = -1;
	char *buf;
  kurl_opt_t opt;
  kstring_t url = { 0, 0, NULL };
  kstring_t res = { 0, 0, NULL };

  kputs(ENSEMBL_API_URL,&url);
  kputs(query,&url);

	kurl_init();
	f = kurl_open(url.s, &opt);
	if (f == 0) {
		fprintf(stderr, "ERROR: fail to open URL\n");
		return "";
	}

	buf = (char*)calloc(l_buf, 1);
	while (rest != 0) {
		int to_read = rest > 0 && rest < l_buf? rest : l_buf;
		l = kurl_read(f, buf, to_read);
		if (l == 0) break;
    kputs(buf,&res);
		//fwrite(buf, 1, l, stdout);
		rest -= l;
	}

	free(buf);
  free(url.s);
	kurl_close(f);
	kurl_destroy();
  return res.s;
}

char *get_sequence(char *chr, unsigned start, unsigned end) {
  kstring_t query = { 0, 0, NULL };
  ksprintf(&query,"/sequence/region/human/%s:%d..%d:1?content-type=text/plain",chr,start,end);
  char* res = api_call(query.s);
  free(query.s);
  return res;
}
