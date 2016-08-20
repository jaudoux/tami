CC=gcc
#CFLAGS=-g -Wall -O2 -Wno-unused-function -lm
CFLAGS=-g -Wall -O2 -Wno-unused-function -lm -lpthread
HEADERS=kstring.h kurl.h api.h kthread.h tam.h intervals.h
OBJECTS=$(HEADERS:.h=.o)
LIBS=-lz -lm -lcurl

all:tami

t/tests_tami: tami

test: t/tests_tami
	./t/tests_tami

tami: tami.c $(HEADERS) $(OBJECTS) dna.h
	$(CC) $(CFLAGS) $(OBJECTS) $< -o $@ $(LIBS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@
