CC=gcc
#CFLAGS=-g -Wall -O2 -Wno-unused-function -lm
CFLAGS=-g -Wall -Wno-unused-function -lm
HEADERS=kstring.h kurl.h api.h
OBJECTS=$(HEADERS:.h=.o)
LIBS=-lz -lm -lcurl

all:tami

t/tests_tami: tami

test: t/tests_tami
	./t/tests_tami

%: %.c $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) $< -o $@ $(LIBS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@
