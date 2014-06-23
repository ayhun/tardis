CC=gcc
CFLAGS = -c -O2 -g -I htslib -I vh
LDFLAGS = htslib/libhts.a vh/libvh.a -lz -lm -lpthread
SOURCES = tardis.c cmdline.c common.c processbam.c config.c
OBJECTS = $(SOURCES:.c=.o)
EXECUTABLE = tardis
INSTALLPATH = /usr/local/bin/

all: libs $(SOURCES) $(EXECUTABLE)
	rm -rf *.o

$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)

.c.o:
	$(CC) $(CFLAGS) $< -o $@

clean: 
	rm -f $(EXECUTABLE) *.o *~

libs: htslib readpair

htslib:
	make -C htslib
readpair:
	make -C vh

install:
	cp tardis $(INSTALLPATH)
