CC=gcc
CFLAGS =  -O2 -g -I htslib -I vhc -I vhsc 
LDFLAGS = htslib/libhts.a vhc/libvhc.a vhsc/libvhsc.a -lz -lm -lpthread
SOURCES = tardis.c cmdline.c common.c processbam.c config.c processfq.c external.c VHtoVCF.c clustering.c
OBJECTS = $(SOURCES:.c=.o)
EXECUTABLE = tardis
INSTALLPATH = /usr/local/bin/

all: $(SOURCES) $(EXECUTABLE)
	rm -rf *.o

$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)

.c.o:
	$(CC) -c $(CFLAGS) $< -o $@

clean: 
	rm -f $(EXECUTABLE) *.o *~

libs: 
	make -C htslib
	make clean -C htslib
	make -C vhc
	make clean -C vhc
	make -C vhsc
	make clean -C vhsc

install:
	cp tardis $(INSTALLPATH)
