TARDIS_VERSION := "0.1-alpha"
TARDIS_UPDATE := "April 28, 2015"
TARDIS_DEBUG := 1
BUILD_DATE := "$(shell date)"
CC=gcc
CFLAGS =  -O2 -g -I htslib -I vhc -I vhsc -DTARDIS_VERSION=\"$(TARDIS_VERSION)\" -DBUILD_DATE=\"$(BUILD_DATE)\" -DTARDIS_UPDATE=\"$(TARDIS_UPDATE)\" -DTARDIS_DEBUG=$(TARDIS_DEBUG)
LDFLAGS = htslib/libhts.a vhc/libvhc.a vhsc/libvhsc.a -lz -lm -lpthread
SOURCES = tardis.c cmdline.c common.c processbam.c config.c processfq.c external.c vh.c
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
	make -C vhc
	make -C vhsc

install:
	cp tardis $(INSTALLPATH)
