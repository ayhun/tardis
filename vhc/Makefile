CC=gcc
CFLAGS = -O2 -g -I ../htslib
LDFLAGS = ../htslib/libhts.a  -lz -lm -lpthread
SOURCES = vh_main.c vh_logger.c vh_divethandler.c vh_hash.c vh_createMaxClusterDeletion.c vh_heap.c vh_gaphandler.c vh_repeathandler.c vh_createMaxClusterInversion.c vh_createMaxClusterInsertion.c vh_maximalCluster.c
	
OBJECTS = $(SOURCES:.c=.o)
EXECUTABLE = libvhc.a
INSTALLPATH = /usr/local/bin/

all: $(SOURCES) $(EXECUTABLE)
	rm -rf *.o

$(EXECUTABLE): $(OBJECTS) 
	ar -rc libvhc.a  $(OBJECTS) 

.c.o:
	$(CC) -c $(CFLAGS) $< -o $@

clean: 
	rm -f $(EXECUTABLE) *.o *~ 
