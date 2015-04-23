CC=gcc
CFLAGS = -O2 -g -w -I ../htslib
LDFLAGS =  -lz -lm -lpthread ../htslib/libhts.a
SOURCES = vh_setcover.c vh_buffer.c

OBJECTS = $(SOURCES:.c=.o)
EXECUTABLE = libvhsc.a

all: $(SOURCES) $(EXECUTABLE)
	rm -rf *.o

$(EXECUTABLE): $(OBJECTS) 
	ar -rc libvhsc.a $(OBJECTS)  

.c.o:
	$(CC) -c $(CFLAGS) $< -o $@

clean: 
	rm -f $(EXECUTABLE) *.o *~ 
