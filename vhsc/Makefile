CC=gcc
CFLAGS = -O2 -g -w
LDFLAGS =  -lz -lm -lpthread
SOURCES = vh_setcover.c vh_buffer.c

OBJECTS = $(SOURCES:.c=.o)
EXECUTABLE = set_cover

all: $(SOURCES) $(EXECUTABLE)
	rm -rf *.o

$(EXECUTABLE): $(OBJECTS) 
	ar -rc libvhsc.a $(OBJECTS)  

.c.o:
	$(CC) -c $(CFLAGS) $< -o $@

clean: 
	rm -f $(EXECUTABLE) *.o *~ 
