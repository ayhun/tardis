CC=gcc
CFLAGS = -O2 -g -w
LDFLAGS =  -lz -lm -lpthread
SOURCES = multiInd_SetCover.c buffer.c

OBJECTS = $(SOURCES:.c=.o)
EXECUTABLE = set_cover
INSTALLPATH = /usr/local/bin/

all: $(SOURCES) $(EXECUTABLE)
	rm -rf *.o

$(EXECUTABLE): $(OBJECTS) 
	ar -rc libvhsc.a $(OBJECTS)  

.c.o:
	$(CC) -c $(CFLAGS) $< -o $@

clean: 
	rm -f $(EXECUTABLE) *.o *~ *.a