CC=gcc
CFLAGS = -c -O2 -g
LDFLAGS = -lz -lm
SOURCES = tardis.c cmdline.c common.c processbam.c
OBJECTS = $(SOURCES:.c=.o)
EXECUTABLE = tardis

all: $(SOURCES) $(EXECUTABLE)
	rm -rf *.o
		
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)

.c.o:
	$(CC) $(CFLAGS) $< -o $@

clean: 
	rm -f $(EXECUTABLE) *.o *~

	
