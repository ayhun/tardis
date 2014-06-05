CC=gcc
CCP=g++
CFLAGS = -c -O2 -g -I htslib
CPFLAGS= -c -g -O3
LDFLAGS = htslib/libhts.a  -lz -lm -lpthread
SOURCES = tardis.c cmdline.c common.c processbam.c config.c vh_logger.cpp vh_divethandler.cpp vh_hash.cpp vh_createMaxClusterDeletion.cpp vh_heap.cpp vh_gaphandler.cpp vh_repeathandler.cpp vh_createMaxClusterInversion.cpp vh_createMaxClusterInsertion.cpp vh_maximalCluster.cpp
OBJECTS = $(SOURCES:.c=.o) $(SOURCES:.cpp=.o)
EXECUTABLE = tardis
INSTALLPATH = /usr/local/bin/

TARDIS:  tardis.o cmdline.o processbam.o config.o vh_processing.o vh_logger.o vh_divethandler.o vh_hash.o vh_heap.o vh_gaphandler.o vh_maximalCluster.o vh_createMaxClusterInversion.o vh_createMaxClusterDeletion.o vh_createMaxClusterInsertion.o
	$(CC) -g -o TARDIS tardis.o cmdline.o processbam.o config.o vh_processing.o vh_logger.o vh_divethandler.o vh_hash.o vh_heap.o vh_gaphandler.o vh_maximalCluster.o vh_createMaxClusterInversion.o vh_createMaxClusterDeletion.o vh_createMaxClusterInsertion.o
tardis.o: tardis.c tardis.h
	$(CC) $(CFLAGS) tardis.c

cmdline.o: cmdline.c cmdline.h
	$(CC) $(CFLAGS) cmdline.c

common.o: common.c common.h
	$(CC) $(CFLAGS) common.c

processbam.o: processbam.c processbam.h
	$(CC) $(CFLAGS) processbam.c

config.o: config.c config.h
	$(CC) $(CFLAGS) config.c

# vh_commandlineparser.o: vh_commandlineparser.cpp vh_commandlineparser.h
# 	$(CCP) $(CFLAGS) vh_commandlineparser.cpp 

vh_logger.o: vh_logger.cpp vh_logger.h
	$(CCP) $(CFLAGS) vh_logger.cpp 

vh_processing.o:	vh_processing.cpp vh_processing.h
	$(CCP) $(CFLAGS) vh_processing.cpp

vh_divethandler.o: vh_divethandler.cpp vh_divethandler.h
	$(CCP) $(CFLAGS) vh_divethandler.cpp

vh_hash.o: vh_hash.cpp vh_hash.h
	$(CCP) $(CFLAGS) vh_hash.cpp

vh_heap.o: vh_heap.cpp vh_heap.h
	$(CCP) $(CFLAGS) vh_heap.cpp

vh_createMaxClusterDeletion.o: vh_createMaxClusterDeletion.cpp vh_maximalCluster.h
	$(CCP) $(CFLAGS) vh_createMaxClusterDeletion.cpp

vh_createMaxClusterInversion.o: vh_createMaxClusterInversion.cpp vh_maximalCluster.h
	$(CCP) $(CFLAGS) vh_createMaxClusterInversion.cpp

vh_createMaxClusterInsertion.o: vh_createMaxClusterInsertion.cpp vh_maximalCluster.h
	$(CCP) $(CFLAGS) vh_createMaxClusterInsertion.cpp

vh_maximalCluster.o: vh_maximalCluster.cpp vh_maximalCluster.h
	$(CCP) $(CFLAGS) vh_maximalCluster.cpp

vh_gaphandler.o: vh_gaphandler.cpp vh_gaphandler.h
	$(CCP) $(CFLAGS) vh_gaphandler.cpp

vh_repeathandler.o: vh_repeathandler.cpp vh_repeathandler.h
	$(CCP) $(CFLAGS) vh_repeathandler.cpp


libs:
	make clean -C htslib
	make -C htslib

install:
	cp tardis $(INSTALLPATH)


# all: $(SOURCES) $(EXECUTABLE)
# 	rm -rf *.o

# $(EXECUTABLE): $(OBJECTS) 
# 	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)

# .c.o:
# 	$(CC) $(CFLAGS) $< -o $@

# .cpp.o:
# 	$(CCP) $(CPFLAGS) $< -o $@

clean: 
	rm -f $(EXECUTABLE) *.o *~
