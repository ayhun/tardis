CC=gcc
CFLAGS= -c -g -O3

VH : vh_main.o vh_commandlineparser.o vh_logger.o vh_divethandler.o vh_hash.o vh_createMaxClusterDeletion.o vh_heap.o vh_gaphandler.o vh_repeathandler.o vh_createMaxClusterInversion.o vh_createMaxClusterInsertion.o vh_maximalCluster.o
	ar -rc libvh.a vh_main.o vh_commandlineparser.o vh_logger.o vh_divethandler.o vh_hash.o vh_createMaxClusterDeletion.o vh_heap.o vh_gaphandler.o vh_repeathandler.o vh_createMaxClusterInversion.o vh_createMaxClusterInsertion.o vh_maximalCluster.o


vh_commandlineparser.o: vh_commandlineparser.c vh_commandlineparser.h
	$(CC) $(CFLAGS) vh_commandlineparser.c 

vh_logger.o: vh_logger.c vh_logger.h
	$(CC) $(CFLAGS) vh_logger.c 

vh_main.o:	vh_main.c vh_main.h
	$(CC) $(CFLAGS) vh_main.c

vh_divethandler.o: vh_divethandler.c vh_divethandler.h
	$(CC) $(CFLAGS) vh_divethandler.c

vh_hash.o: vh_hash.c vh_hash.h
	$(CC) $(CFLAGS) vh_hash.c

vh_heap.o: vh_heap.c vh_heap.h
	$(CC) $(CFLAGS) vh_heap.c

vh_createMaxClusterDeletion.o: vh_createMaxClusterDeletion.c vh_maximalCluster.h
	$(CC) $(CFLAGS) vh_createMaxClusterDeletion.c

vh_createMaxClusterInversion.o: vh_createMaxClusterInversion.c vh_maximalCluster.h
	$(CC) $(CFLAGS) vh_createMaxClusterInversion.c

vh_createMaxClusterInsertion.o: vh_createMaxClusterInsertion.c vh_maximalCluster.h
	$(CC) $(CFLAGS) vh_createMaxClusterInsertion.c

vh_maximalCluster.o: vh_maximalCluster.c vh_maximalCluster.h
	$(CC) $(CFLAGS) vh_maximalCluster.c

vh_gaphandler.o: vh_gaphandler.c vh_gaphandler.h
	$(CC) $(CFLAGS) vh_gaphandler.c

vh_repeathandler.o: vh_repeathandler.c vh_repeathandler.h
	$(CC) $(CFLAGS) vh_repeathandler.c

clean:
	rm -f *.o VH libvh.a *~
