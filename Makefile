CC=g++
CFLAGS= -Wall -g -O3

bin/genotypeSamplesVH: src/genotypeSamplesVH.cpp
	$(CC) $(CFLAGS) -std=c++0x -I tabixpp/htslib/ -I tabixpp/  src/genotypeSamplesVH.cpp src/split.cpp tabixpp/tabix.a tabixpp/htslib/libhts.a -lm -lz -lpthread -o $@

bin/genotype_MultipleSamples2: src/genotype_MultipleSamples2.cpp
	$(CC) $(CFLAGS) -mcmodel=medium -o $@ $^

bin/calculateReadDepthFromBAM: src/calculateReadDepthFromBAM.cpp
	$(CC) $(CFLAGS) -o $@ $^

bin/mergeSamples: src/mergeSamples.cpp
	$(CC) $(CFLAGS) -o $@ $^

bin/spansKnownME: src/spansKnownME.cpp
	$(CC) $(CFLAGS) -o $@ $^
