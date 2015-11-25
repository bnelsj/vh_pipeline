CC=g++
CFLAGS= -Wall -g

bin/genotype_MultipleSamples2: src/genotype_MultipleSamples2.cpp
	$(CC) $(CFLAGS) -mcmodel=medium -o $@ $^

bin/calculateReadDepthFromBAM: src/calculateReadDepthFromBAM.cpp
	$(CC) $(CFLAGS) -o $@ $^

bin/mergeSamples: src/mergeSamples.cpp
	$(CC) $(CFLAGS) -o $@ $^

bin/spansKnownME: src/spansKnownME.cpp
	$(CC) $(CFLAGS) -o $@ $^
