# CC=/usr/bin/g++
CC=g++
#~ CC=clang++
CFLAGS=  -Wall  -Ofast -std=c++11 -march=native -flto -pipe -funit-at-a-time -fomit-frame-pointer   -Wfatal-errors
LDFLAGS=-flto -lpthread


ifeq ($(gprof),1)
CFLAGS=-std=c++0x -pg -O4   -march=native
LDFLAGS=-pg
endif

ifeq ($(valgrind),1)
CFLAGS=-std=c++0x -O4 -g
LDFLAGS=-g
endif



EXEC=TheCreator

all: $(EXEC)


TheCreator: videl.o  
	$(CC) -o $@ $^ $(LDFLAGS)
	
videl.o: videl.cpp  
	$(CC) -o $@ -c $< $(CFLAGS)


clean:
	rm -rf *.o
	rm -rf $(EXEC)


rebuild: clean $(EXEC)
