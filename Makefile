all: pert.o

testDriver: pert.o testDriver.c
	gcc -Wall -std=c99 -pendantic testDriver.c pert.o -lm -o testDriver
	
pert.o: pert.h pert.c
	gcc -Wall -std=c99 -pendantic -c pert.c