CC = g++
CFLAGS = -Wall -O3 --std=c++11

punto.o : punto.cpp punto.h random.o
	$(CC) -c punto.cpp -o punto.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o main.exe seed.out
remove :
	rm *.txt