CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o punto.o
	$(CC) punto.o random.o main.o -o main.exe

main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)

random.o : random.cpp random.h punto.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

punto.o : punto.cpp punto.h
	$(CC) -c punto.cpp -o punto.o $(CFLAGS)

clean :
	rm *.o main.exe seed.out

remove :
	rm *.txt