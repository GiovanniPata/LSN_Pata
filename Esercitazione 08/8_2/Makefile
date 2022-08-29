CC = g++
CFLAGS = -Wall -O3 --std=c++11

SA.exe : SA.o random.o
	$(CC) random.o SA.o -o SA.exe
SA.o : SA.cpp
	$(CC) -c SA.cpp -o SA.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o SA.exe seed.out
remove :
	rm *.txt
