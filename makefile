all: main.o eig.o
	gcc main.o eig.o -o eig -lm -lrt 
debug: main.o eig.o
	gcc main.o eig.o -o eig -lm -lrt -g 
main.o: main.c eig.h
	gcc -c main.c -std=gnu99
eig.o: eig.c eig.h
	gcc -c eig.c -std=gnu99
clean:
	rm -rf *.o
