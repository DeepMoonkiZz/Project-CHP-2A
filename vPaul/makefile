CC=gcc
CFLAGS=-Wall -g
OBJS=mod_function.o mod_matrix.o mod_gradient.o mod_display.o

all: main

main: main.c $(OBJS)
	$(CC) $(CFLAGS) $^ -o $@

mod_function.o: mod_function.c mod_function.h
	$(CC) $(CFLAGS) -c $< -o $@

mod_matrix.o: mod_matrix.c mod_matrix.h
	$(CC) $(CFLAGS) -c $< -o $@

mod_gradient.o: mod_gradient.c mod_gradient.h
	$(CC) $(CFLAGS) -c $< -o $@

mod_display.o: mod_display.c mod_display.h
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f *.o main
