CC=g++ -O4
CFLAGS=-c -std=c++11
INC=-I/opt/local/include
LIBS=-L/opt/local/lib
LDLIBS=-lgsl -lgslcblas -lm

OBJS =  \
A_NumDenom.o \
AdiabLR.o \
main.o

all:main

main:$(OBJS)
	$(CC) $(LIBS)  $(OBJS) -o EOB_ihes.out $(LDLIBS)
main.o: main.cpp
	$(CC) $(CFLAGS) $(INC) main.cpp
A_NumDenom.o: A_NumDenom.cpp
	$(CC) $(CFLAGS) $(INC) A_NumDenom.cpp
AdiabLR.o: AdiabLR.cpp
	$(CC) $(CFLAGS) $(INC) AdiabLR.cpp

clean:
	rm -rf *.o EOB_ihes.out