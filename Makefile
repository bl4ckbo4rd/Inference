EXEC=Indenteur
CC=g++

S=GRAPH

BITS=$(getconf LONG_BIT)
ifeq ($(BITS),64)
	FBITS=-m64
else
	FBITS=
endif

#per compilare con openmp
#FLAGS=$(FBITS) -std=c++11 -O3 -fopenmp
#per compilare senza
FLAGS=$(FBITS) -std=c++11 -O3 

all : Graph

Graph : main.o problems.o inverse.o 
	$(CC) $(FLAGS) -o bin/Graph main.o problems.o inverse.o

main.o : $(S)/main.cpp $(S)/problems.h
	$(CC) -c $(FLAGS) $(S)/main.cpp

problems.o : $(S)/problems.cpp $(S)/problems.h $(S)/inverse.h
	$(CC) -c $(FLAGS) $(S)/problems.cpp

inverse.o : $(S)/inverse.cpp $(S)/inverse.h
	$(CC) -c $(FLAGS) $(S)/inverse.cpp



#per cancellare i file oggetto fai make clean
clean :
	@rm *.o

rmexec :
	@rm $(EXEC)
