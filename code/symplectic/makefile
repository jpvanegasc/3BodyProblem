# Useful flags
DEBUG = -g -Wall -fsanitize=address -fsanitize=leak -fsanitize=undefined
VAL = --track-origins=yes --leak-check=full
CACH = --tool=cachegrind
C_STD = -std=c++11

# Dependecies
DEPS = molecular_dynamics.cpp vector3D.cpp file_handler.cpp

all : main.o 

main.o : main.cpp $(DEPS) *.h
	g++ $(C_STD) -O2 $< $(DEPS) -o $@
	time ./$@
	mkdir first; mkdir second; mkdir third
	mv earth*.txt first
	mv moon*.txt second
	mv kozai*.txt third
	bash organize.sh

.PHONY : animation
animation : animate.cpp *.h $(DEPS) animation.hpp
	g++ $(C_STD) -O2 $< $(DEPS) -o animate.o
	time ./animate.o animate | gnuplot

.PHONY : assembler
assembler : main.cpp $(DEPS) *.h
	gcc $(C_STD) -g -S main.cpp $(DEPS) -o main.S

.PHONY : parallel
parallel : main.cpp
	g++ $(C_STD) -fopenmp -O2 $< $(DEPS) 
	time ./a.out

.PHONY : debug
debug : main.cpp
	g++ $(C_STD) $(DEBUG) $< $(DEPS) 
	./a.out

.PHONY : valgrind
valgrind : main.cpp
	g++ $(C_STD) -g $< $(DEPS) 
	valgrind $(VAL) ./a.out

.PHONY : cachegrind
cachegrind : main.cpp
	g++ $(C_STD) -g $< $(DEPS) 
	valgrind $(CACH) ./a.out

.PHONY : gprof
gprof : main.cpp
	g++ $(C_STD) -Wall -pg $< $(DEPS) 
	./a.out
	gprof a.out gmon.out > analysis.txt
	most analysis.txt

.PHONY : perf
perf : main.cpp
	g++ $(C_STD) -Wall -pg $< $(DEPS) 
	perf record ./a.out ; perf report

.PHONY : clean
clean :
	rm -f *.o *.txt *.out* debug main *.data


