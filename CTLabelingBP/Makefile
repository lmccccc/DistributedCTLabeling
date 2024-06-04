CCOMPILE=mpic++
PLATFORM=Linux-amd64-64
CPPFLAGS= -I$(JAVA_HOME)/include -I$(JAVA_HOME)/include/linux -I ~/pregelplus-code -I ~/blogel-code
LDFLAGS =  -Wno-deprecated -O2

all: clean run

run: run.cpp
	$(CCOMPILE) -fopenmp run.cpp $(CPPFLAGS) $(LDFLAGS) -o run

clean:
	-rm run
