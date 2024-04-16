CCOMPILE=mpic++
PLATFORM=Linux-amd64-64
CPPFLAGS= -I$(JAVA_HOME)/include -I$(JAVA_HOME)/include/linux -I ~/pregelplus-code -I ~/blogel-code
#-I$(HADOOP_HOME)/src/c++/libhdfs 
# LIB = -L$(HADOOP_HOME)/c++/$(PLATFORM)/lib
LDFLAGS =  -Wno-deprecated -O2
#-lhdfs

all: run

run: run.cpp
	$(CCOMPILE) run.cpp $(CPPFLAGS) $(LDFLAGS)  -o run
    # $(LIB) 

clean:
	-rm run
