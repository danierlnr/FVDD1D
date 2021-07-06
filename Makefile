CC := g++
#CFLAGS     := -c -O2 -Wall -s
#CFLAGS     := -c -g -Wall 
CFLAGS     := -c -g 
#CFLAGS     := -g -pg -c -Wall 
#LFLAGS     := -pg -o
LFLAGS     := -o
SOURCES    := $(wildcard src/*.cpp)
LIBPATH    := -L/usr/lib64
#LIBRARIES  := $(wildcard lib/*.a) $(wildcard lib/*.so)
#LIBRARIES  := -L$LIBPATH $(wildcard lib/*.a) -ldl -llapack # For levmar using lapack
LIBRARIES  := $(wildcard lib/*.a)
OBJECTS    := $(SOURCES:.cpp=.o)
EXECUTABLE := FVDD1D

all: $(EXECUTABLE) clean
#all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) $(LIBRARIES)
	$(CC) $(OBJECTS) $(LIBRARIES) $(LFLAGS) $(EXECUTABLE)
	
.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
	
clean:
	rm $(OBJECTS)
