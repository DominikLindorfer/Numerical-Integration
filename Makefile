CC=g++
CFLAGS=-c -Wall -fopenmp -std=c++1y
LDFLAGS=-fopenmp
SOURCES=MultiDim_Gauss.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=MultiDim_Gauss.out

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@


