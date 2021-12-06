SHELL = /bin/bash
CXX = g++
CXX_FLAGS = -std=c++17 -O2 -Wall -funroll-loops
INCLUDES = -I./src/

test: src/test.cpp src/obj.cpp src/mesh.cpp
	$(CXX) $(CXX_FLAGS) $(INCLUDES) -o $@ $^

.PHONY: clean
clean:
	rm -rf *.o test
