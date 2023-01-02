CXX:=g++
CXX_FLAGS:=-std=c++11 -ggdb

BIN:=bin
SRC:=src
INCLUDE:=include

LIBRARIES:=
EXECUTABLE:=main

MKL_LIB:=-llapack -lblas

all:$(BIN)/$(EXECUTABLE)

run: clean all
	clear
	./$(BIN)/$(EXECUTABLE)

$(BIN)/$(EXECUTABLE): $(SRC)/*.cpp
	$(CXX) $(CXX_FLAGS) $(MKL_LIB) -I$(INCLUDE) $^ -o $@ $(LIBRARIES)

clean:
		-rm $(BIN)/*