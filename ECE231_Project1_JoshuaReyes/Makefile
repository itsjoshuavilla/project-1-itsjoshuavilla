# Makefile 
CXX      := g++
CXXFLAGS := -O2 -std=c++17 -Wall -Wextra

TARGET := project1
SRC    := main.cpp

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $< -o $@

run: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(TARGET)

.PHONY: all run clean
