# Compiler
CC = g++

# Compiler flags
CFLAGS = -Wall -std=c++17

# Executable name
EXEC = pic_simulation

# Source files
SOURCES = pic_main.cpp

# Object files
OBJECTS = $(SOURCES:.cpp=.o)

# Header files
HEADERS = pic_main.h pic_math.h pic_structs.h

# Default target
all: $(EXEC)

# Rule to link the program
$(EXEC): $(OBJECTS)
	$(CC) $(CFLAGS) -o $(EXEC) $(OBJECTS)

# Rule for object files
%.o: %.cpp $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

# Clean target
clean:
	rm -f $(EXEC) $(OBJECTS)
