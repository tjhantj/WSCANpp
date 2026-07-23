# Compiler
CXX = g++

# Compile options
CXXFLAGS = -std=c++17 -O2 -Wall -Wextra -fopenmp

# Target executable
TARGET = wscan

# Source files
SRCS = \
	code/main.cpp \
	code/graph.cpp \
	code/clustering.cpp \
	code/similarity.cpp \
	code/args.cpp \
	code/metrics.cpp \
	code/utils.cpp \
	code/perturb_edge.cpp

# Default rule
all: $(TARGET)

# Link
$(TARGET): $(SRCS)
	$(CXX) $(CXXFLAGS) $(SRCS) -o $(TARGET)

# Clean rule
clean:
	rm -f $(TARGET)

# Rebuild rule
rebuild: clean all
