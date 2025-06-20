CXX = g++

CXXFLAGS = -std=c++17 -O3 -g -Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include

LDFLAGS = -L/opt/homebrew/opt/libomp/lib -lomp

TARGET = optimizer

SRCS = MeanShift.cpp utils.cpp

# .cpp to .o
OBJS = $(SRCS:.cpp=.o)

all: $(TARGET)

# link
$(TARGET): $(OBJS)
	$(CXX) -o $@ $(OBJS) $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -f $(TARGET) $(OBJS)

.PHONY: all clean
