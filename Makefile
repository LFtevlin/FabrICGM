CXX := g++
CXXFLAGS := -std=c++17 -O3 -Wall -Wextra

HDF5_INCLUDE := $(shell find /usr/include /usr/local/include -name H5Cpp.h 2>/dev/null | head -n 1)

ifeq ($(HDF5_INCLUDE),)
$(error HDF5 C++ header H5Cpp.h not found)
endif

HDF5_INC_DIR := $(dir $(HDF5_INCLUDE))

HDF5_LIB := $(shell find /usr/lib /usr/local/lib -name libhdf5_cpp.so* 2>/dev/null | head -n 1)

ifeq ($(HDF5_LIB),)
$(error HDF5 C++ library not found)
endif

HDF5_LIB_DIR := $(dir $(HDF5_LIB))

CXXFLAGS += -I$(HDF5_INC_DIR)

LDFLAGS := -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5_cpp -lhdf5 -lfftw3

SOURCES := src/main.cpp src/simulation.cpp src/cooling.cpp src/parameters.cpp src/galaxy_catalogue.cpp src/turbulence.cpp src/functions.cpp src/sampler.cpp

OBJECTS := $(SOURCES:.cpp=.o)

TARGET := FabriCGM

all: check $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(OBJECTS) -o $(TARGET) $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

check:
	@echo "Checking compiler..."
	@command -v $(CXX) >/dev/null || (echo "ERROR: g++ not found"; exit 1)
	@echo "Using compiler:"
	@$(CXX) --version | head -n 1
	@echo ""
	@echo "Found HDF5:"
	@echo "Include: $(HDF5_INC_DIR)"
	@echo "Library: $(HDF5_LIB_DIR)"

clean:
	rm -f $(OBJECTS)
	rm -f $(TARGET)