CXX := g++
CXXFLAGS := -O3 -march=native -flto -DNDEBUG -std=c++17 -Wall -Wextra -Wshadow -Wconversion
LDFLAGS :=

BUILD_DIR := build
SRC_DIR := src

OBJS := $(BUILD_DIR)/oc_dp.o $(BUILD_DIR)/main.o

all: build_dir $(BUILD_DIR)/oc_dp

build_dir:
	mkdir -p $(BUILD_DIR)

$(BUILD_DIR)/oc_dp.o: $(SRC_DIR)/oc_dp.cpp $(SRC_DIR)/oc_dp.h
	$(CXX) $(CXXFLAGS) -c $(SRC_DIR)/oc_dp.cpp -o $@

$(BUILD_DIR)/main.o: $(SRC_DIR)/main.cpp $(SRC_DIR)/oc_dp.h
	$(CXX) $(CXXFLAGS) -c $(SRC_DIR)/main.cpp -o $@

$(BUILD_DIR)/oc_dp: $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $(BUILD_DIR)/oc_dp $(LDFLAGS)

clean:
	rm -rf $(BUILD_DIR)

.PHONY: all clean build_dir
