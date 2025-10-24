CXX := g++
CXXFLAGS := -O3 -march=native -flto -DNDEBUG -std=c++17 -Wall -Wextra -Wshadow -Wconversion

SRC_DIR := src
BUILD_DIR := build
TARGET := oc_closure

SRCS := $(SRC_DIR)/oc_closure.cpp $(SRC_DIR)/main.cpp
OBJS := $(SRCS:%.cpp=$(BUILD_DIR)/%.o)
BIN := $(BUILD_DIR)/$(TARGET)

.PHONY: all clean

all: $(BIN)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)
	mkdir -p $(BUILD_DIR)/$(SRC_DIR)

$(BUILD_DIR)/%.o: %.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BIN): $(OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@
	@echo "✅ Build complete: $(BIN)"

clean:
	rm -rf $(BUILD_DIR)
