# ===========================
# Makefile for Greedy-OC demo
# ===========================
CXX := g++
CXXFLAGS := -O3 -std=c++17 -Wall -Wextra -Wshadow -Wconversion -DNDEBUG

SRC_DIR := src
BUILD_DIR := build
TARGET := greedy_oc

SRCS := $(SRC_DIR)/GreedyOCColumnPlacer.cpp $(SRC_DIR)/main.cpp
OBJS := $(SRCS:%.cpp=$(BUILD_DIR)/%.o)
BIN := $(BUILD_DIR)/$(TARGET)

.PHONY: all run clean

all: $(BIN)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)
	mkdir -p $(BUILD_DIR)/$(SRC_DIR)

$(BUILD_DIR)/%.o: %.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BIN): $(OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@
	@echo "✅ Build complete: $(BIN)"

run: all
	$(BIN)

clean:
	rm -rf $(BUILD_DIR)
