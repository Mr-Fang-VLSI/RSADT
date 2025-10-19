# ============================
# Makefile for OC-MaxT
# ============================

CXX := g++
CXXFLAGS := -O3 -march=native -flto -DNDEBUG -std=c++17 -Wall -Wextra -Wshadow -Wconversion -fopenmp

SRC_DIR := src
BUILD_DIR := build
TARGET := oc_maxt

SRCS := $(SRC_DIR)/lightOCMaxT.cpp $(SRC_DIR)/main.cpp
OBJS := $(SRCS:%.cpp=$(BUILD_DIR)/%.o)
BIN := $(BUILD_DIR)/$(TARGET)

.PHONY: all clean run

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
