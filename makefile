# ==============================
# Makefile for Param-Closure MCF
# ==============================
CXX := g++
CXXFLAGS := -O3 -std=c++17 -Wall -Wextra -Wshadow -Wconversion -DNDEBUG

# conda lemon
CONDA_PREFIX ?= $(shell echo $$CONDA_PREFIX)
LEMON_PREFIX ?= $(CONDA_PREFIX)
INCLUDES := -I$(LEMON_PREFIX)/include
LDFLAGS  := -L$(LEMON_PREFIX)/lib -lemon

SRC_DIR := src
BUILD_DIR := build
TARGET := param_closure

SRCS := $(SRC_DIR)/ParamClosureOCPlacer.cpp $(SRC_DIR)/main.cpp
OBJS := $(SRCS:%.cpp=$(BUILD_DIR)/%.o)
BIN := $(BUILD_DIR)/$(TARGET)

.PHONY: all run clean envcheck

all: envcheck $(BIN)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)
	mkdir -p $(BUILD_DIR)/$(SRC_DIR)

$(BUILD_DIR)/%.o: %.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

$(BIN): $(OBJS)
	$(CXX) $(CXXFLAGS) $^ $(LDFLAGS) -o $@
	@echo "✅ Build complete: $(BIN)"

run: all
	$(BIN)

clean:
	rm -rf $(BUILD_DIR)

envcheck:
	@if [ -z "$(CONDA_PREFIX)" ]; then \
		echo "⚠️  Not inside conda env. Run: conda activate rsad_mcmf"; \
	else \
		echo "✅ Conda env: $(CONDA_PREFIX)"; \
	fi
	@if [ ! -f "$(LEMON_PREFIX)/include/lemon/preflow.h" ]; then \
		echo "⚠️  lemon headers not found. Try: conda install -c conda-forge lemon"; \
	else \
		echo "✅ lemon found in $(LEMON_PREFIX)"; \
	fi
