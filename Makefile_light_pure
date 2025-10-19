# ============================
# Makefile for LightMCMF
# ============================

CXX := g++
CXXFLAGS := -O2 -std=c++17 -Wall -Wextra -Wshadow -Wconversion
CONDA_PREFIX ?= $(shell echo $$CONDA_PREFIX)
LEMON_PREFIX ?= $(CONDA_PREFIX)
INCLUDES := -I$(LEMON_PREFIX)/include
LDFLAGS := -L$(LEMON_PREFIX)/lib -lemon

SRC_DIR := src
BUILD_DIR := build
TARGET := light_mcmf

SRCS := $(SRC_DIR)/lightPureMcmf.cpp $(SRC_DIR)/main.cpp
OBJS := $(SRCS:%.cpp=$(BUILD_DIR)/%.o)
BIN := $(BUILD_DIR)/$(TARGET)

.PHONY: all clean run envcheck

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

envcheck:
	@if [ -z "$(CONDA_PREFIX)" ]; then \
		echo "⚠️  Not inside conda environment."; \
	else \
		echo "✅ Conda env: $(CONDA_PREFIX)"; \
	fi
	@if pkg-config --exists lemon; then \
		echo "✅ LEMON found"; \
	else \
		echo "⚠️  LEMON not found. Try: conda install -c conda-forge lemon"; \
	fi

clean:
	rm -rf $(BUILD_DIR)
