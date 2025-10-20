# ============================
# Makefile for OC-CapT DAG
# ============================

CXX := g++
CXXFLAGS := -O3 -march=native -flto -DNDEBUG -std=c++17 -Wall -Wextra -Wshadow -Wconversion
# 如需调试可换成: -O0 -g -std=c++17 ...

# 使用 conda 环境中的 LEMON（不强依赖本目标）
CONDA_PREFIX ?= $(shell echo $$CONDA_PREFIX)
LEMON_PREFIX ?= $(CONDA_PREFIX)
INCLUDES := -I$(LEMON_PREFIX)/include
LDFLAGS := -L$(LEMON_PREFIX)/lib

SRC_DIR := src
BUILD_DIR := build
TARGET := oc_captdag

SRCS := $(SRC_DIR)/ocCapTDAG.cpp $(SRC_DIR)/main.cpp
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
	$(BIN) 8 20 1 16 3

envcheck:
	@if [ -z "$(CONDA_PREFIX)" ]; then \
		echo "ℹ️  Not inside conda env (OK). If you need LEMON headers: conda activate <env>"; \
	else \
		echo "✅ Conda env: $(CONDA_PREFIX)"; \
	fi
	@if [ -d "$(LEMON_PREFIX)/include/lemon" ]; then \
		echo "✅ LEMON headers found (not required for this target)"; \
	else \
		echo "ℹ️  LEMON headers not found (OK)."; \
	fi

clean:
	rm -rf $(BUILD_DIR)
