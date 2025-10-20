# ============================
# Makefile for OC-CapT-DAG (APT)
# ============================

CXX := g++
MODE ?= release

# Common
CXXSTD := -std=c++17
WARN   := -Wall -Wextra -Wshadow -Wconversion
INCS   :=
LIBS   :=
DEFS   :=

# Flags per mode
ifeq ($(MODE),debug)
  CXXFLAGS := -O0 -g $(CXXSTD) $(WARN) -fno-omit-frame-pointer -D_GLIBCXX_ASSERTIONS
  SAN      := -fsanitize=address,undefined
else
  CXXFLAGS := -O3 -march=native -flto -DNDEBUG $(CXXSTD) $(WARN)
  SAN      :=
endif

SRC_DIR   := src
BUILD_DIR := build
TARGET    := oc_captdag

SRCS := $(SRC_DIR)/ocCapTDAG.cpp $(SRC_DIR)/main.cpp
OBJS := $(SRCS:%.cpp=$(BUILD_DIR)/%.o)
BIN  := $(BUILD_DIR)/$(TARGET)

.PHONY: all clean run

all: $(BIN)

$(BUILD_DIR):
	@mkdir -p $(BUILD_DIR)
	@mkdir -p $(BUILD_DIR)/$(SRC_DIR)

$(BUILD_DIR)/%.o: %.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $(DEFS) $(INCS) -c $< -o $@

$(BIN): $(OBJS)
	$(CXX) $(CXXFLAGS) $^ $(LIBS) $(SAN) -o $@
	@echo "✅ Build complete: $(BIN)  (MODE=$(MODE))"

run: all
	$(BIN) 8 8 1 16 3 0

clean:
	rm -rf $(BUILD_DIR) *.dSYM *.prof *.out run.log
