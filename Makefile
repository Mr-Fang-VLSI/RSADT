CXX := g++
MODE ?= release

ifeq ($(MODE),debug)
  CXXFLAGS := -O0 -g -std=c++17 -Wall -Wextra -Wshadow -Wconversion \
              -fsanitize=address,undefined -fno-omit-frame-pointer -D_GLIBCXX_ASSERTIONS
else
  CXXFLAGS := -O3 -march=native -std=c++17 -Wall -Wextra -Wshadow -Wconversion -DNDEBUG
  # 如需再开 LTO：在完全稳定后把下一行注释解除
  # CXXFLAGS += -flto
endif

SRC_DIR := src
BUILD_DIR := build
TARGET := oc_captdag

SRCS := $(SRC_DIR)/ocCapTDAG.cpp $(SRC_DIR)/main.cpp
OBJS := $(SRCS:%.cpp=$(BUILD_DIR)/%.o)
BIN := $(BUILD_DIR)/$(TARGET)

.PHONY: all clean run

all: $(BIN)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR) $(BUILD_DIR)/$(SRC_DIR)

$(BUILD_DIR)/%.o: %.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BIN): $(OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@
	@echo "✅ Build complete: $(BIN)"

run: all
	$(BIN) 8 8 1 64 2 4

clean:
	rm -rf $(BUILD_DIR) run.log
