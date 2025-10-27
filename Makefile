# Makefile (新增 oc_td 目标，不覆盖你现有的可执行文件)
CXX      := g++
CXXFLAGS := -O3 -march=native -flto -DNDEBUG -std=c++17 -Wall -Wextra -Wshadow -Wconversion -fopenmp
INC      := -Isrc
LDFLAGS  := -fopenmp

BUILD_DIR := build
SRC_DIR   := src

OBJS := $(BUILD_DIR)/src/lightOCShortest.o $(BUILD_DIR)/src/main_td.o
BIN  := $(BUILD_DIR)/oc_td

all: $(BIN)

$(BUILD_DIR)/src/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(BUILD_DIR)/src
	$(CXX) $(CXXFLAGS) $(INC) -c $< -o $@

$(BIN): $(OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

clean:
	rm -rf $(BUILD_DIR)

.PHONY: all clean
