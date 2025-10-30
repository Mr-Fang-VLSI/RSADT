# ===== Compiler & Flags =====
CXX      := g++
CXXFLAGS := -O3 -std=c++17 -Wall -Wextra -DNDEBUG
OMPFLAG  := -fopenmp

SRC_DIR   := src
BUILD_DIR := build

# ---- single column executable ----
SRC_SINGLECOL := $(SRC_DIR)/lightOCShortest.cpp $(SRC_DIR)/main_singlecol.cpp
HDR_SINGLECOL := $(SRC_DIR)/lightOCShortest.h $(SRC_DIR)/momentum_weighter.h

OBJ_SINGLECOL := $(SRC_SINGLECOL:.cpp=.o)
TARGET_SINGLECOL := $(BUILD_DIR)/oc_singlecol

# ===== Default =====
all: $(TARGET_SINGLECOL)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(TARGET_SINGLECOL): $(OBJ_SINGLECOL) | $(BUILD_DIR)
	@echo "[Link SingleCol] -> $@"
	$(CXX) $(CXXFLAGS) $(OMPFLAG) -o $@ $(OBJ_SINGLECOL)

$(SRC_DIR)/%.o: $(SRC_DIR)/%.cpp $(HDR_SINGLECOL)
	@echo "[Compile] $<"
	$(CXX) $(CXXFLAGS) $(OMPFLAG) -c $< -o $@

clean:
	@echo "[Clean]"
	rm -rf $(BUILD_DIR) $(SRC_DIR)/*.o

rebuild: clean all
