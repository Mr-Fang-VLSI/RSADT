# ============================================================
#               Global build configuration
# ============================================================
CXX       := g++
CXXFLAGS  := -O3 -std=c++17 -Wall -Wextra -DNDEBUG
OMPFLAGS  := -fopenmp
SRC_DIR   := src
BUILD_DIR := build

# ============================================================
#                   Source group definitions
# ============================================================

# ---- 1. Single-column TD baseline ----
SINGLE_SRCS := \
  $(SRC_DIR)/lightOCShortest.cpp \
  $(SRC_DIR)/placement_parser.cpp \
  $(SRC_DIR)/main_singlecol.cpp

SINGLE_OBJS := $(SINGLE_SRCS:.cpp=.o)

# ---- 2. Policy framework (multi-weighting strategies) ----
POLICY_SRCS := \
  $(SRC_DIR)/lightOCShortest.cpp \
  $(SRC_DIR)/weighter_policy.cpp \
  $(SRC_DIR)/main_policy_singlecol.cpp

POLICY_OBJS := $(POLICY_SRCS:.cpp=.o)

# ---- 3. Metrics observer (visualization & logging) ----
OBS_SRCS := \
  $(SRC_DIR)/lightOCShortest.cpp \
  $(SRC_DIR)/metrics_observer.cpp \
  $(SRC_DIR)/main_observe.cpp

OBS_OBJS := $(OBS_SRCS:.cpp=.o)

# ============================================================
#                     Build rules
# ============================================================
.PHONY: all clean rebuild rerun8x8

all: $(BUILD_DIR)/oc_singlecol \
     $(BUILD_DIR)/oc_policy_singlecol \
     $(BUILD_DIR)/oc_observe

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# ---- oc_singlecol ----
$(BUILD_DIR)/oc_singlecol: $(SINGLE_OBJS) | $(BUILD_DIR)
	@echo "[Link] -> $@"
	$(CXX) $(CXXFLAGS) $(OMPFLAGS) -o $@ $(SINGLE_OBJS)

# ---- oc_policy_singlecol ----
$(BUILD_DIR)/oc_policy_singlecol: $(POLICY_OBJS) | $(BUILD_DIR)
	@echo "[Link] -> $@"
	$(CXX) $(CXXFLAGS) $(OMPFLAGS) -o $@ $(POLICY_OBJS)

# ---- oc_observe ----
$(BUILD_DIR)/oc_observe: $(OBS_OBJS) | $(BUILD_DIR)
	@echo "[Link] -> $@"
	$(CXX) $(CXXFLAGS) $(OMPFLAGS) -o $@ $(OBS_OBJS)

# ---- Compile rule for all .cpp ----
$(SRC_DIR)/%.o: $(SRC_DIR)/%.cpp
	@echo "[Compile] $<"
	$(CXX) $(CXXFLAGS) $(OMPFLAGS) -c $< -o $@

# ============================================================
#                     Utility targets
# ============================================================
clean:
	@echo "[Clean]"
	rm -rf $(BUILD_DIR) $(SRC_DIR)/*.o

rebuild: clean all

# ---- Batch rerun for 8x8, T=8..20 ----
rerun8x8: $(BUILD_DIR)/oc_singlecol
	@chmod +x ./run_rerun_8x8.sh || true
	./run_rerun_8x8.sh
