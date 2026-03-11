###############################################
#   RSADT + CTR FAST  —  BUILD SYSTEM (FINAL) #
###############################################

CXX       := g++
CXXFLAGS  := -O3 -std=c++17 -Wall -Wextra -fopenmp -DNDEBUG
LDFLAGS   := -fopenmp
INCFLAGS  := -Isrc

SRC_DIR   := src
BUILD_DIR := build

###############################################
#               核心库（无 main）              #
###############################################

LIB_SRCS := \
    lightOCShortest.cpp \
    weighter_policy.cpp \
    multicol_splitter.cpp \
    placement_parser.cpp \
    rowswap_ctr_expander.cpp \
    metrics_observer.cpp \
    strip_placer.cpp

LIB_OBJS := $(addprefix $(BUILD_DIR)/,$(LIB_SRCS:.cpp=.o))

###############################################
#         3 个可执行文件的 main               #
###############################################

# 1) 旧版单列（momentum）入口
MAIN_SINGLE_SRCS := main_singlecol.cpp
MAIN_SINGLE_OBJS := $(addprefix $(BUILD_DIR)/,$(MAIN_SINGLE_SRCS:.cpp=.o))
BIN_SINGLE       := $(BUILD_DIR)/oc_singlecol

# 2) 新 CTR-fast 单列入口
MAIN_POLICY_SRCS := main_policy_singlecol.cpp
MAIN_POLICY_OBJS := $(addprefix $(BUILD_DIR)/,$(MAIN_POLICY_SRCS:.cpp=.o))
BIN_POLICY       := $(BUILD_DIR)/oc_policy_singlecol

# 3) 多列封装入口
MAIN_MULTICOL_SRCS := main_multicol.cpp
MAIN_MULTICOL_OBJS := $(addprefix $(BUILD_DIR)/,$(MAIN_MULTICOL_SRCS:.cpp=.o))
BIN_MULTICOL       := $(BUILD_DIR)/oc_multicol

###############################################
#                  build rules                 #
###############################################

all: $(BIN_SINGLE) $(BIN_POLICY) $(BIN_MULTICOL)

$(BUILD_DIR):
	@mkdir -p $@

# build *.o
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $(INCFLAGS) -c $< -o $@

###############################################
#                 link executables             #
###############################################

$(BIN_SINGLE): $(LIB_OBJS) $(MAIN_SINGLE_OBJS)
	$(CXX) $^ -o $@ $(LDFLAGS)

$(BIN_POLICY): $(LIB_OBJS) $(MAIN_POLICY_OBJS)
	$(CXX) $^ -o $@ $(LDFLAGS)

$(BIN_MULTICOL): $(LIB_OBJS) $(MAIN_MULTICOL_OBJS)
	$(CXX) $^ -o $@ $(LDFLAGS)

###############################################
#                    util                      #
###############################################

clean:
	rm -rf $(BUILD_DIR)

print:
	@echo "LIB_SRCS  = $(LIB_SRCS)"
	@echo "LIB_OBJS  = $(LIB_OBJS)"

###############################################
#                END OF FILE                   #
###############################################
