# Simple Makefile for IDA* single-column placer
CXX ?= g++
CXXFLAGS ?= -O3 -std=gnu++17 -Wall -Wextra -Wno-sign-compare -Wno-unused-parameter
LDFLAGS ?= 

SRC_DIR := .
BUILD_DIR := build

TARGET := $(BUILD_DIR)/ida_star
SRCS := IdaStarColumnPlacer.cpp main.cpp
OBJS := $(patsubst %.cpp,$(BUILD_DIR)/%.o,$(SRCS))
DEPS := $(OBJS:.o=.d)

all: $(TARGET)

$(BUILD_DIR):
	@mkdir -p $(BUILD_DIR)

$(BUILD_DIR)/%.o: %.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -MMD -MP -c $< -o $@

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $@ $(LDFLAGS)

clean:
	rm -rf $(BUILD_DIR)

.PHONY: all clean
-include $(DEPS)
