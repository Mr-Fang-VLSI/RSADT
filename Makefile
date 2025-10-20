# ======== Compiler & Flags ========
CXX      := g++
CXXFLAGS := -O3 -std=c++17 -Wall -Wextra -Wshadow -Wconversion -DNDEBUG -fopenmp
LDFLAGS  := -fopenmp
TARGET   := build/oc_maxt

# ======== Source Files ========
# 注意：networkOCMax.cpp 是 RCDC 的实现文件（实现 networkOCMaxT 类）
SRCS := src/main.cpp src/lightOCMaxT.cpp src/networkOCMaxT.cpp
OBJS := $(SRCS:src/%.cpp=build/%.o)

# ======== Build Rules ========
all: $(TARGET)

$(TARGET): $(OBJS)
	@mkdir -p build
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDFLAGS)
	@echo "✅ Build complete: $(TARGET)"

build/%.o: src/%.cpp | build
	$(CXX) $(CXXFLAGS) -c $< -o $@

build:
	mkdir -p build

clean:
	rm -rf build
	@echo "🧹 Cleaned."

run:
	./build/oc_maxt 8 8 1 16 0 1 3 3

.PHONY: all clean run
