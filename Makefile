# ======== Compiler & Flags ========
CXX      := g++
CXXFLAGS := -O3 -std=c++17 -Wall -Wextra -Wshadow -Wconversion -DNDEBUG -fopenmp
LDFLAGS  := -fopenmp
TARGET   := build/oc_maxt

# ======== Source Files ========
SRCS := src/main.cpp src/lightOCMaxT.cpp src/networkOCMaxT.cpp
OBJS := $(SRCS:src/%.cpp=build/%.o)

# ======== Build Rules ========
all: $(TARGET)

# 链接
$(TARGET): $(OBJS)
	@mkdir -p build
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDFLAGS)
	@echo "✅ Build complete: $(TARGET)"

# 编译 src 下的所有 .cpp 文件到 build 目录
build/%.o: src/%.cpp | build
	$(CXX) $(CXXFLAGS) -c $< -o $@

# 若 build 目录不存在则创建
build:
	mkdir -p build

clean:
	rm -rf build
	@echo "🧹 Cleaned."

run:
	./build/oc_maxt 8 8 1 16 4 1 3 3

.PHONY: all clean run
