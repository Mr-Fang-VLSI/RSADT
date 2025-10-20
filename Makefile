CXX      := g++
CXXFLAGS := -O3 -std=c++17 -Wall -Wextra -Wshadow -Wconversion -DNDEBUG -fopenmp
LDFLAGS  := -fopenmp
TARGET   := build/oc_maxt

SRCS := src/main.cpp src/lightOCMaxT.cpp src/networkOCMaxT.cpp
OBJS := $(SRCS:src/%.cpp=build/%.o)

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

# AddressSanitizer 构建（如需快速定位越界）
asan:
	@mkdir -p build
	$(CXX) -O1 -g -std=c++17 -fsanitize=address -fno-omit-frame-pointer -Wall -Wextra -Wshadow -Wconversion -o $(TARGET)_asan $(SRCS) $(LDFLAGS)
	@echo "🔎 Built ASAN binary: $(TARGET)_asan"

run:
	./build/oc_maxt 8 8 1 20 0 1 3 3

.PHONY: all clean run asan
