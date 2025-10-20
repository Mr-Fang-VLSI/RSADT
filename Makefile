# ===== Greedy (No-T) =====
GREEDY_SRCS := src/greedy_main.cpp src/greedyNoT.cpp
GREEDY_OBJS := $(GREEDY_SRCS:src/%.cpp=build/%.o)
GREEDY_BIN  := build/oc_greedy

$(GREEDY_BIN): $(GREEDY_OBJS)
	@mkdir -p build
	$(CXX) $(CXXFLAGS) -o $@ $(GREEDY_OBJS) $(LDFLAGS)
	@echo "✅ Build complete: $(GREEDY_BIN)"

# 通用规则已有的话可复用；若没有，请保留这一条
build/%.o: src/%.cpp | build
	$(CXX) $(CXXFLAGS) -c $< -o $@

greedy: $(GREEDY_BIN)

run_greedy: $(GREEDY_BIN)
	./build/oc_greedy 8 8 1 3
