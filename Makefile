CXX ?= g++
# 屏蔽 timespec_get 相关的 LIB EXT1；避免 <ctime> 干扰
CXXFLAGS ?= -O3 -std=c++17 -DNDEBUG -D__STDC_WANT_LIB_EXT1__=0
LDFLAGS ?=

OBJS = dp_oc.o main_dp.o

all: dp_model

dp_model: $(OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

dp_oc.o: dp_oc.cpp dp_oc.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

main_dp.o: main_dp.cpp dp_oc.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f *.o dp_model
