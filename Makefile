CXX := g++
CXXFLAGS := -O3 -std=c++17 -Wall -Wextra

PKG_CFLAGS := $(shell pkg-config --cflags ortools 2>/dev/null)
PKG_LIBS   := $(shell pkg-config --libs   ortools 2>/dev/null)

# 如果没有 pkg-config 支持，请导出 OR_TOOLS_DIR=/your/prefix
ifdef OR_TOOLS_DIR
  INC := -I$(OR_TOOLS_DIR)/include
  LIB := -L$(OR_TOOLS_DIR)/lib -lortools -lpthread
else
  INC := $(PKG_CFLAGS)
  LIB := $(PKG_LIBS)
  # 回退：某些发行版只需链接 -lortools
  ifeq ($(LIB),)
    LIB := -lortools -lpthread
  endif
endif

all: rsad_cp

rsad_cp: rsad_cp.cpp
	$(CXX) $(CXXFLAGS) $(INC) $< -o $@ $(LIB)

run: rsad_cp
	./rsad_cp

clean:
	rm -f rsad_cp
