CXX := g++
CXXFLAGS := -O2 -std=gnu++17 -Wall -Wextra -Wpedantic

all: rsad_benders

rsad_benders: benders_grid.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

clean:
	rm -f rsad_benders
