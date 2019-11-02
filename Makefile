.PHONY: all build test clean

all: build

build:
	mkdir -p build
	cd build && cmake .. && make

test: build
	cd build && make test

clean:
	rm -rf build
