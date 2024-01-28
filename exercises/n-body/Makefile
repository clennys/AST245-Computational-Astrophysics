.PHONY: build run run-simple

run: build
	./build/n-body data/data.txt

run-simple: build
	./build/n-body data/data-red.txt

build:
	mkdir -p plots/png
	cmake -B build && cmake --build build --parallel

clean:
	rm -rf build plots
