.PHONY: run build

run: build
	./build/n-body data/data.txt

build:
	mkdir -p plots
	cmake -B build && cmake --build build --parallel

clean:
	rm -rf build plots
