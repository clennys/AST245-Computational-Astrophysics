.PHONY: run build

run: build
	./build/n-body data/data.txt

build:
	cmake -B build && cmake --build build --parallel

clean:
	rm -r build
	rm *.png
