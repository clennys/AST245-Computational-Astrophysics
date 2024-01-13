.PHONY: run build

run: build
	./build/n-body data/data.txt

build:
	cmake -B build && cmake --build build --parallel

clean:
	rm -rf build *.png *.jpg
