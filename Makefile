PYTHON = python
PIP = $(PYTHON) -m pip -v

all: build

build:
	$(PYTHON) setup.py build_ext

install:
	$(PIP) install --no-index --ignore-installed .

clean:
	rm -rf build
	rm -rf src/pydeformation/*.c

.PHONY: all build install clean
