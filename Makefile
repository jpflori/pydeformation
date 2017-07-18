all: cythonized

sage_local:
	@if [ -z $(SAGE_LOCAL) ]; then echo "make should be run from within a Sage shell" >&2; exit 1; fi

CYTHON_SOURCES = $(wildcard *.pyx) $(wildcard *.pxd)

cythonized: $(CYTHON_SOURCES)
	python setup.py build_ext
	cp cythonized/lib.*/*.so ./

clean:
	rm -rf cythonized
	rm -rf *.so

.PHONY: clean sage_local 
