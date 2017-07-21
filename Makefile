# This Makefile is for convenience as a reminder and shortcut for the most used commands

# Package folder
PACKAGE = pydeformation

# change to your sage command if needed
SAGE = sage

all: install test

build:
	$(SAGE) setup.py build_ext

install:
	$(SAGE) setup.py install

pip-install:
	$(SAGE) -pip install --upgrade --no-index -v .

pip-uninstall:
	$(SAGE) -pip uninstall .

pip-develop:
	$(SAGE) -pip install --upgrade -e .

sdist:
	$(SAGE) setup.py sdist

test:
	$(SAGE) setup.py test

coverage:
	$(SAGE) -coverage $(PACKAGE)/*

doc:
	cd docs && $(SAGE) -sh -c "make html"

doc-pdf:
	cd docs && $(SAGE) -sh -c "make latexpdf"

clean: clean-doc
	rm -rf build dist *.egg-info
	rm -rf $(PACKAGE)/*.c

clean-doc:
	cd docs && $(SAGE) -sh -c "make clean"

.PHONY: all build install test coverage sdist pip-install pip-uninstall pip-develop clean clean-doc doc doc-pdf
