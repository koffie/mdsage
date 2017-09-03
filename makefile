# This Makefile is for convenience as a reminder and shortcut for the most used commands

# Package folder
PACKAGE = mdsage

# change to your sage command if needed
ifndef SAGE
	SAGE = /Applications/sagedev/sage
	export SAGE
endif

all: install test

install:
	$(SAGE) -pip install $(PIP_ARGS) --ignore-installed -v .

uninstall:
	$(SAGE) -pip uninstall .

develop:
	$(SAGE) -pip install --upgrade -e .

test:
	$(SAGE) setup.py test

coverage:
	$(SAGE) -coverage $(PACKAGE)/*

doc:	install
	cd docs && $(SAGE) -sh -c "make html"

doc-pdf:
	cd docs && $(SAGE) -sh -c "make latexpdf"

clean:	clean-test clean-doc

clean-test: clean-build
	find $(PACKAGE) -name "*.pyc" -type f -print0 -delete | xargs -0 echo "removing"
	find $(PACKAGE) -name "*.so" -type f -print0 -delete | xargs -0 echo "removing"

clean-build:
	find $(PACKAGE) -name "*.c" -type f -print0 -delete | xargs -0 echo "removing"
	rm -rf build/*

clean-doc:
	cd docs && $(SAGE) -sh -c "make clean"

.PHONY: all install develop test coverage doc doc-pdf
