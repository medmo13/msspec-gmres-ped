include options.mk

.PHONY: sdist frontend pybinding results clean


sdist: dist/msspec-$(VERSION).tar.gz


frontend: $(INSTALL_PREFIX)/bin/msspec


dist/msspec-$(VERSION).tar.gz: msspec/VERSION
	@echo "Creating Python source distribution..."
	@+$(INSIDE_VENV) pip install build && python -m build


$(INSTALL_PREFIX)/bin/msspec: msspec.sh.template msspec/VERSION
	@echo "Installing frontend command..."
	@mkdir -p $(dir $@)
	@cat $< | sed -e 's#__VENV_PATH__#$(VENV_PATH)#' > $@
	@chmod 755 $@


pybinding:
	@echo "Building Python binding for phagen and spec..."
	@+$(MAKE) -C msspec/phagen/fortran all
	@+$(MAKE) -C msspec/spec/fortran all
	@echo "$(VERSION)" > msspec/VERSION


results: msspec/results.txt


msspec/results.txt: pybinding
	@echo "Generating results for unittests"
	@coverage run --source=./ --omit=msspec/es/*,msspec/msspecgui/* msspec/create_tests_results.py
#	create the html coverage report
	@coverage html -d ../doc/source/htmlcov
	@rm .coverage


tests: pybinding
	@echo "Runing unittests"
	@python -m msspec.tests 1>/dev/null


clean::
	@echo "Cleaning all..."
	@find ./ -type f -name '*.pyc' -exec rm -f {} +
	@find ./ -type d -name '__pycache__' -exec rm -rf {} +
	@+$(MAKE) -C msspec/spec/fortran clean
	@+$(MAKE) -C msspec/phagen/fortran clean
#	remove previous sdist
	@rm -rf dist
	@rm -rf *.egg*
	@rm -f msspec/VERSION


help:
	@echo "SRC Makefile help message"

