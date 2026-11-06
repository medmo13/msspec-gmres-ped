include src/options.mk


.PHONY: pybinding install devel venv doc clean _attrdict


pybinding:
	@+$(INSIDE_VENV) $(MAKE) -C src pybinding


venv:
ifeq ($(NO_VENV),0)
	@virtualenv --python=$(PYTHON_EXE) --prompt="(msspec-$(VERSION)) " $(VENV_PATH)
	$(INSIDE_VENV) python -m ensurepip --upgrade
endif

# wget https://bootstrap.pypa.io/get-pip.py && \
# python get-pip.py                         && \
# rm -f get-pip.py
# pip install --upgrade setuptools          && \
# pip install -r src/pip.freeze             && \


install: venv pybinding wx
	@+$(INSIDE_VENV) $(MAKE) -C src sdist
	@+$(INSIDE_VENV) $(MAKE) -C src frontend
	@+$(INSIDE_VENV) pip install src/dist/msspec-$(VERSION)*.whl
	@echo "Do not forget to check that $(INSTALL_PREFIX)/bin is set in your \$$PATH"


devel: VENV_PATH = ./_venv
devel: venv pybinding wx
	@$(INSIDE_VENV) pip install -e src/


light: VENV_PATH = ./_venv
light: venv
	@$(INSIDE_VENV) pip install src/

nogui: VENV_PATH = ./_venv
nogui: venv pybinding
	@$(INSIDE_VENV) pip install -e src/


_attrdict:
	# Check if virtualenv python version > 3.3.0
	# If so, install the patched version of attrdict used to build the version 4.2.0 of wxPython
	@$(INSIDE_VENV) if `python -c "import sys; exit(sys.version_info > (3,3))"`; then \
		pip install --no-cache attrdict; \
	else \
		pip install thirdparty/attrdict-2.0.1.tar.gz; \
	fi


_build_wx/wxPython.target: _attrdict
	@$(INSIDE_VENV) echo "Building wxPython for your `python --version 2>&1` under Linux $(DISTRO_RELEASE)..."
	# Create a folder to build wx into
	@mkdir -p _build_wx
	# download the wheel or the source if it cannot find a wheel
	$(INSIDE_VENV) cd _build_wx && pip download -f https://extras.wxpython.org/wxPython4/extras/linux/gtk3/$(DISTRO_RELEASE) wxPython
	# Build the source if a tar.gz was downloaded
	@$(INSIDE_VENV) cd _build_wx && \
          if [ -e wxPython*.tar.gz ]; then \
            tar -x --skip-old-files -vzf wxPython*.tar.gz; \
            cd `ls -d wxPython*/`; \
            pip install requests sip; \
            python build.py dox etg --nodoc sip build bdist_wheel; \
            ln -sf `readlink -f dist/wxPython*.whl` ../; \
          fi;
	# Finally touch a dummy file to avoid phony target
	@touch _build_wx/wxPython.target


wx: _build_wx/wxPython.target
	# Install the wheel
	@$(INSIDE_VENV) cd _build_wx && pip install wxPython*.whl


doc: VENV_PATH = ./_venv
doc: venv
	@echo "Building pdf and html documentation..."
	@$(INSIDE_VENV) pip install sphinx
	# @+$(INSIDE_VENV) $(MAKE) -C doc/ latexpdf
	@+$(INSIDE_VENV) $(MAKE) -C doc/ html


clean::
	@+$(MAKE) -C src/ clean


help:
	@echo "Top Makefile help message"

