PYTHON   = python
PYMAJ    = 3
PYMIN    = 5

FC       = gfortran
F2PY     = f2py --f77exec=$(FC) --f90exec=$(FC)

NO_VENV  = 0
DEBUG    = 0
VERBOSE  = 0
# Absolute path: f2py's meson backend builds in a temporary directory, where a
# relative -I./build does not resolve.
BUILDDIR ?= $(abspath build)
INSTALL_PREFIX = $(HOME)/.local


################################################################################
# GFORTRAN CONFIGURATION                                                       #
################################################################################
GFORTRAN_FFLAGS      = -O2 -ffast-math
GFORTRAN_FFLAGS_DBG  = -g -Wall -Wextra -Warray-temporaries -Wconversion 
GFORTRAN_FFLAGS_DBG += -fbacktrace -ffree-line-length-0 -fcheck=all
GFORTRAN_FFLAGS_DBG += -ffpe-trap=zero,overflow,underflow -finit-real=nan
################################################################################

################################################################################
# IFORT CONFIGURATION                                                          #
################################################################################
IFORT_FFLAGS     = 
IFORT_FFLAGS_DBG =
################################################################################

################################################################################
# F2PY CONFIGURATION                                                           #
################################################################################
F2PYFLAGS     = --opt=-O2 -llapack
F2PYFLAGS_DBG = --debug-capi --debug
################################################################################



################################################################################
# /!\ DO NOT EDIT BELOW THAT LINE (unlesss you know what you're doing...)      #
# CORE CONFIGURATION                                                           #
################################################################################
VERSION:=$(shell git describe|sed 's/-\([0-9]\+\)-.*/.dev\1/g')
VENV_PATH := $(INSTALL_PREFIX)/src/msspec_venv_$(VERSION)


# Checking Python path and version
PYTHON_EXE := $(shell command -v $(PYTHON) 2> /dev/null)

ifndef PYTHON_EXE
  $(error Unable to find the $(PYTHON) executable!)
endif

ifeq ($(shell $(PYTHON_EXE) -c "import sys; exit(sys.version_info >= ($(PYMAJ),$(PYMIN)))"; echo $$?),0)
   $(error Python version >= $(PYMAJ).$(PYMIN) is needed!)
endif

# Checking your distribution type and release
LSB_RELEASE := $(shell command -v lsb_release 2> /dev/null)
DISTRO_RELEASE =

ifdef LSB_RELEASE
  DISTRO_RELEASE = $(shell echo $(shell $(LSB_RELEASE) -s -i)-$(shell $(LSB_RELEASE) -s -r) | tr [:upper:] [:lower:] )
endif

ifndef DISTRO_RELEASE
  $(warning Unable to guess your OS name and version)
endif

# Define the INSIDE_VENV variable
ifeq ($(NO_VENV), 1)
  INSIDE_VENV =
else
  INSIDE_VENV = . $(VENV_PATH)/bin/activate && 
endif




ifeq ($(VERBOSE),0)
  OUPUT_REDIRECTION := 1>/dev/null 2>/dev/null
  MAKEFLAGS = -s --no-print-directory
else
  OUPUT_REDIRECTION:=
  MAKEFLAGS += --debug=b
endif

PREFIX=
SUFFIX=

ifeq ($(FC),gfortran)
  PREFIX = GFORTRAN
endif

ifeq ($(FC),IFORT)
  PREFIX = IFORT
endif

ifeq ($(DEBUG),1)
  SUFFIX = _DBG
endif

FFLAGS = $($(PREFIX)_FFLAGS$(SUFFIX))

OBJS = $(addprefix $(BUILDDIR)/, $(patsubst %.f,%.o, $(filter-out $(MAIN_F), $(SRCS))))

.PHONY: clean obj all info

help:


info:
	@echo "Compilation information for msspec $(VERSION)"
	@echo "FC    : $(FC)"
	@echo "FFLAGS:"
	@for flag in $(FFLAGS); do echo "\t$$flag"; done


all: $(SO)


obj: $(OBJS)


clean::
	@if test -d $(abspath $(BUILDDIR)); then \
	    echo "Removing $(abspath $(BUILDDIR))..."; \
	    rm -r $(abspath $(BUILDDIR)); \
	 fi
	@if test x$(SO) != x; then \
           if test -f $(SO); then \
             echo "Removing $(SO)..."; rm $(SO); \
           fi; \
         fi


$(BUILDDIR)/%.o: %.f
	@echo "Compiling $@..."
	mkdir -p $(basename $@)
	$(FC) $(FFLAGS) -J $(BUILDDIR) -I $(BUILDDIR) -fPIC -o $@ -c $^ $(OUPUT_REDIRECTION)


$(SO): $(OBJS) $(MAIN_F)
	@echo "building Python binding $@..."      
	mkdir -p $(BUILDDIR)
	$(F2PY) $(F2PYFLAGS) -I$(BUILDDIR) -m $(basename $@) -c $(OBJS) $(MAIN_F) $(OUPUT_REDIRECTION)
	mv $(basename $@).*.so $@
