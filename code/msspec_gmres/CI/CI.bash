#!/bin/bash

SCRIPT_NAME=$(basename $0)
TOPLEVEL=$(git rev-parse --show-toplevel)
VERSION=$(git describe|sed 's/-\([[:digit:]]\+\)-.*/\.post\1/')
(return 0 2>/dev/null) && sourced=1 || sourced=0

pythonpath_prepend() {
    if [ -d "$1" ]; then
        PYTHONPATH=${PYTHONPATH//":$1:"/:} #delete all instances in the middle
        PYTHONPATH=${PYTHONPATH/%":$1"/} #delete any instance at the end
        PYTHONPATH=${PYTHONPATH/#"$1:"/} #delete any instance at the beginning
        PYTHONPATH="$1${PYTHONPATH:+":$PYTHONPATH"}" #prepend $1 or if $PATH is empty set to $1
	export PYTHONPATH
    fi
}

# Setup or upgrade the Python virtual environment
setup_virtualenv(){
    VENV_PATH="$1"
    virtualenv --python=python3 --system-site-packages "$VENV_PATH"
    source "$VENV_PATH/bin/activate"
    pip install --upgrade pip
    pip install --upgrade numpy ase h5py lxml pint terminaltables pycairo 
    pip install --upgrade sphinx
    status=$?
    deactivate
    [ $status -eq 0 ] || exit 1
}

activate_msspec_python3_virtualenv() {
    VENV_PATH="$1"
    # do not activate virtualenv if already activated
    if [ x"${VIRTUAL_ENV}" == x ]; then
	# We are not in a virtualenv
        if [ x"${VENV_PATH}" != x ]; then
            # activate the specified venv
            export MSSPEC_ROOT="${TOPLEVEL}"
            pythonpath_prepend "$MSSPEC_ROOT/src"
	    #PYTHONPATH="$PYTHONPATH":"$MSSPEC_ROOT"/src
            source "$VENV_PATH/bin/activate"
        else
            # error
            echo "Error: No venv specified!!"
        fi
    else
        # We are in a virtualenv so load only variables
        export MSSPEC_ROOT="${TOPLEVEL}"
        pythonpath_prepend "$MSSPEC_ROOT/src"
    fi
}


make_pybinding() {
	echo "Creating Python library for Fortran code, please wait... "
	activate_msspec_python3_virtualenv "$1"
	cd ${TOPLEVEL}/src/ 
	make pybinding || exit 1
}


make_selfex() {
	echo "Creating the self-extractible setup program... "
	cd ${TOPLEVEL}
        make selfex || exit 1
        mv pymsspec-${VERSION}.setup package/
}

make_pdf() {
    DEST=${TOPLEVEL}/package
    mkdir -p ${DEST}
    # first activate the virtual environment
    activate_msspec_python3_virtualenv "$1"
    # then build the documentation.
    cd "${TOPLEVEL}/doc"
    make --no-print-directory latexpdf || exit 1
    # remove previous pdf
    rm -rf "${DEST}"/*.pdf
    # copy the produced pdf file to the final location
    cp "build/latex/MsSpec-python.pdf" "${TOPLEVEL}/package/MsSpec-${VERSION}.pdf"
}

make_html () {
    make_pdf "$1"
    # first activate the virtual environment
    activate_msspec_python3_virtualenv "$1"
    # then build the documentation.
    cd "${TOPLEVEL}/doc"
    make --no-print-directory html || exit 1
}

make_tests() {
    echo "generate results for unit tests..."
    activate_msspec_python3_virtualenv "$1"
    make_pybinding
    python -c "from msspec.tests import create_tests_results; create_tests_results()" || exit 1
}


create_package() {
    echo "Creating the package..."
    #make_tests "$1"
    make_selfex "$1"
    cleanup
}

cleanup() {
    echo "cleaning up..."
   # remove results.txt file
    rm -rf ${TOPLEVEL}/src/msspec/results.txt
}

show_help () {
    echo "Usage: $SCRIPT_NAME [OPTIONS]"
    echo "List of possible options:"
    echo "    -i   Initialize a Python virtual environnment in FOLDER"
    echo "    -d   Build the html documentation within the Python "
    echo "         virtual environnment in FOLDER"
}



if [ "$sourced" == 0 ]; then

	while getopts "hi:d:t:p:" option; do
	    case $option in
		i) setup_virtualenv "$OPTARG"
		   ;;
		d) make_html "$OPTARG"
		   ;;
		t) make_tests "$OPTARG"
		   ;;
		p) create_package "$OPTARG"
		   ;;
		h) show_help
		   ;;
		*) show_help
		   ;;
	    esac
	done

fi
