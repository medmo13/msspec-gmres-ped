Introduction
============

This is the Python MsSpec version with support for Python 3.x and dynamic memory allocation for Phagen and Spec


Installation
============

Requirements
------------

To compile MsSpec shared libraries, you need:

- Python >= 3.6
- gfortran compiler
- cairo library (runtime files and header files)

Depending on your OS and Python version, you may need all the requirements (runtime libraries
and header files) to build wxPython (please consult your distribution's package list for the 
appropriate package names) and their dependencies:

- python-dev (for your version of python)
- gtk (preferably version 3, but depends on your needs)
- gstreamer
- gstreamer-plugins-base
- glut
- libwebkitgtk (matching your gtk version)
- libjpeg
- libpng
- libtiff
- libsdl
- libnotify
- libsm

You also need a tool to create and manage Python virtual environments. For now only
*virtualenv* is supported. So if you plan to install msspec with the one-step command,
you need to install it.



Fetching the code
-----------------

You first need to clone the repository with the command below:

```Bash
git clone https://git.ipr.univ-rennes1.fr/epsi/msspec_python3.git
```

If you want to pull the devel branch (unstable):

```Bash
git clone --branch devel https://git.ipr.univ-rennes1.fr/epsi/msspec_python3.git
```

This will create a msspec\_python3 folder with a local copy of the repository.



One step install
----------------

The command below will do all the stuff to install msspec in its dedicated Python
virtual environment provided all the above requirements are fullfilled:


```Bash
cd msspec_python3
make install
```

This will:

- Create a Python virtual environment in your '~/.local/src' folder
- Install all the required Python packages in this virtual environment.
- Install (or even maybe build) wxPython for your OS and Python version in this virtual environment.
- Build the msspec shared librarires.
- Install the Python msspec package in this virtual environment.
- Install a frontend script called 'msspec' in the '~/.local/bin' folder.

You can tweak the process by specifying space separated  key=value pairs on the command line.
Allowed keys are:

- PYTHON, to give the executable name of the python interpreter to use rather than the default (python)
- VERBOSE, set to 1 if you want the install process to be more verbose
- NO\_VENV, set to 1 to not create a virtual environment
- VENV\_PATH, specify the path to a virtual environment root folder
- DEBUG, set to 1 to add debugging symbols
- BUILDDIR, to set a different name than the default to store building artifacts (./build)
- FC, the Fortran compiler to use (gfortran)
- F2PY

To install msspec in developpment mode, use:
```Bash
cd msspec_python3
make devel
```

This will do the same steps as above except that the virtual environment will be located in 
the "\_venv" folder inside the "msspec\_python3" folder and the frontend script will not be
installed.


Working with msspec
-------------------

Use the frontend command to process a Python script with msspec.

```Bash
msspec -p myscript.py
```

This will execute the myscript.py inside the msspec virtual environment.

You can also load a hdf5 file with the '-l' option:

```Bash
msspec -l results.hdf5
```

You can load an iPython console with "msspec -i". It is a good way to modify the virtual environment
of MsSpec since you can add or remove packages withe the "%pip" command.

If you want to activate the MsSpec virtual environment in your current shell session:

```Bash
eval $(msspec -e)
```

You can then run scripts using the 'python' command directly.
Type in 'deactivate' to go back to your standard shell.

See the online [tutorials](https://msspec.cnrs.fr/tutorials/index.html) for more information on how to use msspec.


Uninstall MsSpec
----------------

To remove MsSpec, use the command:

```Bash
msspec -u
```

This will simply remove the virtual environment folder in '~/.local/src' and the frontend script in '~/.local/bin'.


When Installed in "devel mode", nothing is created outside the msspec\_python3 folder, so you just have to remove 
that folder to get rid of msspec on your computer


Manual control over the install
-------------------------------

If for any reason you want to control each step of the building process (for example you do not have *virtualenv* but you
are using *conda*), here is how:

- Create a virtual environment or use an existing one and activate it
- Install the list of requirements that are listed in the src/pip.freeze file
- Install wxPython
- Build the msspec shared librarires with the command "make pybinding"
- Install the msspec Python package with "pip install -e src/"

The last step can be replaced by adding the src/ path to your PYTHONPATH, but bear in mind that the package will then
be also visible to other virtual environments no matter they fit the msspec requirements.
