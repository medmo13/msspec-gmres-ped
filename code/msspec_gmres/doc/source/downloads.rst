
.. |LINUXSCRIPT| replace:: https://git.ipr.univ-rennes1.fr/epsi/msspec_python3/raw/branch/devel/utils/dockerized/linux/msspec


##########################
Download and install notes
##########################


************
Installation
************


There are 2 ways to install MsSpec. You can either:
	- Use a Docker image. This is, by far, the most straightforward and easy way to work with MsSpec.
	- Compile your own version for your system.



1. Using a Docker image
-----------------------

You first need `Docker <https://www.docker.com>`_ to be installed on your system. 
This is really straightforward. Please see the `Docker documentation
<https://docs.docker.com/engine/install/>`_ for more information depending on 
your OS.


* **For Linux**, 
  
  * Download :download:`this script <../../utils/dockerized/linux/msspec>` and
    make it executable (with the *chmod u+x msspec* command)
 
  * Place it in a location known to your $PATH (or add this location to your $PATH if needed).
    The *~/.local/bin* folder is a good choice.

* **For Windows**

  * You need a running X server. Download and install `VcXsrv <https://sourceforge.net/projects/vcxsrv/>`_

  * Install the small :download:`MsSpec frontend <../../utils/dockerized/windows/msspec_setup.exe>`

  * While not necessary, it is also a good idea to use a better terminal application than the default one. 
    `Windows Terminal <https://www.microsoft.com/fr-fr/p/windows-terminal/9n0dx20hk701#activetab=pivot:overviewtab>`_ 
    may be a good choice.

* **For Mac**
 
  * *To be documented*


Open a terminal in Linux, or a powershell in Windows and type in::

  msspec

The first time you run the command, it will download the msspec image (it may takes several minutes or half an hour 
depending on your internet connexion).
The command will automatically create a new container and start it.
As the command was entered without any argument on the command-line, the help message should be printed on the screen

.. code-block:: text

	Usage: 1) msspec -p [PYTHON OPTIONS] SCRIPT [ARGUMENTS...]
	       2) msspec [-l FILE | -i | -h]
	       3) msspec [bash | reset]

	Form (1) is used to launch a script
	Form (2) is used to load a hdf5 data file
	Form (3) is used to control the Docker container/image.

	List of possible options:
	    -p   Pass every arguments after this option to the msspec
		 virtual environment Python interpreter.
	    -i   Run the interactive Python interpreter within msspec
		 virtual environment.
	    -l   Load and display a *.hdf5 data file in a graphical
		 window.
	    -v   Print the version.
	    -h   Show this help message.

	    bash        This command starts an interactive bash shell in the
			MsSpec container.
	    reset       This command removes the MsSpec container (but not the
			image). Changes made in the container will be lost and
			any new call to msspec will recreate a new fresh container.


2. Compile your own version
---------------------------

To install MsSpec this way, follow the instructions `here <https://git.ipr.univ-rennes1.fr/epsi/msspec_python3>`_


***************************
Running your Python scripts
***************************

You can run your MsSpec Python scripts (e.g. *my_script.py*) by typing in::

	msspec -p my_script.py

This command is equivalent to activating the Python virtual environment MsSpec is intsalled in 
and to run the Python interpreter of that environment (everything that follows the -p option is
passed to the python command).

You can also launch the Interactive Python (iPython) of MsSpec::

	msspec -i

Inside this interactive session, you can run a script with the command::

	%run my_script.py

You can interact with your filesystem with the classical *cd*, *ls*, *cp*, *rm*... commands.
and you can edit your script with::

	%ed my_script.py


.. warning::

	**If using the Docker image of MsSpec in Linux**, your home folder on the host machine is bind mounted
        in the same location in the container and your UID and GID are also set so that creating files within
        your home file hierarchy is totally transparent.

	**If using the Docker image of MsSpec in Windows**, the drive containing your "My Documents" folder on the
        host machine is bind mounted on the container at the root of the filesystem. For example, if your 
        "My documents" folder on Windows are in 'D:\\Home\\Bob\\MyDocuments', it will be available in the container as
        '/D/Home/Bob/MyDocuments'. It has two consequences:

	#. The msspec command will fail if running on another drive than the one where is located "My Documents"

	#. You have to specify filenames with the Unix slashes. For example if you want to run the script located
           in *.\\results\\my_script.py*, you will have to enter *msspec -p ./results/my_script.py*


**************
Uninstallation
**************

* **Under Linux**, type in::

	msspec -u

* **Under Windows**, simply `uninstall the application from the Settings page
  <https://support.microsoft.com/en-us/windows/uninstall-or-remove-apps-and-programs-in-windows-10-4b55f974-2cc6-2d2b-d092-5905080eaf98>`_.
  A command window will pop-up and you will have to answer 'y' to remove MsSpec.

* **Under Mac OS**, *to be documented.*
