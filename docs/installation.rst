Getting Started
**********

Overview of Requirements
============

* **python>=3.7**

  * numpy (any version except 1.10.*)
  * scipy
  * astropy (3.2.3)
  * h5py
  * scikit-learn
  * six
  * p_tqdm


* **compilers**
  
  * gcc
  * gfortran
  

* **Additional Packages (with Instructions Below)**
  
  * git  <http://git-scm.com/>
  * powderday <https://github.com/dnarayanan/powderday.git>
  * yt <http://yt-project.org>
  * FSPS <https://code.google.com/p/fsps/source/checkout>
  * python-fsps <https://dfm.io/python-fsps/current/>
  * Hyperion <http://www.hyperion-rt.org/>
  * Hyperion Dust Files <http://docs.hyperion-rt.org/en/stable/dust/dust.html>

Installation
============
    


Manual Installation
--------------

What follows is a self-contained installation manual, though for
problematic installs of any of the sub packages, it's definitely best
to visit the main docs on the main software site (which are always
linked below in each subsection).

.. _python:

python
--------------

`powderday <https://github.com/dnarayanan/powderday.git>`_ should work with python >=3.5 though is ideal with 3.6 (and some issues have been noted that may relate to python 3.7).
Please file an issue if you encounter one.

We very strongly recommend that the user set up a new python environment for the
`powderday <https://github.com/dnarayanan/powderday.git>`_ installation to avoid software conflicts.   This could look something like (assuming a ``conda`` installation of python)::

  >conda create --name pd_environment
  >source activate pd_environment

(And then when you want to exit the environment, you can type)::

  >source deactivate pd_environment

Then, whenever you're in the ``pd_environment``, everything you
install will remain contained to that particular installation of
python, and not conflict with any of your other installed packages.

.. _powderday:


powderday
--------------

Simply clone the latest and greatest from the repo::

  >git clone https://github.com/dnarayanan/powderday.git

To install, `cd` into the cloned repository and run the usual::

  >python setup.py install


.. _yt:



.. _Hyperion:

Hyperion
--------------

`Hyperion <http://www.hyperion-rt.org>`_ is the main work horse of
`powderday <https://github.com/dnarayanan/powderday.git>`_.  The full
directions for installation are well-described on the main
`Installation page for Hyperion
<http://docs.hyperion-rt.org/en/stable/installation/installation.html>`_.
Here, we summarize the installation which should get most users
through without any real difficulty.

There are two ways to install `Hyperion <http://www.hyperion-rt.org>`_.  The first is via ``conda``::

  >conda install -c conda-forge hyperion

Please note, though, that there is an issue with six no longer being
bundled with astropy that was fixed here:
https://github.com/hyperion-rt/hyperion/issues/219.  This said, at the
time of the last update of these docs (July 10th, 2020), this has not translated to the conda installation, meaning you will need to manually update all of the files listed here:

https://github.com/hyperion-rt/hyperion/issues/219#issuecomment-600036854  by replacing::

  >#from astropy.extern import six
  >import six

(for example, the files might be located in a location like:)::

  >home/desika.narayanan/miniconda3/envs/pd_test/lib/python3.6/site-packages/hyperion/filter/filter.py
  
The second and manual way to install `Hyperion
<http://www.hyperion-rt.org>`_ follows (note, for the manual installation you don't have to worry about the six replacement above):


#. First clone the main repository.::

     >git clone https://github.com/hyperion-rt/hyperion.git

#. Make sure that you have the correct modules loaded on your cluster.
   This will require a compiler, openmpi and HDF5.  For example, on
   the University of Florida HiPerGator supercomputing system, I would
   have::

   >module load intel/2018.1.163 openmpi/4.0.3 hdf5/1.10.1
     
#. Install the python module::

   >cd hyperion
   >python setup.py install


#. Ensure that if you type::

   >hyperion

it returns a sensible output.  It should return something along the lines of::

  >usage: hyperion [-h] [-f] [-m n_cores] input output
  >hyperion: error: too few arguments

If it can't find `Hyperion <http://www.hyperion-rt.org>`_, check the
the path that is near one of the last lines of the setup.py
installation (that is something associated with the number 755) and
make sure it's in your path.  It's most likely to be a python binaries
directory.

#. Install the submodules manually::

   >git submodule init
   >git submodule update

#. Install the Fortran binaries::

     > ./configure

or::

  >./configure --prefix=$HOME/local

or some such path if you aren't administrator on your computer.  Note
for this step you'll need your compilers, MPI and HDF5 installations
active (so, on a supercomputer you might need to load these modules
such as [for example, on the University of Florida HiPerGator
supercomputer])::

  >module load  1) intel/2018.1.163   2) openmpi/4.0.3   3) hdf5/1.10.1

of course please be careful of mixing and matching compilers, and
ensuring that you have the same compilers loaded for all
installations.
  
#. Compile the code::

   > make
   > make install
   

Note this will take a while!  Make sure this works by typing at the command line::

  >hyperion_sph

which should return something like::

  >Usage: hyperion_sph [-f] input_file output_file


  .. _Hyperion_dust:

Hyperion Dust Files
--------------

Unless you've written your own dust files, you will likely want to use
the pre-compiled dust files developed by Tom Robitaille (though don't
ship with `Hyperion <http://www.hyperion-rt.org>`_ due to their size).
To install these download them here:
http://docs.hyperion-rt.org/en/stable/dust/dust.html.  Then to
install::

  >tar -xvzf hyperion-dust-xxx.tar.gz
  >cd hyperion-dust-0.1.0
  >python setup.py build_dust

If you want to use the PAH model in `powderday
<https://github.com/dnarayanan/powderday.git>`_, you'll additionally need
these files in the same dust directory.  To download, click on the link,
then click 'raw' on the right side of each page.

1. https://github.com/hyperion-rt/paper-galaxy-rt-model/blob/master/dust/big.hdf5
2. https://github.com/hyperion-rt/paper-galaxy-rt-model/blob/master/dust/vsg.hdf5
3. https://github.com/hyperion-rt/paper-galaxy-rt-model/blob/master/dust/usg.hdf5

Please note the caveat that the PAH files are generated using some
approxmations described in `Robitaille et
al. <http://www.aanda.org/articles/aa/abs/2012/09/aa19073-12/aa19073-12.html>`_,
and we encourage the user of these PAH files to read this paper,
especially section 3.4.2.


yt
--------------------

Next we need `yt <http://yt-project.org>`_ - to install this, clone the source and install::

  >git clone https://github.com/yt-project/yt
  >cd yt
  >pip install -e .

Note, it is important to install this *after*  `Hyperion <http://www.hyperion-rt.org>`_.  This is because  if you used the conda installation of `Hyperion <http://www.hyperion-rt.org>`_ , then `yt <http://yt-project.org>`_ 3.x ships with it and auto-installs. However, powderday is no longer compatible with `yt <http://yt-project.org>`_ 3.x.



.. _fsps:

fsps
--------------

`fsps <https://code.google.com/p/fsps/source/checkout>`_ can be checked out with::
  
  > git clone https://github.com/cconroy20/fsps

and directions to the installation are in the `Manual <https://www.cfa.harvard.edu/~cconroy/ FSPS_files/MANUAL.pdf>`_.

To explicitly compile::

  make clean
  make
  
Finally, the SPS_HOME variable must be set in your environment to point to the FSPS/src directory.  For example, if your environment is bash, in your .bashrc set something along the lines of::
   
  >export SPS_HOME=/Users/desika/fsps/

Note that the same compilers used for `Hyperion
<http://www.hyperion-rt.org>`_ and `yt <http://yt-project.org>`_ need
to be used here.  An easy way to do this is in the Makefile to set F90=$(FC)


.. _python-fsps:

python-fsps
--------------

To install::

  >git clone --recursive https://github.com/dfm/python-fsps.git
  >cd python-fsps
  >python setup.py install

`python-fsps <https://dfm.io/python-fsps/current/>`_  will be installed automatically by the `powderday` setup.py script.
  
You can test the installation by opening python and typing::

>import fsps







Troubleshooting your Installation
============

  .. _fsps installation issues:

fsps Installation Issues
---------------
* One possibility can be that there are issues in compiling
   src/autosps.f90.  One solution is to replace RETURN with STOP in
   line 21.



  .. _python-fsps installation issues:

python-fsps installation issues
--------------
* With intel compilers (e.g., on the University of Florida HiPerGator system) you should try::
     
   >CC=icc F90=ifort python setup.py install

*  `python-fsps <https://dfm.io/python-fsps/current/>`_ can't find f2py
   
   f2py is a numpy package that is sometimes named f2py2.7 by numpy.
   At the same time, `python-fsps
   <https://dfm.io/python-fsps/current/>`_ expects it
   to be called f2py (as it sometimes is; for example in Anaconda).
   So, you might need to locate f2py (it ships with `yt
   <http://yt-project.org>`_, so if you for example use the `yt
   <http://yt-project.org>`_ python) you need to link the following
   files::

   >cd /Users/desika/yt-x86_64/bin
   >ln -s f2py2.7 f2py

   and::

   >cd /Users/desika/yt-x86_64/lib/python2.7/site-packages
   >ln -s numpy/f2py/ f2py

   This should hopefully fix it.


* Issues with 'f2py' in the  `python-fsps
   <https://dfm.io/python-fsps/current/>`_ installation:

   Numpy has made some changes to f2py in the 1.10.x version of numpy.
   The easiest fix is to use a non 1.10.* version of numpy (thanks to
   Ben Johnson for finding this).

*  `python-fsps <https://dfm.io/python-fsps/current/>`_ has mysterious
installation failures.  Often this has to do with a bad `FSPS
<https://github.com/cconroy20/fsps>`_ compilation. Even if it seems
like `FSPS <https://github.com/cconroy20/fsps>`_ has compiled, it may
not actually execute properly if the correct compilers aren't set in
the MakeFile.  Thanks to Ena Choi for pointing this one out.

  .. _hyperion installation issues:


Hyperion Installation Issues
---------------

  .. _yt installation issues:

   
yt Installation Issues
---------------

* If you have trouble with this installation, you may want to unload
your openmpi module that you previously had loaded for the `Hyperion
<http://www.hyperion-rt.org>`_ install.



* Another common trick to help the installation is to install with::

   >LDSHARED="icc -shared" CC=icc pip install -e .


System Specific Installation Notes
============

HiPerGator at the University of Florida
--------------

[1] The first set of instructions for the University of Florida
HiPerGator3.0 facility is to employ intel compilers, and to compile
everything manually.  This allows the greatest flexibility, as well as
the ability to use private forks of individual codes.

First, load up the compilers that we'll use throughout::

  >module load intel/2018.1.163
  >module load openmpi/4.0.3
  >module load hdf5/1.10.1
  >module load git

yt::

  >cd $HOME
  >git clone https://github.com/yt-project/yt
  >cd yt
  >pip install -e .



fsps and python-fsps

The development version of python-fsps now includes the Fortran FSPS source code::

  >cd $HOME
  >git clone --recursive https://github.com/dfm/python-fsps.git

then in your .bashrc set the analog to::
  
  >export SPS_HOME=$HOME/python-fsps/src/fsps/libfsps
  
  >cd python-fsps
  >CC=icc F90=ifort python setup.py install



hyperion::

  >cd $HOME
  >git clone https://github.com/hyperion-rt/hyperion.git
  >cd hyperion
  >pip install .
  >git submodule init
  >git submodule update

  >./configure --prefix=$HOME/local

  >make
  >make install

hyperion dust::

  >cd $HOME
  >wget http://pypi.python.org/packages/source/h/hyperion-dust/hyperion-dust-0.1.0.tar.gz
  >tar -xzvf hyperion-dust-0.1.0.tar.gz
  >cd hyperion-dust-0.1.0
  >python setup.py build_dust

  
powderday::

  >git clone https://github.com/dnarayanan/powderday.git
  >cd powderday
  >conda install numpy scipy cython h5py matplotlib psutil joblib six astropy scikit-learn ipython
  >python setup.py install

[2] The second set of instructions use gcc, but a manual installation of everything. Thanks to Prerak Garg for these.::

First, load up the compilers that we'll use throughout::

  >module load gcc/9.3.0 openmpi/4.1.1 libz/1.2.11 hdf5/1.10.1 git/2.30.1

  
yt::

  >cd $HOME
  >git clone https://github.com/yt-project/yt
  >cd yt
  >pip install -e .


fsps and python-fsps

The development version of python-fsps now includes the Fortran FSPS source code::

  >cd $HOME
  >git clone --recursive https://github.com/dfm/python-fsps.git


then in your .bashrc set the analog to::
  
  >export SPS_HOME=$HOME/python-fsps/src/fsps/libfsps
  
  >cd python-fsps
  >CC=gcc F90=gfortran F77=gfortran python setup.py install



hyperion::

  >cd $HOME
  >git clone https://github.com/hyperion-rt/hyperion.git
  >cd hyperion
  >pip install .
  >git submodule init
  >git submodule update

  >./configure --prefix=$HOME/local

  >make
  >make install

hyperion dust::

  >cd $HOME
  >wget http://pypi.python.org/packages/source/h/hyperion-dust/hyperion-dust-0.1.0.tar.gz
  >tar -xzvf hyperion-dust-0.1.0.tar.gz
  >cd hyperion-dust-0.1.0
  >python setup.py build_dust

  
powderday::

  >git clone https://github.com/dnarayanan/powderday.git
  >conda install numpy scipy cython h5py matplotlib psutil joblib six astropy scikit-learn ipython
  >cd powderday
  >python setup.py install




  
  

[3] The third set of instructions use gcc, and the conda installation
of `Hyperion <http://www.hyperion-rt.org>`_.  Thanks to Paul Torrey
for these.::

  >module load openmpi/4.1.1 libz/1.2.11 hdf5/1.10.1 conda/4.12.0 git/2.30.1 gcc
  >conda install -c conda-forge hyperion
  >python -c "import hyperion" (just to ensure no errors thrown)
  >hyperion (just to ensure command is found)
  >python -m pip install fsps
  >[set $SPS_HOME variable in .bashrc)
  >cd $HOME
  >git clone https://github.com/dnarayanan/powderday.git
  >cd powderday
  >python setup.py install

then fix import six line in the equivalent of all of these::

  >vi /home/paul.torrey/.conda/envs/pd_gcc/lib/python3.8/site-packages/hyperion/model/model.py
  >vi /home/paul.torrey/.conda/envs/pd_gcc/lib/python3.8/site-packages/hyperion/util/validator.py 
  >vi /home/paul.torrey/.conda/envs/pd_gcc/lib/python3.8/site-packages/hyperion/conf/conf_files.py
  >vi /home/paul.torrey/.conda/envs/pd_gcc/lib/python3.8/site-packages/hyperion/filter/filter.py
  >vi /home/paul.torrey/.conda/envs/pd_gcc/lib/python3.8/site-packages/hyperion/dust/dust_type.py
  >vi /home/paul.torrey/.conda/envs/pd_gcc/lib/python3.8/site-packages/hyperion/model/model_output.py
  >vi /home/paul.torrey/.conda/envs/pd_gcc/lib/python3.8/site-packages/hyperion/densities/flared_disk.py
  >vi /home/paul.torrey/.conda/envs/pd_gcc/lib/python3.8/site-packages/hyperion/densities/alpha_disk.py
  >vi /home/paul.torrey/.conda/envs/pd_gcc/lib/python3.8/site-packages/hyperion/densities/bipolar_cavity.py
  >vi /home/paul.torrey/.conda/envs/pd_gcc/lib/python3.8/site-packages/hyperion/densities/ulrich_envelope.py
  >vi /home/paul.torrey/.conda/envs/pd_gcc/lib/python3.8/site-packages/hyperion/densities/power_law_envelope.py 
  >vi /home/paul.torrey/.conda/envs/pd_gcc/lib/python3.8/site-packages/hyperion/densities/ambient_medium.py
  >vi /home/paul.torrey/.conda/envs/pd_gcc/lib/python3.8/site-packages/hyperion/model/sed.py
  >vi /home/paul.torrey/.conda/envs/pd_gcc/lib/python3.8/site-packages/hyperion/model/image.py
  >vi /home/paul.torrey/.conda/envs/pd_gcc/lib/python3.8/site-packages/hyperion/grid/yt3_wrappers.py
