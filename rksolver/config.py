import sys, sysconfig
import copy
import numpy
from distutils.extension import Extension
from distutils.util import get_platform
from distutils.dist import Distribution
from distutils.command.install_lib import install_lib


def get_rksolver_object_name():
  """Returns the name of the main rksolver shared library object"""

  ext_suffix = sysconfig.get_config_var('SOABI')

  if ext_suffix is None:
    filename = '_rksolver.so'
  else:
    filename = '_rksolver.{0}.so'.format(ext_suffix)

  return filename


def get_shared_object_path():
  """Returns the name of the distutils build directory"""

  install_lib_command = install_lib(Distribution())
  install_lib_command.initialize_options()
  install_lib_command.finalize_options()

  directory = install_lib_command.build_dir

  return directory


def get_rksolver():
  """Returns the path and name of the main shared library object"""

  return get_shared_object_path() + '/' + get_rksolver_object_name()



class configuration:
  """User-defined build configuration options for OpenRK
  Configuration options may be set using compile time flags. To view a
  list of these options, run 'python setup.py install --help' in the
  console. The default configuration options are shown below and should
  only be revised by developers familiar with the code and its configuration
  management system.
  """

  #############################################################################
  #                               User Options
  #############################################################################

  # Only supports GCC as the default compiler right now ??????
  # Default C++ compiler for the main openmoc module is GCC
  cc = 'g++'

  # Compile using ccache (for developers needing fast recompilation)
  with_ccache = False
  
  # Supported C++ compilers: 'gcc'
  cpp_compilers = list()

  # List of C/C++/CUDA distutils.extension objects which are created based
  # on which flags are specified at compile time.
  extensions = list()

  # List of the packages to install - only openmoc is guaranteed to be built
  # while the others will be built based on which flags are specified
  # at compile time
  packages = ['rksolver']

  #############################################################################
  #                                 Source Code
  #############################################################################

  # Dictionary of source code files to compile for each extension module
  sources = dict()

  sources['gcc'] = ['rksolver/rksolver_wrap.cpp',
                    'src/diffSolver.cpp',
                    'src/transientSolver.cpp',
                    'src/pkeSolver.cpp',
                    'src/ftSolver.cpp',
                    'src/Sparse.cpp',
                    'src/XSdata.cpp',
                    'src/Mesh.cpp',
                    'src/utils.cpp',
                    'src/Solutions.cpp']


  #############################################################################
  #                                Compiler Flags
  #############################################################################

  # A dictionary of the compiler flags to use for each compiler type
  compiler_flags = dict()

  compiler_flags['gcc'] = ['-c', '-O3', '-ffast-math', '-fopenmp',
                           '-std=c++0x', '-fpic']


  #############################################################################
  #                                 Linker Flags
  #############################################################################

  # A dictionary of the linker flags to use for each compiler type
  linker_flags = dict()

  if (get_platform()[:6] == 'macosx'):
    linker_flags['gcc'] = ['-fopenmp', '-dynamiclib', '-lpython2.7',
                           '-Wl,-install_name,' + get_rksolver_object_name()]
  else:
    linker_flags['gcc'] = ['-fopenmp', '-shared',
                           '-Wl,-soname,' + get_rksolver_object_name()]


  #############################################################################
  #                               Shared Libraries
  #############################################################################

  # A dictionary of the shared libraries to use for each compiler type
  shared_libraries = dict()

  shared_libraries['gcc'] = ['stdc++', 'gomp', 'dl','pthread', 'm']


  #############################################################################
  #                              Library Directories
  #############################################################################

  # A dictionary of the library directories to use for each compiler type
  # if not set in the LD_LIBRARY_PATH environment variable
  library_directories = dict()

  usr_lib = sys.exec_prefix + '/lib'

  library_directories['gcc'] = [usr_lib]


  #############################################################################
  #                              Include Directories
  #############################################################################

  # A dictionary of the include directories to use for each compiler type
  # for header files not found from paths set in the user's environment
  include_directories = dict()

  include_directories['gcc'] = list()

  ###########################################################################
  #                                 SWIG Flags
  ###########################################################################

  # A list of the flags for SWIG
  swig_flags = ['-c++', '-python', '-keyword'] #, '-keyword', '-py3']


  def setup_extension_modules(self):
    """Sets up the C/C++/CUDA extension modules for this distribution.
    Create list of extensions for Python modules within the openmoc
    Python package based on the user-defined flags defined at compile time.
    """

    try:
      numpy_include = numpy.get_include()

    except AttributeError:
      numpy_include = numpy.get_numpy_include()

    # Add the NumPy include directory to the include directories
    # list for each type of compiler
    for cc in self.include_directories.keys():
      self.include_directories[cc].append(numpy_include)

    # The main openmoc extension (defaults are gcc and single precision)
    self.extensions.append(
      Extension(name = '_rksolver',
                sources = copy.deepcopy(self.sources[self.cc]),
                library_dirs = self.library_directories[self.cc],
                libraries = self.shared_libraries[self.cc],
                extra_link_args = self.linker_flags[self.cc],
                include_dirs = self.include_directories[self.cc],
                swig_opts = self.swig_flags + ['-D' + self.cc.upper()]))
