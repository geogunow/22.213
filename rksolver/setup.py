from distutils.core import setup
from distutils.command.build_ext import build_ext
from distutils.command.build_py import build_py
from distutils.command.install import install
from distutils.errors import DistutilsOptionError
import os, string
import config

# Instantiate a configurations class for this OpenMOC build
config = config.configuration()

class custom_install(install):
  """Defines the compile time options for OpenRK.
  This class derives from the distutils.command.install class. Distutils
  provides a set of flags which may be invoked by the user at compile
  time. The custom_install class adds to that list a series of options
  which are useful in compiling a specific distribution of OpenRK.
  To view a list of available compile time options, simply type the following
  from a console:
  python setup.py install --help
  The custom_install class extracts user-defined options from the command
  line and uses them to populate the config.configurations class. The
  configurations class then uses these options to generate a list of
  Python C/C++/CUDA extension objects which are delivered to the distutils
  setup method at the end of this script.
  Developers wishing to extend OpenRK's functionality with new compilation
  options will need to extend this class and the configurations class.
  """

  # The user options for a customized OpenMOC build
  user_options = [
    ('cc=', None, "Compiler (gcc) for main openmoc module"),
    ('with-gcc', None, "Build openmoc.gnu modules using GNU compiler"),
    ('with-ccache', None, "Build with ccache for rapid recompilation")
  ]

  # Include all of the default options provided by distutils for the
  # install command parent class
  user_options += install.user_options

  # Set some compile options to be boolean switches
  boolean_options = ['with-gcc',
                     'with-ccache']


  # Include all of the boolean options provided by distutils for the
  # install command parent class
  boolean_options += install.boolean_options


  def initialize_options(self):
    """Set the default OpenRK build options
    The default installation is invoked by following console command:
    python setup.py install
    This will build the main openmoc C/C++ Python extension using the
    GCC compiler with single precision. No additional modules will be
    build with Intel or IBM compilers, or with double precision.
    """

    # Run the install command parent class' initialize_options method
    install.initialize_options(self)

    # Default compiler and precision level for the main openmoc module
    self.cc = 'gcc'
    
    # By default, do not build openmoc.gnu.single, openmoc.intel.double, etc
    # extension modules
    self.with_gcc = False
    self.with_ccache = False
    
  def finalize_options(self):
    """Extract options from the flags invoked by the user at compile time.
    This method performs error checking of the options specified by
    the user at compile time, and initialize the config.configurations
    class instance. The method conclude with a call to the
    configurations.setup_extension_modules class method which builds
    a list of C/C++/CUDA extension modules to be passed to the distutils
    setup method at the end of this script.
    """

    # Run the install command parent class' finalize_options method
    install.finalize_options(self)

    config.with_ccache = self.with_ccache
    
    # Check that the user specified a supported C++ compiler
    if self.cc not in ['gcc']:
      raise DistutilsOptionError \
            ('Must supply the -cc flag with one of the supported ' +
             'C++ compilers: gcc')
    else:
      config.cc = self.cc

    # Build the openmoc.gnu.single and/or openmoc.gnu.double
    # extension module(s)
    if self.with_gcc:
      config.cpp_compilers += ['gcc']

    # Build a list of the C/C++/CUDA extension modules to be built
    # for this distribution
    config.setup_extension_modules()

def customize_compiler(self):
  """Inject redefined _compile method into distutils
  This method enables us to choose compilers based on the macros defined
  in the compiler flags (ie, '-DGNU', '-DCUDA', etc), or on the
  source extension (ie, *.cpp, *.cu, etc.).
  Adapted from Robert McGibbon's CUDA distutils setup provided in open source
  form here: https://github.com/rmcgibbo/npcuda-example
  """

  # Save reference to the default _compile method
  super_compile = self._compile

  # Redefine the _compile method. This gets executed for each
  # object but distutils doesn't have the ability to change compilers
  # based on source extension, so we add that functionality here
  def _compile(obj, src, ext, cc_args, extra_postargs, pp_opts):

    self.set_executable('compiler_so', 'ccache gcc')

    postargs = config.compiler_flags['gcc']

    # Now call distutils-defined _compile method
    super_compile(obj, src, ext, cc_args, postargs, pp_opts)

  # Inject our redefined _compile method into the class
  self._compile = _compile

def customize_linker(self):
  """Inject redefined link method into distutils
  This method enables us to choose the linker based on the name
  of the output shared library filename (ie, _openmoc_intel_single.so)
  Adapted from Robert McGibbon's CUDA distutils setup provided in open source
  form here: https://github.com/rmcgibbo/npcuda-example
  """

  # Save references to the default link method
  super_link = self.link

  # Redefine the link method. This gets executed to link each extension
  # module. We add the functionality to choose the compiler for linking
  # based on the name of the extension module
  def link(target_desc, objects, output_filename,
           output_dir=None, libraries=None,
           library_dirs=None, runtime_library_dirs=None,
           export_symbols=None, debug=0, extra_preargs=None,
           extra_postargs=None, build_temp=None, target_lang=None):

    # If compiling different extensions of openmoc using different compilers
    # and/or floating point precision levels, we must remove autogenerated
    # files from distutils for each subsequent extension. If the user is
    # compiling multiple modules at once (ie, openmoc.gnu.single and
    # openmoc.intel.single) we have to ensure that only one of the objects
    # openmoc_gnu_single.o or openmoc_intel_single.o is specified at the
    # link stage. Unfortunately, distutils enumerates all object files
    # compiled up to this point at the final linker stage which leads to
    # 'previously defined...' errors
    for obj in objects[:]:

      if 'rksolver' in obj:

        if 'gnu' in output_filename and 'gnu' not in obj:
          objects = [o for o in objects if o is not obj]
        elif 'gnu' not in output_filename and 'gnu' in obj:
          objects = [o for o in objects if o is not obj]

    # If the linker receives -fopenmp as an option, then the objects
    # are built by a GNU compiler
    if '-fopenmp' in extra_postargs:
      self.set_executable('linker_so', 'g++')
      self.set_executable('linker_exe', 'g++')

    # Now call distutils-defined link method
    super_link(target_desc, objects,
               output_filename, output_dir, libraries,
               library_dirs, runtime_library_dirs,
               export_symbols, debug, extra_preargs,
               extra_postargs, build_temp)

  # Inject our redefined link method into the class
  self.link = link


# Run the customize_compiler to inject redefined and customized _compile and
# link methods into distutils
class custom_build_ext(build_ext):
  """Customizes distutils to work with different compiler types
  This class derives from the distutils.command.build_ext command class.
  It extends build_ex by creates customized compile and link methods
  which can accommodate different compiler types and options.
  """

  def build_extensions(self):
    customize_compiler(self.compiler)
    customize_linker(self.compiler)

    os.system('swig {0} -o '.format(str.join(' ', config.swig_flags)) + \
                'rksolver/rksolver_wrap.cpp rksolver/rksolver.i')

    build_ext.build_extensions(self)


# Run the distutils setup method for the complete build
dist = setup(name = 'rksolver',
      version = '0.1',

      # Set the C/C++/CUDA extension modules built in setup_extension_modules()
      # in config.py based on the user-defined flags at compile time
      ext_modules = config.extensions,

      # Extract all of the Python packages for OpenRK
      # (ie, openmoc.log, openmoc.materialize, etc)
      packages = config.packages,

      # Inject our custom compiler and linker triggers
      cmdclass={ 'build_ext': custom_build_ext,
                 'install': custom_install},

)

# Rerun the build_py to setup links for C++ extension modules created by SWIG
# This prevents us from having to install twice
build_py = build_py(dist)
build_py.ensure_finalized()
build_py.run()
