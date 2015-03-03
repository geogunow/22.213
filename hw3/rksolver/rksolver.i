/* File: rksolver.i */
%module rksolver

%{
  #define SWIG_FILE_WITH_INIT
  #include "../src/diffSolver.h"
  #include "../src/utils.h"
  #include "../src/Mesh.h"
  #include "../src/XSdata.h"
  #include "../src/Sparse.h"
  #include "../src/Solutions.h"
  
  #define printf PySys_WriteStdout

  /* Exception helpers */
  static int swig_c_error_num = 0;
  static char swig_c_err_msg[1024];

  const char* err_occurred(void) {
    if (swig_c_error_num) {
      swig_c_error_num = 0;
      return (const char*)swig_c_err_msg;
    }
    return NULL;
  }

  void set_err(const char *msg) {
    swig_c_error_num = 1;
    strncpy(swig_c_err_msg, msg, 1024);
  }

%}

%exception {
  try {
    $function
  } catch (const std::exception &e) {
      SWIG_exception(SWIG_RuntimeError, e.what());
  }
}


%include "numpy.i"
%include "std_vector.i"

/* SWIG template for exporting vectors of doubles */
namespace std {
  %template(DoubleVector) vector<double>;
}

%init %{
  import_array();
%}


%apply (int* IN_ARRAY1, int DIM1) {(int * npts, int n_npts)};
%apply (double* IN_ARRAY1, int DIM1) {(double * widths, int n_widths), 
        (double * xs, int ng), (double * vals, int len), 
        (double * timeArray, int n_steps};

/* Typemap for setting transient object:
 * Transient.setTimes (double * timeArray, int n_interp,
 *          double * calcTimes, int n_steps)
 * Transient.setMeshVector ( Mesh * meshArray, int n_interp)
 * method - allows user to declare a transient from python
 */
%typemap(in) (Mesh ** meshArray, int n_interp) {
    $2 = PySequence_Length($input); // num mesh
    $1 = (Mesh**) malloc(($2) * sizeof(Mesh)); // mesh array
    
    /* loop through mesh */
    for (int i = 0; i < $2; i++) {
        PyObject* o = PyList_GetItem($input, i);
        void *p1 = 0;
        SWIG_ConvertPtr(o, &p1, SWIGTYPE_p_Mesh, 0 | 0);
        $1[i] = (Mesh *) p1;
    }
}

/* Typemap for Mesh::setMaterials (XSdata** materials, int n_materials)
 * method - allows user to pass in Python lists of XSdata for each node
 */
%typemap(in) (XSdata** materials, int n_materials) {

  if (!PyList_Check($input)) {
    PyErr_SetString(PyExc_ValueError,"Expected a Python list of integers "
                    "for the Lattice cells");
    return NULL;
  }

  $2 = PySequence_Length($input);  // num_materials
  $1 = (XSdata**) malloc(($2) * sizeof(XSdata*)); // XSdata

  /* Loop over materials */
  for (int i = 0; i < $2; i++) {

     /* Extract the value from the list at this location and convert
       * SWIG wrapper to pointer to underlying C++ class instance */
      PyObject* o = PyList_GetItem($input, i);
      void *p1 = 0;
      SWIG_ConvertPtr(o, &p1, SWIGTYPE_p_XSdata, 0 | 0);
      $1[i] = (XSdata*) p1;
  }
}

/* Typemap for XSdata::setSigs (double ** xs, int ng1, int ng2)
 * method - allows user to pass scattering cross sections
 */
%typemap(in) (double ** xs, int ng1, int ng2) {

  if (!PyList_Check($input)) {
    PyErr_SetString(PyExc_ValueError,"Expected a Python list of integers "
                    "for the Lattice cells");
    return NULL;
  }

  $3 = PySequence_Length(PyList_GetItem($input,0)); // ng2
  $2 = PySequence_Length($input);  // ng1
  $1 = (double**) malloc(($2) * sizeof(double*)); // xs

  /* Loop over materials */
  for (int i = 0; i < $2; i++) {

    /* Get the inner list in the nested list */
    PyObject* outer_list = PyList_GetItem($input, i);

    /* Check that the length of this list is the same as the length
     * of the first list */
    if (PySequence_Length(outer_list) != $3) {
      PyErr_SetString(PyExc_ValueError, "Size mismatch. Expected $2 x $3 "
                      "elements for Lattice\n");
        return NULL;
    }

    double * temp = (double *) malloc( ($3) * sizeof(double));
    
    for (int j=0; j < $3; j++){
    

     /* Extract the value from the list at this location and convert
       * SWIG wrapper to pointer to underlying C++ class instance */
      PyObject* o = PyList_GetItem(outer_list, j);
     
      /* If value is a number, cast it as an int and set input array value */
        if (PyNumber_Check(o)) {
          temp[j] = (double) PyFloat_AsDouble(o);
        }
        else {
          free($1);
          PyErr_SetString(PyExc_ValueError,"Expected a matrix of numbers "
                          "for cross-section values\n");
          return NULL;
        }
    }
    $1[i] = temp;
  }
}


%include <exception.i>
%include ../src/diffSolver.h
%include ../src/Mesh.h
%include ../src/XSdata.h
%include ../src/utils.h
%include ../src/Sparse.h
%include ../src/Solutions.h

#define printf PySys_WriteStdout
