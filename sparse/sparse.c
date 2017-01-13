
//#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include "Python.h"
#include "numpy/arrayobject.h"
#include "sparse.h"


// Convert spherical coordinates to Cartesian coordinates.
inline double sph_to_cart_x(double k, double theta_k, double phi_k){
    return k * sin(theta_k) * cos(phi_k);
}
inline double sph_to_cart_y(double k, double theta_k, double phi_k){
    return k * sin(theta_k) * sin(phi_k);
}
inline double sph_to_cart_z(double k, double theta_k){
    return k * cos(theta_k);
}


inline double window_tophat(double kx, double ky, double kr, 
                     double x0, double y0, double r0,
                     double sgn, double w, double dr)
{
    /*
    Fourier-space tophat window function.
    */
    // FIXME: Missing factor of (Delta r) . w^2
    // FIXME: This is the real part only!
    return gsl_sf_sinc(kx * w) * gsl_sf_sinc(ky * w) * gsl_sf_sinc(kr * dr) 
           * cos(sgn * (kx*x0 + ky*y0 + kr*r0));
}


double angle_avg_window(double q, double kx, double ky, double kr,
                        double kxp, double kyp, double krp,
                        double w, double dr)
{
    /*
    Calculate the Fourier angle-averaged window function for a given (k, k') 
    pair.
    */
    double dth = M_PI / (double)NTHETA;
    double dphi = 2.*M_PI / (double)NPHI;
    
    int i, j;
    double qth, qphi;
    double qx, qy, qr;
    double wfn;
    
    // FIXME
    double x0 = 0.;
    double y0 = 0.;
    double r0 = 0.;
    
    // Loop over theta and phi, doing simple box integration
    wfn = 0.;
    for (i=0; i < NTHETA; i++){
        // Get r coord (depends only on theta)
        qth = (double)i * dth;
        qr = sph_to_cart_z(q, qth);
        
        // Loop over phi
        for (j=0; j < NPHI; j++){
            qphi = (double)j * dphi;
            
            // Get x and y coords
            qx = sph_to_cart_x(q, qth, qphi);
            qy = sph_to_cart_y(q, qth, qphi);
            
            // Add to integral
            wfn += window_tophat(kx - qx, ky - qy, kr - qr,
                                 x0, y0, r0, +1., w, dr)
                 * window_tophat(kxp - qx, kyp - qy, krp - qr, 
                                 x0, y0, r0, -1., w, dr)
                 * sin(qth);
        } // end phi loop
    } // end theta loop
    
    // Multiply by element size and return
    return dth * dphi * wfn;
}


static char docstring_mode_cpl_fn[] = 
  "Create a realisation of galaxy physical properties on top of an input halo\n"
  "catalogue.\n\n"
  "Parameters\n"
  "----------\n"
  "mhalo, z : array_like (must have dtype==np.float32)\n"
  "  Arrays containing the halo mass [M_sun] and redshift of the objects in \n"
  "  the halo catalogue.\n\n"
  "params : dict (optional)\n"
  "  Dictionary of ghost model parameters. Uses default values for parameters\n"
  "  that are not specified. Uses all default values if argument not passed.\n\n"
  "Returns\n"
  "-------\n"
  "mstar, sfr, passive : array_like (float, float, bool)\n"
  "  Realisations of the stellar mass [M_sun], star-formation rate [M_sun/yr]\n"
  "  and whether the galaxy is passive (True) or star-forming (False).";

static PyObject* mode_cpl_fn(PyObject *self, PyObject *args){
    /*
    Calculate the mode-coupling function in (|k|, |k|') by integrating P(k) 
    over the Fourier-space window functions.
    */
    
    // Input Fourier vectors and parameters
    PyObject *arg_q;
    PyObject *arg_k, *arg_theta_k, *arg_phi_k;
    PyObject *arg_kp, *arg_theta_kp, *arg_phi_kp;
    double arg_w, arg_dr;
    //PyObject *arg_params;
    
    // arg_params is optional, so initialise as empty dictionary first
    //arg_params = PyDict_New();
    
    // Get halo catalogue arrays from input arguments
    if (!PyArg_ParseTuple(args, "OOOOOOOdd", &arg_q,
                          &arg_k, &arg_theta_k, &arg_phi_k,
                          &arg_kp, &arg_theta_kp, &arg_phi_kp,
                          &arg_w, &arg_dr)){
        PyErr_SetString(PyExc_RuntimeError, 
                        "Failed to interpret input arguments.");
        return NULL;
    }
    
    // Convert input arguments to numpy arrays and expose C pointers to data
    PyArrayObject *np_q        = (PyArrayObject*)PyArray_FROM_O(arg_q);
    PyArrayObject *np_k        = (PyArrayObject*)PyArray_FROM_O(arg_k);
    PyArrayObject *np_theta_k  = (PyArrayObject*)PyArray_FROM_O(arg_theta_k);
    PyArrayObject *np_phi_k    = (PyArrayObject*)PyArray_FROM_O(arg_phi_k);
    PyArrayObject *np_kp       = (PyArrayObject*)PyArray_FROM_O(arg_kp);
    PyArrayObject *np_theta_kp = (PyArrayObject*)PyArray_FROM_O(arg_theta_kp);
    PyArrayObject *np_phi_kp   = (PyArrayObject*)PyArray_FROM_O(arg_phi_kp);
    
    // Get length of input arrays
    int N = (int)PyArray_DIM(np_k, 0);
    int NQ = (int)PyArray_DIM(np_q, 0);
    
    // Create new ndarrays and provide data access pointers for C code
    int ndim = 2;
    npy_intp shape[2] = {N, NQ};
    PyObject *np_wfn = PyArray_SimpleNew(ndim, shape, NPY_DOUBLE);
    if (np_wfn == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "Building output arrays failed.");
        Py_XDECREF(np_wfn);
        return NULL;
    }
    
    // Get references to data structures for galaxy properties
    double *q        = (double*)PyArray_DATA(np_q);
    double *k        = (double*)PyArray_DATA(np_k);
    double *theta_k  = (double*)PyArray_DATA(np_theta_k);
    double *phi_k    = (double*)PyArray_DATA(np_phi_k);
    double *kp       = (double*)PyArray_DATA(np_kp);
    double *theta_kp = (double*)PyArray_DATA(np_theta_kp);
    double *phi_kp   = (double*)PyArray_DATA(np_phi_kp);
    double *wfn      = (double*)PyArray_DATA(np_wfn);
    
    double kx, ky, kr, kxp, kyp, krp;
    
    // Loop over (k, k') pairs
    for(int i=0; i < N; i++){
        printf("\t(%d / %d) k = %3.3e, k' = %3.3e\n", i+1, N, k[i], kp[i]);
        
        // Get Cartesian representations of (k, k') vectors
        kx = sph_to_cart_x(k[i], theta_k[i], phi_k[i]);
        ky = sph_to_cart_y(k[i], theta_k[i], phi_k[i]);
        kr = sph_to_cart_z(k[i], theta_k[i]);
        
        kxp = sph_to_cart_x(kp[i], theta_kp[i], phi_kp[i]);
        kyp = sph_to_cart_y(kp[i], theta_kp[i], phi_kp[i]);
        krp = sph_to_cart_z(kp[i], theta_kp[i]);
        
        // Loop over values of |q|
        #pragma omp parallel for
        for(int j=0; j < NQ; j++){
            *(wfn + NQ*i + j) = angle_avg_window(q[j], kx, ky, kr, 
                                                 kxp, kyp, krp, arg_w, arg_dr);
        } // end j loop
    } // end i loop
    
    // Clean-up references
    //Py_DECREF(arg_params);
    Py_DECREF(np_q);
    Py_DECREF(np_k); Py_DECREF(np_theta_k); Py_DECREF(np_phi_k);
    Py_DECREF(np_kp); Py_DECREF(np_theta_kp); Py_DECREF(np_phi_kp);
    //Py_DECREF(dtype_mhalo); Py_DECREF(dtype_z);
    
    // Construct tuple of arrays to be returned
    PyObject *winfn = Py_BuildValue("O", np_wfn);
    return winfn;
}


////////////////////////////////////////////////////////////////////////////////
// Define public methods and initialisation routine
////////////////////////////////////////////////////////////////////////////////

static struct PyMethodDef methods[] = {
    {"mode_cpl_fn", mode_cpl_fn, METH_VARARGS, 
        docstring_mode_cpl_fn},
    {NULL, NULL, 0, NULL} // Sentinel block
};

PyMODINIT_FUNC initsparse(void){
    (void)Py_InitModule("sparse", methods);
    import_array();
}
