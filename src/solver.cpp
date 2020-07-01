#include "solver.h"
#include "message.h"

#include <cvode/cvode.h>               /* prototypes for CVODE fcts., consts.  */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype      */


#include <iostream>

using namespace std;


void initial_values(int dimension, const InitialValue* init_values, N_Vector y, N_Vector abstol)
{
  for (int i = 0; i < dimension; ++i) {
    auto no_of_particles = init_values[i].number_of_particles;
    auto temp = init_values[i].temperature_in_ev;
    NV_Ith_S(abstol, i) = 1e-5;
    NV_Ith_S(abstol, dimension + i) = 1e-3;
    NV_Ith_S(y,i) = no_of_particles;
    NV_Ith_S(y, dimension + i) = 1.5 * temp * no_of_particles;
  }
}

static int check_retval(void *returnvalue, const char *funcname, int opt)
{
  int *retval;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && returnvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1);
  }

  /* Check if retval < 0 */
  else if (opt == 1) {
    retval = (int *) returnvalue;
    if (*retval < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n",
	      funcname, *retval);
      return(1); 
    }
  }

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && returnvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  return(0);
}


Result* do_solve(EBITChargeBreedingSimulation* simulation) 
{
  auto diff_eq_params = simulation->diff_eq_parameters;
  auto dimension = diff_eq_params->no_dimensions;

  auto y = N_VNew_Serial(dimension * 2);
  auto abstol = N_VNew_Serial(dimension * 2);
  realtype reltol = 1e-4;
  
  initial_values(dimension, diff_eq_params->initial_values, y, abstol);

  auto cvode_mem = CVodeCreate(CV_BDF);
  CVodeInit(cvode_mem, f, simulation->start_time, y);
  CVodeSVtolerances(cvode_mem, reltol, abstol);

  auto A = SUNDenseMatrix(dimension*2, dimension*2);
  if(check_retval((void *)A, "SUNDenseMatrix", 0)) throw "err";
  auto LS = SUNLinSol_Dense(y, A);
  if(check_retval((void *)LS, "SUNLinSol_Dense", 0)) throw "SUNLinSol_Dense";
    
  CVodeSetLinearSolver(cvode_mem, LS, A);
  auto iout = 0;
  auto tout = simulation->end_time;
    
  realtype t;
  CVodeSetUserData(cvode_mem, diff_eq_params);
  auto retval = CVode(cvode_mem, 0.1, y, &t, CV_NORMAL);
  return NULL;
}



