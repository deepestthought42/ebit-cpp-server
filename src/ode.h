#ifndef ODE_H
#define ODE_H

#include <cvode/cvode.h>               /* prototypes for CVODE fcts., consts.  */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype      */


#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <tuple>
#include <vector>
#include <limits>
#include <algorithm>

#include "ode.h"
#include "ebit-ode-messages.pb.h"


typedef std::vector<double> state_type;

typedef ::google::protobuf::RepeatedPtrField<::EbitODEMessages::MatrixValue> matrix;
typedef ::google::protobuf::RepeatedPtrField<::EbitODEMessages::InitialValue> initValue;
typedef ::google::protobuf::RepeatedPtrField<::EbitODEMessages::Nuclide> nuclides;

typedef ::google::protobuf::RepeatedField<double> times;


struct ebit_ode 
{
  const double *qV_e;
  const double *qV_t;
  const double *A;
  const double *phi;
  const double *q;

  const double *source_n;
  const double *source_kt;
  
  const double *Xi_ij;
  const double *eta_ij;
  const double *CX_ij;

  double *tau;
  double *beta;

  unsigned int no_dimensions;
  double min_N;

  double (*ion_ion_overlap)(double beta1, double beta2) = 
    [](double beta1, double beta2) { return 1.0; };
  double (*electron_beam_ion_overlap)(double kt) = [](double kt) { return 1.0; };
  double (*one_over_ioncloud_vol)(double kt) = [](double kt) { return 1.0; };
  double (*heat_capacity)(double kt) = [](double kt) { return 1.5; };


  public : ~ebit_ode() {}
  ebit_ode(const EbitODEMessages::DiffEqParameters &p);

  //void operator()(const state_type &x, state_type &dxdt, const double);
};

int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);

#endif //ODE_H

