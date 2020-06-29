#ifndef MESSAGE_H
#define MESSAGE_H

typedef struct {
  unsigned int A;
  unsigned int Z;
  unsigned int q;
  unsigned int i;
} Nuclide;

typedef struct {
  double value;
  unsigned int row;
  unsigned int column;
} MatrixValue;


typedef struct {
  unsigned int index;
  double number_of_particles;
  double temperature_in_ev;
} InitialValue;


typedef struct {
  unsigned int no_dimensions;
  
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
  

  double (*ion_ion_overlap)(double beta1, double beta2) = [](double beta1, double beta2) { return 1.0; };
  double (*electron_beam_ion_overlap)(double kt) = [](double kt) { return 1.0; };
  double (*one_over_ioncloud_vol)(double kt) = [](double kt) { return 1.0; };
  double (*heat_capacity)(double kt) = [](double kt) { return 1.5; };

  double V_0;
  double r_e;
  double r_dt;
  double l_dt;

  double min_N;
  InitialValue* initial_values;
  double initial_temperature;
  
} DiffEqParameters;

typedef struct {
  double start_time;
  double end_time;

  unsigned int no_timesteps;
  double* timesteps;

  unsigned int no_nuclides;
  Nuclide* nuclides;
  DiffEqParameters* diff_eq_parameters;
  
} EBITChargeBreedingSimulation;


typedef struct {
  unsigned int A;
  unsigned int Z;
  unsigned int q;
  unsigned int i;

  unsigned int no_values;
  double* values;
} ValuesPerNuclide;

enum ReturnCode 
{
 Default,
 Success,
 MaxIters,
 DtLessThanMin,
 Unstable,
 InitialFailure,
 ConvergenceFailure,
 Failure
};

typedef struct {
  ReturnCode return_code;
  unsigned int no_values;
  double* times;
  ValuesPerNuclide* n;
  ValuesPerNuclide* kT;
} Result;

#endif /* MESSAGE_H */
