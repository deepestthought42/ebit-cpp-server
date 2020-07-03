#ifndef MESSAGE_H
#define MESSAGE_H
#include <vector>

template <class T>
T *create_array(unsigned int dimension, T default_value)
{
  T* array = new T[dimension];
  for (auto i = 0; i < dimension; ++i)
    array[i] = default_value;
  return array;
}


extern "C" {

struct DiffEqParameters;


double inactive_ion_ion_overlap(DiffEqParameters *userdata, double beta1,
                                double beta2);
double inactive_electron_beam_ion_overlap(DiffEqParameters *userdata,
                                          double kt);
double inactive_one_over_ioncloud_vol(DiffEqParameters *userdata, double kt);
double inactive_heat_capacity(DiffEqParameters *userdata, double kt);

typedef struct {
  unsigned int index = 0;
  double number_of_particles = 0.0;
  double temperature_in_ev = 0.0;
} InitialValue;

struct DiffEqParameters {
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

  double V_0;
  double r_e;
  double r_dt;
  double l_dt;

  double min_N;
  double initial_temperature;
  InitialValue *initial_values;


  double (*ion_ion_overlap)(DiffEqParameters *, double beta1,
                            double beta2) = inactive_ion_ion_overlap;
  double (*electron_beam_ion_overlap)(DiffEqParameters *, double kt) =
      inactive_electron_beam_ion_overlap;
  double (*one_over_ioncloud_vol)(DiffEqParameters *, double kt) = inactive_one_over_ioncloud_vol;
  double (*heat_capacity)(DiffEqParameters *, double kt) = inactive_heat_capacity;

  DiffEqParameters(unsigned int _no_dimensions)
      : no_dimensions(_no_dimensions)
  {
    qV_e = create_array(_no_dimensions, 0.0);
    qV_t = create_array(_no_dimensions, 0.0);
    A = create_array(_no_dimensions, 0.0);
    phi = create_array(_no_dimensions, 0.0);
    q = create_array(_no_dimensions, 0.0);
    source_n = create_array(_no_dimensions, 0.0);
    source_kt = create_array(_no_dimensions, 0.0);


    V_0 = 0.0; 
    r_e = 0.0;
    r_dt = 0.0;
    l_dt = 0.0;

    min_N = 0.0;
    initial_temperature = 0.0;
    
    
    Xi_ij = create_array(_no_dimensions * _no_dimensions, 0.0);
    eta_ij = create_array(_no_dimensions * _no_dimensions, 0.0);
    CX_ij = create_array(_no_dimensions * _no_dimensions, 0.0);

    InitialValue empty = {0, 0.0, 0.0};
    initial_values = create_array(_no_dimensions, empty);
  }

  ~DiffEqParameters() {
    delete[] qV_e;
    delete[] qV_t;
    delete[] A;
    delete[] phi;
    delete[] q;
    delete[] source_n;
    delete[] source_kt;

    delete[] Xi_ij;
    delete[] eta_ij;
    delete[] CX_ij;
    delete[] initial_values;
  }
};

struct EBITChargeBreedingSimulation {
  double start_time;
  double end_time;

  unsigned int no_timesteps;
  double *timesteps;
  DiffEqParameters *diff_eq_parameters;

  EBITChargeBreedingSimulation(unsigned int _no_timesteps,
                               DiffEqParameters *_params)
      : no_timesteps(_no_timesteps), diff_eq_parameters(_params) {
    timesteps = new double[_no_timesteps];
  }

  ~EBITChargeBreedingSimulation() { delete[] timesteps; }
};

struct ValuesPerNuclide {
  unsigned int A;
  unsigned int Z;
  unsigned int q;
  unsigned int i;

  unsigned int no_values;
  double *values;

  ValuesPerNuclide(unsigned int _no_values) : no_values(_no_values) {
    values = new double[_no_values];
  }
  ~ValuesPerNuclide() { delete[] values; }
};

enum ReturnCode {
  Default,
  Success,
  MaxIters,
  DtLessThanMin,
  Unstable,
  InitialFailure,
  ConvergenceFailure,
  Failure
};

struct Result {
  ReturnCode return_code;
  unsigned int no_timesteps;
  double *times;
  unsigned int no_nuclides;

  // using a std::vector is necessary to be able to initialize the
  // values with _no_timesteps
  ValuesPerNuclide *n;
  ValuesPerNuclide *kT;
  std::vector<ValuesPerNuclide> *kt_storage;
  std::vector<ValuesPerNuclide> *n_storage;

  Result(unsigned int _no_timesteps, unsigned int _no_nuclides)
      : no_timesteps(_no_timesteps), no_nuclides(_no_nuclides) {
    times = new double[_no_timesteps];
    n_storage = new std::vector<ValuesPerNuclide>(
        _no_nuclides, ValuesPerNuclide(_no_timesteps));
    kt_storage = new std::vector<ValuesPerNuclide>(
        _no_nuclides, ValuesPerNuclide(_no_timesteps));
    n = n_storage->data();
    kT = kt_storage->data();
  }
  ~Result() {
    delete n_storage;
    delete kt_storage;
  }
};
}
#endif /* MESSAGE_H */
