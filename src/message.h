#ifndef MESSAGE_H
#define MESSAGE_H
#include <vector>

extern "C" {

  struct DiffEqParameters;

  double inactive_ion_ion_overlap(DiffEqParameters* userdata, double beta1, double beta2);
  double inactive_electron_beam_ion_overlap(DiffEqParameters* userdata, double kt);
  double inactive_one_over_ioncloud_vol(DiffEqParameters* userdata, double kt);
  double inactive_heat_capacity(DiffEqParameters* userdata, double kt);


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
    InitialValue* initial_values;
    double initial_temperature;

    double (*ion_ion_overlap)(DiffEqParameters*, double beta1, double beta2) = inactive_ion_ion_overlap;
    double (*electron_beam_ion_overlap)(DiffEqParameters*, double kt) = inactive_electron_beam_ion_overlap;
    double (*one_over_ioncloud_vol)(DiffEqParameters*, double kt) = one_over_ioncloud_vol;
    double (*heat_capacity)(DiffEqParameters*, double kt) = heat_capacity;
    
    DiffEqParameters(unsigned int _no_dimensions) 
    : no_dimensions(_no_dimensions) 
    {
      qV_e = new double[_no_dimensions];
      qV_t = new double[_no_dimensions];
      A = new double[_no_dimensions];
      phi = new double[_no_dimensions];
      q = new double[_no_dimensions];
      source_n = new double[_no_dimensions];
      source_kt = new double[_no_dimensions];

      Xi_ij = new double[_no_dimensions * _no_dimensions];
      eta_ij = new double[_no_dimensions * _no_dimensions];
      CX_ij = new double[_no_dimensions * _no_dimensions];	    
      initial_values = new InitialValue[_no_dimensions];
    }

    ~DiffEqParameters()
    {
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

  struct EBITChargeBreedingSimulation 
  {
    double start_time;
    double end_time;

    unsigned int no_timesteps;
    double* timesteps;
    DiffEqParameters* diff_eq_parameters;

  EBITChargeBreedingSimulation(unsigned int _no_timesteps, DiffEqParameters* _params) 
    : no_timesteps(_no_timesteps), diff_eq_parameters(_params)
    {
      timesteps = new double[_no_timesteps];
    }

    ~EBITChargeBreedingSimulation() { delete[] timesteps; }
  }; 


  struct ValuesPerNuclide
  {
    unsigned int A;
    unsigned int Z;
    unsigned int q;
    unsigned int i;

    unsigned int no_values;
    double* values;

  ValuesPerNuclide(unsigned int _no_values) : no_values(_no_values)
    {
      values = new double[_no_values];
    }
    ~ValuesPerNuclide() { delete[] values; }
  };

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

  struct Result {
    ReturnCode return_code;
    unsigned int no_timesteps;
    double* times;
    unsigned int no_nuclides;
    ValuesPerNuclide* n;
    ValuesPerNuclide* kT;
    std::vector<ValuesPerNuclide>* kt_storage;
    std::vector<ValuesPerNuclide>* n_storage;

    Result(unsigned int _no_timesteps, unsigned int _no_nuclides)
      : no_timesteps(_no_timesteps), no_nuclides(_no_nuclides)
    {
      times = new double[_no_timesteps];
      n_storage = new std::vector<ValuesPerNuclide>(_no_nuclides, ValuesPerNuclide(_no_timesteps));
      kt_storage = new std::vector<ValuesPerNuclide>(_no_nuclides, ValuesPerNuclide(_no_timesteps));
      n = n_storage->data();
      kT = kt_storage->data();
    }
  };
}
#endif /* MESSAGE_H */
