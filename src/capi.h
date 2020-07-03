#ifndef _CAPI_H_
#define _CAPI_H_

#include "message.h"

extern "C" {
  DiffEqParameters* create_DiffEqParameters(unsigned int);
  void free_DiffEqParameters(DiffEqParameters*);

  EBITChargeBreedingSimulation* create_EBITChargeBreedingSimulation(unsigned int, DiffEqParameters*);
  void free_EBITChargeBreedingSimulation(EBITChargeBreedingSimulation*);

  ValuesPerNuclide* create_ValuesPerNuclide(unsigned int);
  void free_ValuesPerNuclide(ValuesPerNuclide*);
  
  Result* solve_ode(const EBITChargeBreedingSimulation *simulation);
  void free_result(Result* answer);
}

#endif /* _CAPI_H_ */
















