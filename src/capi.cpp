#include "capi.h"
#include "solver.h"


DiffEqParameters* create_DiffEqParameters(unsigned int dimensions) 
{
  return new DiffEqParameters(dimensions);
}

void free_DiffEqParameters(DiffEqParameters* params) 
{ 
  delete params; 
}


EBITChargeBreedingSimulation* 
create_EBITChargeBreedingSimulation(unsigned int no_dimensions, DiffEqParameters* params) 
{
  return new EBITChargeBreedingSimulation(no_dimensions, params);
}

void free_EBITChargeBreedingSimulation(EBITChargeBreedingSimulation* sim) 
{ 
  delete sim; 
}


ValuesPerNuclide* create_ValuesPerNuclide(unsigned int no_values) 
{
  return new ValuesPerNuclide(no_values);
}

void free_ValuesPerNuclide(ValuesPerNuclide* values) 
{ 
  delete values; 
}

void free_Result(Result* result) 
{
  delete result;
}


Result* solve_ode(EBITChargeBreedingSimulation* simulation)
{
  return do_solve(simulation);
}
