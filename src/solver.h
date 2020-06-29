#ifndef SOLVER_H
#define SOLVER_H


#include "message.h"
#include "ode.h"

extern "C" {
    void solve_ode(const EBITChargeBreedingSimulation *simulation, Result **result);
    void free_answer(Result* answer);
}


#endif /* SOLVER_H */
