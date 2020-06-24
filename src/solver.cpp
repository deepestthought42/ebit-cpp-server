#include "solver.h"

#include <cvode/cvode.h>               /* prototypes for CVODE fcts., consts.  */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype      */


#include <iostream>

using namespace std;





void publish(const nuclides& m_nuclides, const state_type &x, double time){
    std::cout << "time: " << time << " ";
    for (auto n = m_nuclides.begin(); n < m_nuclides.end(); ++n) {
	std::cout << ", N(A=" << n->a() << ",Z=" << n->z() << ",q=" << n->q() << "+): "
		  << x[n->i() - 1];
    }
    std::cout << endl;
}

void write_state::operator()(const state_type &x, double time) const {
	publish(m_nuclides, x, time);
}





bool save_state(const state_type &x, EbitODEMessages::Result &result,
                const nuclides &nuclides, int no_dimensions, double last_time, 
		double only_if_larger = 0.0) {
    if (last_time < only_if_larger)
	return false;
    
    result.add_times(last_time);
    for (auto n = nuclides.begin(); n < nuclides.end(); ++n) {
	auto i = n->i() - 1;
	result.mutable_n(i)->add_values(x[i]);
	result.mutable_kt(i)->add_values(x[no_dimensions + i]);
    }

    return true;
}

void push_back_state_and_time::operator()(const state_type &x, double t) {
  // publish(m_nuclides, x, t);
  auto saved = save_state(x, m_result, m_nuclides, m_no_dimensions, t,
                          m_times.Get(m_last_time_i));
  if (saved && m_last_time_i + 1 < m_times.size())
    ++m_last_time_i;
}

EbitODEMessages::Result *prepare_result(const nuclides &nuclides) {
    auto result = new EbitODEMessages::Result();

    for (auto p = nuclides.begin(); p < nuclides.end(); ++p) {
	auto n = result->add_n();
	auto kt = result->add_kt();

	n->set_allocated_nuclide(new EbitODEMessages::Nuclide(*p));
	kt->set_allocated_nuclide(new EbitODEMessages::Nuclide(*p));
    }
    return result;
}

void initial_values(int dimension, const initValue &init_values, N_Vector y, N_Vector abstol)
{
  
  for (auto ptr = init_values.begin(); ptr < init_values.end(); ++ptr) {
    auto i = ptr->index() - 1;
    NV_Ith_S(abstol, i) = 1e-5;
    NV_Ith_S(abstol, dimension + i) = 1e-3;
    NV_Ith_S(y,i) = ptr->number_of_particles();
    NV_Ith_S(y, dimension + i) = 1.5 * ptr->temperature_in_ev() * ptr->number_of_particles();
  }
}

static int check_retval(void *returnvalue, const char *funcname, int opt)
{
  int *retval;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && returnvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  /* Check if retval < 0 */
  else if (opt == 1) {
    retval = (int *) returnvalue;
    if (*retval < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n",
	      funcname, *retval);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && returnvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  return(0);
}


EbitODEMessages::Result* do_solve(ebit_ode &ode,
				  const EbitODEMessages::SolverParameters &solver_params,
				  const EbitODEMessages::DiffEqParameters &diff_params,
				  const EbitODEMessages::ProblemParameters& problem_params,
				  const nuclides &nuclides) {
    std::cout << "Start solving problem of size: " << diff_params.no_dimensions()
	      << std::endl;

    std::cout << "t1: " << problem_params.time_span().stop()
	      << std::endl;
    
    auto y = N_VNew_Serial(diff_params.no_dimensions() * 2);
    auto abstol = N_VNew_Serial(diff_params.no_dimensions() * 2);
    realtype reltol = 1e-4;
    
    
    auto dimension = diff_params.no_dimensions();
    
    initial_values(dimension, diff_params.initial_values(), y, abstol);

    auto cvode_mem = CVodeCreate(CV_BDF);
    CVodeInit(cvode_mem, f, problem_params.time_span().start(), y);
    CVodeSVtolerances(cvode_mem, reltol, abstol);

    auto A = SUNDenseMatrix(dimension*2, dimension*2);
    if(check_retval((void *)A, "SUNDenseMatrix", 0)) throw "err";
    auto LS = SUNLinSol_Dense(y, A);
    if(check_retval((void *)LS, "SUNLinSol_Dense", 0)) throw "SUNLinSol_Dense";
    
    CVodeSetLinearSolver(cvode_mem, LS, A);
    auto iout = 0;
    auto tout = problem_params.time_span().stop();
    
    realtype t;
    CVodeSetUserData(cvode_mem, &ode);
    auto retval = CVode(cvode_mem, 0.1, y, &t, CV_NORMAL);

    auto result = prepare_result(nuclides);
    result->set_return_code(EbitODEMessages::Success);
    return result;
}

void solve_ode(const char *msg_buffer, unsigned int size, char **answer_buffer,
               unsigned int *answer_size) {

    EbitODEMessages::Message msg;
    if (!msg.ParseFromArray(msg_buffer, size))
	throw "Couldn't parse message.";

    auto problem = msg.ode_problem();
    auto diff_params = problem.diff_eq_parameters();
    auto ode = ebit_ode(diff_params);

    auto answer = new EbitODEMessages::Message();

    auto result = do_solve(ode, problem.solver_parameters(), diff_params, 
			   problem.problem_parameters(), problem.nuclides());
    result->set_allocated_problem(&problem);
    result->set_start_time(0.0);
    result->set_stop_time(0.0);

    answer->set_msg_type(EbitODEMessages::ODEResult);
    answer->set_allocated_ode_result(result);

    *answer_size = answer->ByteSize();
    *answer_buffer = static_cast<char *>(malloc(*answer_size));

    answer->SerializeToArray(*answer_buffer, *answer_size);
}
