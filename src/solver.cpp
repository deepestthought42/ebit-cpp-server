#include "solver.h"


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

state_type &create_init_values(int dimension, const initValue &init_values) {
  auto result = new state_type(2 * dimension);
  for (auto ptr = init_values.begin(); ptr < init_values.end(); ++ptr) {

    auto i = ptr->index() - 1;
    (*result)[i] = ptr->number_of_particles();
    (*result)[dimension + i] = 1.5 * ptr->temperature_in_ev() * ptr->number_of_particles();
  }
  return *result;
}


EbitODEMessages::Result* do_solve(const ebit_ode &ode,
				  const EbitODEMessages::SolverParameters &solver_params,
				  const EbitODEMessages::DiffEqParameters &diff_params,
				  const EbitODEMessages::ProblemParameters& problem_params,
				  const nuclides &nuclides) {
    using namespace boost::numeric::odeint;

    std::cout << "Start solving problem of size: " << diff_params.no_dimensions()
	      << std::endl;

    auto saveat = solver_params.saveat();
    auto x = create_init_values(diff_params.no_dimensions(),
				diff_params.initial_values());

    typedef runge_kutta_fehlberg78<state_type> error_stepper_type;
    typedef controlled_runge_kutta<error_stepper_type> controlled_stepper_type;
    controlled_stepper_type controlled_stepper;

    auto result = prepare_result(nuclides);

    auto ptr = saveat.begin();
    double last_time = *ptr++;

  
    integrate_adaptive(controlled_stepper, ode, x, 
		       problem_params.time_span().start(), 
		       problem_params.time_span().stop(), 1e-1
		       //write_state(nuclides)
		       ,push_back_state_and_time(*result, nuclides, 
		       				diff_params.no_dimensions(),
		       				solver_params.saveat())
	);


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
