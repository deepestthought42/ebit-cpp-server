#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <tuple>
#include <vector>
#include <limits>

#include <boost/numeric/odeint.hpp>

#include "ode.h"
#include <ebit-ode-messages.pb.h>

using namespace std;

typedef std::vector<double> state_type;

typedef ::google::protobuf::RepeatedPtrField<::EbitODEMessages::MatrixValue>
matrix;
typedef ::google::protobuf::RepeatedPtrField<::EbitODEMessages::InitialValue>
initValue;
typedef ::google::protobuf::RepeatedPtrField<::EbitODEMessages::Nuclide>
nuclides;

typedef ::google::protobuf::RepeatedField<double> times;

double *create_matrix(int dimension, const matrix &sparse) {
    double *result = new double[dimension * dimension];
    for (auto ptr = sparse.begin(); ptr < sparse.end(); ++ptr) {
	// to keep it in sync with julia, decf 1 on array offsets
	auto i = ptr->row() - 1;
	auto j = ptr->column() - 1;
	result[dimension * i + j] = ptr->value();
    }

    return result;
}

state_type &create_init_values(int dimension, const initValue &init_values) {
    auto result =
	new state_type(2 * dimension); // hopefully this initializes to zero
    for (auto ptr = init_values.begin(); ptr < init_values.end(); ++ptr) {
	auto i = ptr->index() - 1;
	(*result)[i] = ptr->number_of_particles();
	(*result)[dimension + i] = ptr->temperature_in_ev();
    }
    return *result;
}

class ebit_ode {
    const double *qV_e;
    const double *qV_t;
    const double *A;
    const double *phi;
    const double *qVe_over_Vol_x_kT;
    const double *source;

    const double *Xi_ij;
    const double *dN_ij;
    const double *CX_ij;

    unsigned int no_dimensions;
    double min_N;

public:
    ~ebit_ode() {}
    ebit_ode(const EbitODEMessages::DiffEqParameters &p) {
	no_dimensions = p.no_dimensions();
	qV_e = p.qve().data();
	qV_t = p.qvt().data();
	A = p.mass_number().data();
	phi = p.spitzer_divided_by_overlap().data();
	qVe_over_Vol_x_kT = p.qve_over_vol_x_kt().data();
	source = p.source_terms().data();
	min_N = p.minimum_n();
	Xi_ij = create_matrix(no_dimensions, p.inverted_collision_constant());
	dN_ij = create_matrix(no_dimensions, p.rate_of_change_divided_by_n());
	CX_ij = create_matrix(no_dimensions, p.dcharge_ex_divided_by_n_times_tau());
    }
    
    void operator()(const state_type &x, state_type &dxdt, const double) {
	auto tau = [&](int i) { return std::max(0.0, x[no_dimensions + i]); };

	for (int i = 0; i < no_dimensions; ++i) {
	    double R_esc_sum_j = 0.0;
	    double R_exchange_sum_j = 0.0;
	    dxdt[i] = source[i];
	    dxdt[i + no_dimensions] = 0.0;

	    for (int j = 0; j < no_dimensions; ++j) {

		if (x[i] > min_N && x[j] > min_N && tau(j) > 0.0 && tau(i) > 0.0) {
		    double f_ij = std::min((tau(i) * qV_e[j]) / (tau(j) * qV_e[i]), 1.0);
		    double n_j = x[j] * qVe_over_Vol_x_kT[j] / tau(j);
		    double arg = (tau(i) / A[i] + tau(j) / A[j]);
		    double Sigma = Xi_ij[i * no_dimensions + j] * n_j * pow(arg, -1.5);
		    R_esc_sum_j += f_ij * Sigma;
		    R_exchange_sum_j += f_ij * Sigma * (tau(j) - tau(i));

		    dxdt[i] += CX_ij[i * no_dimensions + j] * x[j] * sqrt(tau(j));
		}

		dxdt[i] += dN_ij[i * no_dimensions + j] * x[j];
	    }

	    dxdt[no_dimensions + i] += R_exchange_sum_j;
	    if (x[i] > min_N) {
		double R_esc = 3 / sqrt(3) * R_esc_sum_j * (tau(i) / qV_t[i]) *
		    exp(-qV_t[i] / tau(i));
		dxdt[no_dimensions + i] += (std::min(qV_e[i] / tau(i), 1.0) * phi[i]) -
		    (tau(i) + qV_t[i]) * R_esc;
		dxdt[i] -= x[i] * R_esc;
	    }
	}
    }
};

struct write_state {
    void operator()(const state_type &x, double time) const {
	std::cout << "time: " << time << "\t" << x[0] << "\t" << x[1] << "\n";
    }
};

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


class push_back_state_and_time {
    EbitODEMessages::Result &m_result;
    const nuclides &m_nuclides;
    int m_no_dimensions;
    int m_last_time_i = 0;
    const times& m_times;
public:
    push_back_state_and_time(EbitODEMessages::Result &result, const nuclides &nuclides, 
			     int no_dimensions, const times& times)
	: m_result(result), m_nuclides(nuclides), 
	  m_no_dimensions(no_dimensions), m_times(times)
	{}

    void operator()(const state_type &x, double t) {
	auto saved = save_state(x, m_result, m_nuclides, 
				m_no_dimensions, t,
				m_times.Get(m_last_time_i));
	if (saved && m_last_time_i + 1 < m_times.size())
	    ++m_last_time_i;
    }
};



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

    typedef runge_kutta_cash_karp54<state_type> error_stepper_type;
    typedef controlled_runge_kutta<error_stepper_type> controlled_stepper_type;
    controlled_stepper_type controlled_stepper;

    auto result = prepare_result(nuclides);

    auto ptr = saveat.begin();
    double last_time = *ptr++;

  
    integrate_adaptive(controlled_stepper, ode, x, 
		       problem_params.time_span().start(), 
		       problem_params.time_span().stop(), 1.0,
		       push_back_state_and_time(*result, nuclides, 
						diff_params.no_dimensions(),
						solver_params.saveat()));


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
