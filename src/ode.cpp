
#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <tuple>

#include <boost/numeric/odeint.hpp>

#include "ode.h"
#include <ebit-ode-messages.pb.h>


using namespace std;

typedef std::vector<double> state_type;

typedef ::google::protobuf::RepeatedPtrField< ::EbitODEMessages::MatrixValue > matrix;
typedef ::google::protobuf::RepeatedPtrField< ::EbitODEMessages::InitialValue > initValue;


double* create_matrix(int dimension,  const matrix& sparse)
{
    double* ret_val = new double[dimension*dimension]; 
    for (auto ptr = sparse.begin(); ptr < sparse.end(); ++ptr){
	// to keep it in sync with julia, decf 1 on array offsets
	auto i = ptr->row() - 1;
	auto j = ptr->column() - 1;
	ret_val[dimension*i + j] = ptr->value();
    }
	
    return ret_val;
}

state_type& create_init_values(int dimension,  const initValue& init_values)
{
    auto ret_val = new state_type(2*dimension); // hopefully this initializes to zero
    for (auto ptr = init_values.begin(); ptr < init_values.end(); ++ptr){
	auto i = ptr->index() - 1;
	(*ret_val)[i] = ptr->number_of_particles();
	(*ret_val)[dimension + i] = ptr->temperature_in_ev();
    }
    return *ret_val;
}



class ebit_ode {
    const double* qV_e;
    const double* qV_t;
    const double* A;
    const double* phi;
    const double* qVe_over_Vol_x_kT;
    const double* source;

    const double* Xi_ij;
    const double* dN_ij;
    const double* CX_ij;

    unsigned int no_dimensions;
    double min_N;


public:
    ~ebit_ode() { }
    ebit_ode(const EbitODEMessages::DiffEqParameters& p)
	{
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
    void operator()(const state_type& x, state_type& dxdt, const double)
	{
	    auto tau = [&](int i){ return std::max(0.0, x[no_dimensions+ i]); };
	    
        
	    for (int i = 0; i < no_dimensions; ++i) {
		double R_esc_sum_j = 0.0;
		double R_exchange_sum_j = 0.0;
		dxdt[i] = source[i];
		dxdt[i + no_dimensions] = 0.0;


		for (int j = 0; j < no_dimensions; ++j) {
			
		    if (x[i] > min_N && x[j] > min_N 
			&& tau(j) > 0.0 
			&& tau(i) > 0.0) 
		    {
			double f_ij = std::min((tau(i)*qV_e[j]) / (tau(j)*qV_e[i]), 1.0);
			double n_j = x[j] * qVe_over_Vol_x_kT[j] / tau(j);
			double arg = (tau(i)/A[i] + tau(j)/A[j]);
			double Sigma = Xi_ij[i*no_dimensions + j] * n_j * pow(arg, -1.5);
			R_esc_sum_j += f_ij * Sigma;
			R_exchange_sum_j +=  f_ij * Sigma * (tau(j) - tau(i));

			dxdt[i] += CX_ij[i*no_dimensions + j]*x[j]*sqrt(tau(j));
		    }

		    dxdt[i] += dN_ij[i*no_dimensions + j]*x[j];

		}

		dxdt[i] += R_exchange_sum_j;
		if (x[i] > min_N){
		    double R_esc = 3/sqrt(3) * R_esc_sum_j 
			* ( tau(i) / qV_t[i] ) 
			* exp( -qV_t[i] / tau(i) );
		    dxdt[i] += ( std::min( qV_e[i] / tau(i), 1.0) * phi[i] ) 
			- ( tau(i) + qV_t[i] ) * R_esc;
		    dxdt[i] -= x[i] * R_esc;
		}
	    }
	}
};



struct push_back_state_and_time
{
    std::vector< state_type >& m_states;
    std::vector< double >& m_times;

    push_back_state_and_time( std::vector< state_type > &states , std::vector< double > &times )
	: m_states( states ) , m_times( times ) { }

    void operator()( const state_type &x , double t )
	{
	    m_states.push_back( x );
	    m_times.push_back( t );
	}
};


struct write_state
{
    void operator()( const state_type &x , double time ) const
	{
	    std::cout << "time: " << time << "\t" << x[0] << "\t" << x[1] << "\n";
	}
};

EbitODEMessages::Result* do_solve(const ebit_ode& ode, 
				  const EbitODEMessages::SolverParameters& solver_params, 
				  const EbitODEMessages::DiffEqParameters& diff_params)
{
    using namespace boost::numeric::odeint;
    
    std::cout << "Start solving problem of size: " << diff_params.no_dimensions() << std::endl;
    
    auto saveat = solver_params.saveat();
    auto x = create_init_values(diff_params.no_dimensions(), diff_params.initial_values());

    typedef runge_kutta_cash_karp54< state_type > error_stepper_type;
    typedef controlled_runge_kutta< error_stepper_type > controlled_stepper_type;
    controlled_stepper_type controlled_stepper;


    auto ret_val = new EbitODEMessages::Result();

    auto ptr = saveat.begin();
    double last_time = *ptr++;
    for (; ptr < saveat.end(); ++ptr){
	double to_time = *ptr;
	integrate_adaptive( controlled_stepper , ode , x , 
			    last_time , to_time , (to_time - last_time));
	last_time = to_time;
    }
    
    return ret_val;
}

void solve_ode(const char* msg_buffer, unsigned int size, 
	       char** answer_buffer, unsigned int* answer_size)
{


    EbitODEMessages::Message msg;
    if (!msg.ParseFromArray(msg_buffer, size))
	throw "Couldn't parse message.";
    auto problem = msg.ode_problem();
    auto diff_params = problem.diff_eq_parameters();
    auto ode = ebit_ode(diff_params);
    
    auto result = do_solve(ode, problem.solver_parameters(), diff_params);
    result->set_allocated_problem(&problem);
}


