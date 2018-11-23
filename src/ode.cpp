
#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>

#include "ode.h"
#include <ebit-ode-messages.pb.h>


typedef double* state_type;


class ebit_ode {
    double* qV_e;
    double* qV_t;
    double* A;
    double* phi;
    double* qVe_over_Vol_x_kT;
    double* source;

    double* Xi_ij;
    double* dN_ij;
    double* CX_ij;

    unsigned int no_dimensions;
    double min_N;


public:
    void operator()(const state_type& x, state_type& dxdt, const double)
	{
	    double* N = x;
	    double* tau = x + no_dimensions;

	    double* dN = dxdt;
	    double* dtau = dxdt + no_dimensions;
	    
	    for (int i = 0; i < no_dimensions; ++i){
		tau[i] = std::max(0.0, tau[i]);
	    }
        
	    for (int i = 0; i < no_dimensions; ++i) {
		double R_esc_sum_j = 0.0;
		double R_exchange_sum_j = 0.0;
		dN[i] = source[i];
		dtau[i] = 0.0;


		for (int j = 0; j < no_dimensions; ++j) {
			
		    if (N[i] > min_N && N[j] > min_N && tau[j] > 0.0 && tau[i] > 0.0) {
			double f_ij = std::min((tau[i]*qV_e[j])/(tau[j]*qV_e[i]), 1.0);
			double n_j = N[j] * qVe_over_Vol_x_kT[j] / tau[j];
			double arg = (tau[i]/A[i] + tau[j]/A[j]);
			double Sigma = Xi_ij[i*no_dimensions + j] * n_j * pow(arg, -1.5);
			R_esc_sum_j += f_ij * Sigma;
			R_exchange_sum_j +=  f_ij * Sigma * (tau[j] - tau[i]);

			dN[i] += CX_ij[i*no_dimensions + j]*N[j]*sqrt(tau[j]);
		    }

		    dN[i] += dN_ij[i*no_dimensions + j]*N[j];

		}

		dtau[i] += R_exchange_sum_j;
		if (N[i] > min_N){
		    double R_esc = 3/sqrt(3) * R_esc_sum_j * ( tau[i] / qV_t[i] ) * exp( -qV_t[i] / tau[i] );
		    dtau[i] += ( std::min( qV_e[i] / tau[i], 1.0) * phi[i] ) - ( tau[i] + qV_t[i] ) * R_esc;
		    dN[i] -= N[i] * R_esc;
		}
            
	    }

	}

};






void solve_ode(char* msg_buffer, char** )
{
    EbitODEMessages::Message msg;
    msg.ParseFromString(std::string(msg_buffer));
}
