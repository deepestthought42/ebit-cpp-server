#if defined(SUNDIALS_EXTENDED_PRECISION)
#define GSYM "Lg"
#define ESYM "Le"
#define FSYM "Lf"
#else
#define GSYM "g"
#define ESYM "e"
#define FSYM "f"
#endif


#include "ode.h"
/*
double *create_matrix(int dimension, const matrix &sparse) {
  double *result = new double[dimension * dimension];

  for (int i = 0; i < 2 * dimension; ++i)
    result[i] = 0.0;
  for (auto ptr = sparse.begin(); ptr < sparse.end(); ++ptr) {
    auto i = ptr->row() - 1;
    auto j = ptr->column() - 1;
    result[dimension * i + j] += ptr->value();
    
  }
  return result;
}
*/

int f(realtype t, N_Vector y, N_Vector ydot, void *user_data) {
  auto e = reinterpret_cast<ebit_ode*>(user_data);
  auto Ntau = [&](int i) { return NV_Ith_S(y, e->no_dimensions + i); };
  auto N = [&](int i) { return NV_Ith_S(y, i); };

  auto dN = [&](int i, double incf) { NV_Ith_S(ydot,i) += incf; };
  auto dNtau = [&](int i, double incf) { NV_Ith_S(ydot, i + e->no_dimensions) += incf; };

  auto set_dN = [&](int i, double set_point) { NV_Ith_S(ydot, i) = set_point; };
  auto set_dNtau = [&](int i, double set_point) { NV_Ith_S(ydot, i + e->no_dimensions) = set_point; };

  auto CX = [&](int i, int j) { return e->CX_ij[e->no_dimensions * i + j]; };
  auto eta = [&](int i, int j) { return e->eta_ij[e->no_dimensions * i + j]; };
  auto Xi = [&](int i, int j) { return e->Xi_ij[e->no_dimensions * i + j]; };

  for (int i = 0; i < e->no_dimensions; ++i) {
    if (N(i) > e->min_N) {
      e->tau[i] = Ntau(i) / (1.5 * N(i));
      e->beta[i] = e->q[i] / e->tau[i];
    } else {
      e->tau[i] = 0.0;
      e->beta[i] = 1.0;
    }
  }

  for (int i = 0; i < e->no_dimensions; ++i) {
    double nu = 0.0;
    double R_esc_sum_j = 0.0;
    double R_exchange_sum_j = 0.0;

    set_dN(i, e->source_n[i]);
    set_dNtau(i, e->source_kt[i]);

    for (int j = 0; j < e->no_dimensions; ++j) {

      if (N(i) > e->min_N && N(j) > e->min_N && e->tau[j] > 0.0 && e->tau[i] > 0.0) {

        double f_ij = e->ion_ion_overlap(e->beta[i], e->beta[j]);
        double n_j = N(j) * e->one_over_ioncloud_vol(e->beta[j]);
        double arg = (e->tau[i] / e->A[i] + e->tau[j] / e->A[j]);
        double Sigma = Xi(i, j) * n_j * pow(arg, -1.5);
        nu += f_ij * Sigma;
        R_exchange_sum_j += f_ij * Sigma * (e->tau[j] - e->tau[i]);

        dN(i, CX(i, j) * N(j) * sqrt(e->tau[j]));
        dNtau(i, CX(i, j) * Ntau(j) * sqrt(e->tau[j]));
      }

      dN(i, eta(i, j) * N(j));
      dNtau(i, eta(i, j) * Ntau(j));
    }
 
    dNtau(i, R_exchange_sum_j * N(i));
    if (N(i) > e->min_N) {
      double omega = e->qV_t[i] / e->tau[i];
      double R_esc = (3 / sqrt(2)) * nu * exp(-omega) / omega;
      dN(i, -N(i) * R_esc);
      dNtau(i, -N(i) * (e->qV_t[i] + e->tau[i]) * R_esc);
    }
  }

  return 0;
}



