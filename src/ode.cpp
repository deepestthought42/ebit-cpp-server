#include "ode.h"

double *create_matrix(int dimension, const matrix &sparse) {
  double *result = new double[dimension * dimension];

  for (int i = 0; i < 2 * dimension; ++i)
    result[i] = 0.0;
  for (auto ptr = sparse.begin(); ptr < sparse.end(); ++ptr) {
    // to keep it in sync with julia, decf 1 on array offsets
    auto i = ptr->row() - 1;
    auto j = ptr->column() - 1;
    result[dimension * i + j] += ptr->value();
    // std::cout << "R_{" << i << "," << j << "}: " << ptr->value() <<
    // std::endl;
  }
  // std::cout << endl;

  return result;
}


ebit_ode::ebit_ode(const EbitODEMessages::DiffEqParameters &p) {
  no_dimensions = p.no_dimensions();
  tau = new double[p.no_dimensions()];
  beta = new double[p.no_dimensions()];
  qV_e = p.qve().data();
  qV_t = p.qvt().data();
  A = p.mass_number().data();
  phi = p.spitzer_divided_by_overlap().data();
  q = p.q().data();
  source_n = p.source_terms_n().data();
  source_kt = p.source_terms_kt().data();
  min_N = p.minimum_n();
  Xi_ij = create_matrix(no_dimensions, p.inverted_collision_constant());
  eta_ij = create_matrix(no_dimensions, p.rate_of_change_divided_by_n());
  CX_ij = create_matrix(no_dimensions, p.dcharge_ex_divided_by_n_times_tau());
}

void ebit_ode::operator()(const state_type &x, state_type &dxdt, const double) {
  auto Ntau = [&](int i) { return std::max(0.0, x[no_dimensions + i]); };
  auto N = [&](int i) { return std::max(0.0, x[i]); };

  auto dN = [&](int i, double incf) { dxdt[i] += incf; };
  auto dNtau = [&](int i, double incf) { dxdt[i + no_dimensions] += incf; };

  auto set_dN = [&](int i, double set_point) { dxdt[i] = set_point; };
  auto set_dNtau = [&](int i, double set_point) {
    dxdt[i + no_dimensions] = set_point;
  };

  auto CX = [&](int i, int j) { return CX_ij[no_dimensions * i + j]; };
  auto eta = [&](int i, int j) { return eta_ij[no_dimensions * i + j]; };
  auto Xi = [&](int i, int j) { return Xi_ij[no_dimensions * i + j]; };

  for (int i = 0; i < no_dimensions; ++i) {
    if (N(i) > min_N) {
      tau[i] = Ntau(i) / (1.5 * N(i));
      beta[i] = q[i] / tau[i];
    } else {
      tau[i] = 0.0;
      beta[i] = 1.0;
    }
  }

  for (int i = 0; i < no_dimensions; ++i) {
    double nu = 0.0;
    double R_esc_sum_j = 0.0;
    double R_exchange_sum_j = 0.0;

    set_dN(i, source_n[i]);
    set_dNtau(i, source_kt[i]);

    for (int j = 0; j < no_dimensions; ++j) {

      if (N(i) > min_N && N(j) > min_N && tau[j] > 0.0 && tau[i] > 0.0) {

        double f_ij = ion_ion_overlap(beta[i], beta[j]);
        double n_j = N(j) * one_over_ioncloud_vol(beta[j]);
        double arg = (tau[i] / A[i] + tau[j] / A[j]);
        double Sigma = Xi(i, j) * n_j * pow(arg, -1.5);
        nu += f_ij * Sigma;
        R_exchange_sum_j += f_ij * Sigma * (tau[j] - tau[i]);

        dN(i, CX(i, j) * N(j) * sqrt(tau[j]));
        dNtau(i, CX(i, j) * Ntau(j) * sqrt(tau[j]));
      }

      dN(i, eta(i, j) * N(j));
      dNtau(i, eta(i, j) * Ntau(j));
    }

    dNtau(i, R_exchange_sum_j * N(i));
    if (N(i) > min_N) {
      double omega = qV_t[i] / tau[i];
      double R_esc = (3 / sqrt(2)) * nu * exp(-omega) / omega;
      dN(i, -N(i) * R_esc);
      dNtau(i, -N(i) * (qV_t[i] + tau[i]) * R_esc);
    }
  }
}
