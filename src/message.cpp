
#include "message.h"

double inactive_ion_ion_overlap(DiffEqParameters* userdata, double beta1, double beta2) { return 1.0; }
double inactive_electron_beam_ion_overlap(DiffEqParameters* userdata, double kt) { return 1.0; }
double inactive_one_over_ioncloud_vol(DiffEqParameters* userdata, double kt) { return 1.0; }
double inactive_heat_capacity(DiffEqParameters* userdata, double kt) { return 1.0; }

