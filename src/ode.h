#ifndef ODE_H
#define ODE_H


#include <ebit-ode-messages.pb.h>


extern "C" {
    void solve_ode(char* msg, char** answer);
}

#endif /* ODE_H */
