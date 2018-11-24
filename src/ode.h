#ifndef ODE_H
#define ODE_H


#include <ebit-ode-messages.pb.h>


extern "C" {
    void solve_ode(const char* msg, unsigned int size, 
		   char** answer, unsigned int* answer_size);
    void free_answer(char* answer);
}

#endif /* ODE_H */
