#ifndef SOLVER_H
#define SOLVER_H


#include "ebit-ode-messages.pb.h"
#include "ode.h"

extern "C" {
    void solve_ode(const char* msg, unsigned int size, 
		   char** answer, unsigned int* answer_size);
    void free_answer(char* answer);
}



/* typedef std::vector<double> state_type; */

/* typedef ::google::protobuf::RepeatedPtrField<::EbitODEMessages::MatrixValue> matrix; */
/* typedef ::google::protobuf::RepeatedPtrField<::EbitODEMessages::InitialValue> initValue; */
/* typedef ::google::protobuf::RepeatedPtrField<::EbitODEMessages::Nuclide> nuclides; */

/* typedef ::google::protobuf::RepeatedField<double> times; */


class write_state {
    const nuclides& m_nuclides;
public:
    write_state(const nuclides& nuclides) : m_nuclides(nuclides) {}
void operator()(const state_type &x, double time) const; 
};


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

    void operator()(const state_type &x, double t); 
};

#endif /* SOLVER_H */
