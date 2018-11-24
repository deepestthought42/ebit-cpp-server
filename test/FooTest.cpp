#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wglobal-constructors"
#pragma clang diagnostic ignored "-Wweak-vtables"

#include <iostream>
#include <stdio.h>
#include <fstream>
#include "gtest/gtest.h"
#include "ode.h"

/**
Template unit test case fixture class. Make a copy of this.

@author Alex Tsui
@date 2011-12-22
*/
class ODETest : public ::testing::Test
{
protected:
    // You can remove any or all of the following functions if its body
    // is empty.

    ODETest() {
        // You can do set-up work for each test here.
    }

    virtual ~ODETest() {
        // You can do clean-up work that doesn't throw exceptions here.
    }

    // If the constructor and destructor are not enough for setting up
    // and cleaning up each test, you can define the following methods:

    virtual void SetUp() {
        // Code here will be called immediately after the constructor (right
        // before each test).
    }

    virtual void TearDown() {
        // Code here will be called immediately after each test (right
        // before the destructor).
    }

    // Objects declared here can be used by all tests in the test case for Foo.
};

// Tests that Foo does Xyz.
TEST_F(ODETest, ReadProto) 
{
    std::ifstream in("/home/renee/tmp/leigh_talk.proto");
    std::string contents((std::istreambuf_iterator<char>(in)), 
			 std::istreambuf_iterator<char>());

    solve_ode(contents.c_str(), contents.length(), nullptr, nullptr);
}





#pragma clang diagnostic pop







