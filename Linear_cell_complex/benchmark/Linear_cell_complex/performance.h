//=============================================================================

#ifndef PERFORMANCE_TEST_H
#define PERFORMANCE_TEST_H


//== INCLUDES =================================================================

#include <cstdlib>
#include <iostream>
#include <assert.h>
#include <sys/resource.h>

#include "Stop_watch.h"


//== CLASS DEFINITION =========================================================


class Performance_test
{
public:

    Performance_test() {}

    void run(const char* input, const char* output)
    {
        graphene::StopWatch timer;

        timer.start();
        if (!read_mesh(input)) { std::cerr << "read error\n"; exit(1); }
        timer.stop();
        std::cout << "Read mesh   : " << timer << std::endl;

       timer.start();
       int c;
       for (int i=0; i<100; ++i)
           c = circulator_test();
       timer.stop();
       assert(c==0);
       std::cout << "Circulator  : " << timer << std::endl;

       timer.start();
       for (int i=0; i<1000; ++i)
           barycenter_test();
       timer.stop();
       std::cout << "Barycenter  : " << timer << std::endl;

       timer.start();
       for (int i=0; i<100; ++i)
           normal_test();
       timer.stop();
       std::cout << "Normals     : " << timer << std::endl;

       timer.start();
       for (int i=0; i<100; ++i)
           smoothing_test();
       timer.stop();
       std::cout << "Smoothing   : " << timer << std::endl;

       timer.start();
       subdivision_test();
       timer.stop();
       std::cout << "Subdivision : " << timer << std::endl;

       timer.start();
       collapse_test();
       timer.stop();
       std::cout << "Collapse    : " << timer << std::endl;

       timer.start();
       if (!write_mesh(output)) { std::cerr << "write error\n"; exit(1); }
       timer.stop();
       std::cout << "Write mesh  : " << timer << std::endl;

       std::cout << std::endl;
    }


protected:

    virtual bool read_mesh(const char* _filename) = 0;
    virtual bool write_mesh(const char* _filename) = 0;
    virtual int  circulator_test()  { return 0; }
    virtual void barycenter_test()  {}
    virtual void normal_test()      {}
    virtual void smoothing_test()   {}
    virtual void subdivision_test() {}
    virtual void collapse_test()    {}
};


//=============================================================================
#endif
//=============================================================================
