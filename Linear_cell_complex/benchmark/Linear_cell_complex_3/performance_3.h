//=============================================================================

#ifndef PERFORMANCE_TEST_3_H
#define PERFORMANCE_TEST_3_H


//== INCLUDES =================================================================

#include <cstdlib>
#include <iostream>
#include <assert.h>
#include <sys/resource.h>

#include "Stop_watch.h"


//== CLASS DEFINITION =========================================================


class Performance_test_3
{
public:

  Performance_test_3() {}

  void run(const char* input, const char* output)
  {
    graphene::StopWatch timer;

    timer.start();
    if (!read_mesh(input)) { std::cerr << "read error\n"; exit(1); }
    timer.stop();
    std::cout << "Read mesh   : " << timer << std::endl;
    display_info();

    timer.start();
    int c;
    for (int i=0; i<7; ++i)
      c = circulator_test();
    timer.stop();
    assert(c==0);
    std::cout << "Circulator: "<<c<<"; time: " << timer << std::endl;

    for (int i=0; i<7; ++i)
      c = circulator2_test();
    timer.stop();
    std::cout << "Circulator2: "<<c<<"; time: " << timer << std::endl;

    timer.start();
    for (int i=0; i<7; ++i)
      barycenter_test(i==0);
    timer.stop();
    std::cout << "Barycenter  : " << timer << std::endl;

    timer.start();
    for (int i=0; i<7; ++i)
      smoothing_test();
    timer.stop();
    std::cout << "Smoothing   : " << timer << std::endl;

    timer.start();
    split_tet_test();
    timer.stop();
    std::cout << "Split tetrahedra : " << timer << std::endl;
    display_info();

    timer.start();
    collapse_test(7);
    timer.stop();
    std::cout << "Collapse    : " << timer << std::endl;
    display_info();

    std::cout << std::endl;
  }


protected:

  virtual bool read_mesh(const char*) = 0;
  virtual bool write_mesh(const char*) = 0;
  virtual int  circulator_test()  { return 0; }
  virtual int  circulator2_test()  { return 0; }
  virtual void barycenter_test(bool)  {}
  virtual void smoothing_test()   {}
  virtual void split_tet_test()   {}
  virtual void collapse_test()    {}
  virtual void collapse_test(unsigned int n)    {}
  virtual void display_info()     {}
};


//=============================================================================
#endif
//=============================================================================
