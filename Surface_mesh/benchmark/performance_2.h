//=============================================================================

#ifndef PERFORMANCE_2_TEST_H
#define PERFORMANCE_2_TEST_H


//== INCLUDES =================================================================

#include <cstdlib>
#include <iostream>
//#include <sys/resource.h>
#include <CGAL/assertions.h>
#include <CGAL/Timer.h>
#include "Stop_watch.h"


//== CLASS DEFINITION =========================================================


class Performance_test_2
{
public:

  Performance_test_2() {}

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
    for (int i=0; i<100; ++i)
      c = circulator_test();
    timer.stop();
    CGAL_assertion(c==0);
    std::cout << "Circulator: "<<c<<"; time: " << timer << std::endl;

    timer.start();
    for (int i=0; i<1000; ++i)
      barycenter_test(i==0);
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
    display_info();

    timer.start();
    collapse_test();
    timer.stop();
    std::cout << "Collapse    : " << timer << std::endl;
    display_info();

    timer.start();
    lindstrom_test(input);
    timer.stop();
    std::cout << "Lindstrom   : " << timer << std::endl;
    display_info();

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
  virtual void barycenter_test(bool)  {}
  virtual void normal_test()      {}
  virtual void smoothing_test()   {}
  virtual void subdivision_test() {}
  virtual void lindstrom_test(const char* _filename)=0;
  virtual void collapse_test()    {}
  virtual void display_info()     {}
};


//=============================================================================
#endif
//=============================================================================
