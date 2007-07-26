//! \file examples/Benchmark/simple.cpp
// Measure the performance of computing sqrt(3.14159265358979323846) in double

#include <CGAL/basic.h>

#ifndef CGAL_USE_BOOST_PROGRAM_OPTIONS
#include <iostream>
int main()
{
  std::cout << "Sorry, this example needs boost program options ..."
            << std::endl;
  return 0;
}
#else

#include <math.h>

#include <CGAL/Benchmark/Benchmark.hpp>
#include <CGAL/Benchmark/Option_parser.hpp>

#if (defined _MSC_VER)
#define M_PI 3.14159265358979323846
#endif

class Bench_sqrt {
private:
  double n;

public:
  Bench_sqrt() : n(M_PI) {}
  int init(void) { return 0; }
  void clean(void) {}
  void sync(void) {}
  void op(void) { sqrt(n); }

  void set(double u) { n = u; }
};

namespace po = boost::program_options;
namespace cb = CGAL::benchmark;

struct Help_exception {};

int main(int argc, char * argv[])
{
  po::options_description opts("Options");
  opts.add_options()("help,h", "print help message");
  cb::Option_parser bench_opts;
  opts.add(bench_opts.get_opts());
  po::variables_map var_map;

  try {
    po::store(po::command_line_parser(argc, argv).options(opts).run(), var_map);
    po::notify(var_map);
    if (var_map.count("help")) {
      std::cout << opts << std::endl;
      throw Help_exception();
    }
  } catch(Help_exception & /* e */) {
    return 0;
  } catch(std::exception & e) {
    std::cerr << e.what() << std::endl;
    return 1;
  }

  cb::Benchmark<Bench_sqrt> bench("Square root", bench_opts.get_seconds());
  bench();
  return 0;
}

#endif
