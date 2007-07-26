//! \file examples/Benchmark/bbox.cpp
// Handling intersection points with multiplicity between conic arcs.
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

#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>

#include <CGAL/Benchmark/Benchmark.hpp>
#include <CGAL/Benchmark/Option_parser.hpp>

typedef CGAL::Quotient<CGAL::MP_Float>          NT;
typedef CGAL::Cartesian<NT>                     Kernel;

namespace po = boost::program_options;

struct Help_exception {};

template <class Kernel>
class Bench_bbox : public Kernel {
public:
  int init(void) { return 0; }
  void clean(void) {}
  void sync(void) {}
  void op(void) {}
};

typedef Bench_bbox<Kernel>                      My_bench_bbox;

namespace po = boost::program_options;
namespace cb = CGAL::benchmark;

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

  cb::Benchmark<My_bench_bbox> bench("Leftturn", bench_opts.get_seconds());
  bench();
  return 0;
}

#endif
