#include <math.h>

#define CGAL_BEGIN_NAMESPACE namespace CGAL {
#define CGAL_END_NAMESPACE }

#include <CGAL/Bench.h>
#include <CGAL/Bench_option_parser.h>

class Bench_sqrt {
private:
  double n;

public:
  Bench_sqrt() : n(M_PI) {}
  int init(void) { return 0; }
  void clean(void) {}
  void sync(void) {}
  void op(void) { sqrt(n); }
};

namespace po = boost::program_options;

struct Help_exception {};

int main(int argc, char * argv[])
{
  po::options_description opts("Options");
  opts.add_options()("help,h", "print help message");
  CGAL::Bench_option_parser bench_opts;
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

  CGAL::Bench<Bench_sqrt> bench("Square root", bench_opts.get_seconds());
  bench();
  return 0;
}
