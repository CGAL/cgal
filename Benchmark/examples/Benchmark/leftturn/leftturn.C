#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>

#include <CGAL/Bench.h>
#include <CGAL/Bench.C>
#include <CGAL/Bench_parse_args.h>
#include <CGAL/Bench_parse_args.C>

typedef CGAL::Quotient<CGAL::MP_Float>          NT;
typedef CGAL::Cartesian<NT>                     Kernel;

template <class Kernel>
class Bench_leftturn : public Kernel {
private:
  typename Kernel::Point_2 p, q, r;
  typename Kernel::Leftturn_2 leftturn;
public:
  Bench_leftturn()
  {
    leftturn = leftturn_2_object();
    p = typename Kernel::Point_2(0,0);
    q = typename Kernel::Point_2(1,1);
    r = typename Kernel::Point_2(0,1);
  }
  int init(void) { return 0; }
  void clean(void) {}
  void sync(void) {}
  void op(void) { leftturn(p, q, r); }
};

int main(int argc, char * argv[])
{
  CGAL::Bench_parse_args parseArgs(argc, argv);
  int rc = parseArgs.parse();
  if (rc > 0) return 0;
  if (rc < 0) return rc;
  int seconds = parseArgs.getSeconds();
  CGAL::Bench_base::printHeader();
  CGAL::Bench<Bench_leftturn<Kernel> > benchInst("Leftturn", seconds, false);
  benchInst();
  return 0;
}
