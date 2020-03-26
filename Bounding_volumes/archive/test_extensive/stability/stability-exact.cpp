// Description: Computes the exact miniball of some sets of balls in
// very "bad" (i.e. degenerate) position and saves the result to files
// data/results.* (where * ranges from 0 to 279).
//
// Usage: ./stability-exact from to
// Note: You must call the program from within the directory it resides in,
// and from <= to must be numbers in {0,...,279}.

// #include <CGAL/MP_Float.h>
// #include <CGAL/Quotient.h>
// #include <CGAL/Lazy_exact_nt.h>
// #include <CGAL/double.h>
// #include <CGAL/Filtered_exact.h>
// #include <CGAL/leda_real.h>

#include <CGAL/Gmpq.h>
#include "stability.h"

// All exact computations will be done with the following number-type:
typedef CGAL::Gmpq ExactFT;
// typedef leda_real ExactFT;
// typedef CGAL::Filtered_exact<double,leda_real> ExactFT;
// typedef CGAL::Lazy_exact_nt<CGAL::Quotient<CGAL::MP_Float> > ExactFT;

// A functor to convert pairs (as returned by the miniball package)
// to doubles:
class PairToDouble {
private:
  double sqrt;

public:
  typedef double Result;

  PairToDouble(const ExactFT& disc) :
    sqrt(std::sqrt(to_double(disc))) {}

  const double operator()(const std::pair<ExactFT,ExactFT> p) {
    return to_double(p.first)+sqrt*to_double(p.second);
  }
};

// A helper routine to write a pair (as returned by the miniball package)
// to a stream:
void writePair(const std::pair<ExactFT,ExactFT>& p,
               std::ostream& out,const ExactFT& disc) {
  out << ' ' << p.first;
  if (p.second >= 0)
    out << '+';
  out << p.second << "*sqrt(" << disc << ')';
}

template<int d,int D>
void atomicTest(const int no,const std::string& name,Perturbation perturb,
                const std::vector< Ball<double,D> >& S,
                const BuildFlag) {
  // open the file containing the exact radius and center:
  std::string filename("data/results.");
  filename += toStr(no);
  std::ofstream fd(filename.c_str());
  std::ofstream ed((filename+".exact").c_str());

  // compute the miniball using exact arbitrary-precision arithmetic:
  typedef BallTraits<D,double,ExactFT,CGAL::Tag_false> Traits;
  typedef CGAL::Min_sphere_of_spheres_d<Traits> Minsphere;
  Minsphere mb(S.begin(),S.end());

  // verify result:
  if (!mb.is_valid()) {
    std::cout << std::endl
              << "ERROR!  Computed ball is *not* valid!" << std::endl
              << "Please contact the author kf@iaeth.ch!" << std::endl;
    exit(-1);
  }

  // save exact coordinates and exact radius:
  using std::setprecision;
  fd << D << std::endl;
  ed << D << std::endl;
  PairToDouble p2d(mb.discriminant());
  fd << setprecision(18) << p2d(mb.radius());
  writePair(mb.radius(),ed,mb.discriminant());
  typename Minsphere::Cartesian_const_iterator it=mb.center_cartesian_begin();
  for (int i=0; i<D; ++i) {
    fd << ' ' << setprecision(18) << p2d(*it);
    writePair(*it,ed,mb.discriminant());
  }
  std::cout << "    [" << setprecision(18) << p2d(mb.radius()) << ']';

  // write control information:
  fd << std::endl
     << "Name:         " << name << " (test #" << no << ')' << std::endl
     << "Perturbation: " << static_cast<int>(perturb) << std::endl;
  ed << std::endl
     << "Name:         " << name << " (test #" << no << ')' << std::endl
     << "Perturbation: " << static_cast<int>(perturb) << std::endl;
}

int main(int argnr,char **argv) {
  using std::cout;
  using std::endl;

  // determine which tests to run:
  if (argnr == 3) {
    from = strTo<int>(argv[1]);
    to   = strTo<int>(argv[2]);
    assert(from <= to);
  }

  cout << "---------------------------------------------------------" << endl
       << "Building the miniball testsuite" << endl
       << "---------------------------------------------------------" << endl
       << endl;

  mainTest<BuildFlag>();

  cout << "---------------------------------------------------------" << endl
       << "Done" << endl
       << "---------------------------------------------------------" << endl
       << endl;
}
