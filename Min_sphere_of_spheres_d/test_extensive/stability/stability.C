// Description: Computes the miniball of some sets of balls in very
// "bad" (i.e. degenerate) position and compares it to the exact
// miniball (which has been precomputed).  The default algorithm
// (Default_algorithm) is checked with and without square-roots for
// doubles and floats.
//
// Usage: ./stability
// Note: You must call the program from within the directory it resides in.
//
// Compatibility: works with or without CGAL installed

#include "stability.h"

// We save the minmal, maximal and average relative error in the radii:
double mine = 100.0, maxe = 0.0, avge = 0.0;
int nr = 0;

double rerr(double exact,double approx) {
  const double e = std::abs((exact-approx)/exact);
  mine = std::min(mine,e);
  maxe = std::max(maxe,e);
  avge += e;
  ++nr;
  return e;
}

template<int d,int D>
void atomicTest(const int test,const std::string& name,Perturbation perturb,
		const std::vector< Ball<double,D> >& S,
		const RunFlag) {
  // open the file containing the exact radius and center:
  std::string filename("data/results.");
  filename += toStr(test);
  std::ifstream ed(filename.c_str());
  
  // read number of coordinates and exact radius:
  int no;
  ed >> no;
  double exact;
  ed >> exact;

  // compute miniball using our different methods:
  typedef CGAL::Tag_true UseSqrt;
  typedef CGAL::Tag_false AvoidSqrt;

  typedef BallTraits<D,double,double,UseSqrt> DUT;
  typedef BallTraits<D,double,double,AvoidSqrt> DAT;
  typedef BallTraits<D,double,float,UseSqrt> FUT;
  typedef BallTraits<D,double,float,AvoidSqrt> FAT;
  CGAL::Min_sphere_of_spheres_d<DUT> dumb(S.begin(),S.end());
  CGAL::Min_sphere_of_spheres_d<DAT> damb(S.begin(),S.end());
  CGAL::Min_sphere_of_spheres_d<DUT> fumb(S.begin(),S.end());
  CGAL::Min_sphere_of_spheres_d<DAT> famb(S.begin(),S.end());

  std::cout << "    Double with sqrts: " << std::setw(50) 
	    << std::setprecision(30) << rerr(exact,dumb.radius())
	    << std::endl;
  std::cout << "    Double w/o sqrts:  " << std::setw(50) 
	    << std::setprecision(30) << rerr(exact,damb.radius())
	    << std::endl;
  std::cout << "    Float with sqrts:  " << std::setw(50) 
	    << std::setprecision(30) << rerr(exact,fumb.radius())
	    << std::endl;
  std::cout << "    Float w/o sqrts:   " << std::setw(50) 
	    << std::setprecision(30) << rerr(exact,famb.radius())
	    << std::endl;
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
       << "Min_sphere_of_spheres_d testsuite" << endl
       << "---------------------------------------------------------" << endl
       << "This program computes the miniball of balls in various" << endl
       << "degenerate configurations.  What it outputs to the screen" << endl
       << "are the relative errors of the computed radii (compared" << endl
       << "to the exact miniball).  Please consult the documentation" << endl
       << "for the details." << endl
       << "---------------------------------------------------------" << endl
       << "Running testsuite (this may take a while)..." << endl
       << "READ THE FINAL 'TEST RESULTS' line(s)..." << endl 
       << "---------------------------------------------------------" << endl
       << endl;

  mainTest<RunFlag>();

  cout << "---------------------------------------------------------" << endl
       << "TEST RESULTS" << endl
       << "---------------------------------------------------------" << endl
       << "Number of miniball computations: " << nr << endl
       << "Minimal relative error in the radius: " << mine << endl
       << "Maximal relative error in the radius: " << maxe << endl
       << "Average relative error in the radius: " << avge/nr << endl
       << "---------------------------------------------------------" << endl
       << endl;
}
