// Description: This program extensively tests the interface of the
// package Min_sphere_of_spheres_d.
//
// Compatibility: works with or without CGAL

#include <iostream>
#include <vector>
#include <CGAL/Min_sphere_of_spheres_d.h>
#include "Rational.h"

// The program will work with numbers of type FieldType.
// You can change FieldType to leda_real, for instance.
typedef Rational FieldType;

template<typename FT>
CGAL::Tag_true is_exact(const FT&) {
  return CGAL::Tag_true();
}

CGAL::Tag_false is_exact(const float) {
  return CGAL::Tag_false();
}

CGAL::Tag_false is_exact(const double) {
  return CGAL::Tag_false();
}

double uniform() {  // a (platform independent) random number generator
  static int lastNo = 230575L;
  const int a = 16807L, m = 2147483647L, q = 127773L, r = 2836L;
  int gamma = a * (lastNo % q) - r * (lastNo / q);
  if (gamma > 0)
    lastNo = gamma;
  else
    lastNo = gamma + m;
  return static_cast<double>(lastNo)/m;
}

template<typename FT,int D>
class Ball {
private:
  FT c[D];
  FT r;

public:
  Ball() {}

  template<typename InputIterator>
  Ball(InputIterator begin,const FT& radius) : r(radius) {
    std::copy(begin,begin+D,c);
  }

  const FT& radius() const {
    return r;
  }

  typedef const FT *Coordinate_iterator;
  Coordinate_iterator center() const {
    return c;
  }
};

template<int _Dim,typename _FT,typename _Sqrt,typename _Alg>
struct BallTraits {
  typedef Ball<_FT,_Dim> Sphere;
  static const int D = _Dim;
  typedef _Alg Algorithm;
  typedef _Sqrt Use_square_roots;
  typedef _FT FT;

  typedef typename Sphere::Coordinate_iterator Coordinate_iterator;
  static Coordinate_iterator begin(const Sphere& b) {
    return b.center();
  }
  static const FT& radius(const Sphere&b) {
    return b.radius();
  }
};

void checkCondition(bool c,const char *msg) {
  if (!c) {
    std::cout << msg << std::endl;
    std::cout.flush();
    exit(-1);
  }
}

void checkRelError(const double exact,
		   const double approx,
		   const double tol,
		   const char *msg) {
  using std::abs;
  if (CGAL_MINIBALL_NTS abs(exact) <= tol) 
    return;

  const double relerr = CGAL_MINIBALL_NTS abs(exact-approx) /
                        CGAL_MINIBALL_NTS abs(exact);
  if (relerr > tol) {
    std::cout << "Relative error " << relerr << " too large." << std::endl;
    checkCondition(false,msg);
  }
}

template<int D,typename Sphere,typename FT>
void compare(const FT& tol,const Sphere& m1,const Sphere& m2,
	     const CGAL::Tag_true is_exact) {
  typedef typename Sphere::Result Pair;
  typedef typename Sphere::Coordinate_iterator CIt;

  // check radii:
  const Pair r1 = m1.radius(), r2 = m2.radius();
  checkCondition(r1 == r2,"Radii do not match.");

  // check coordinates:
  CIt c1 = m1.center_begin();
  CIt c2 = m2.center_begin();
  for (int i=0; i<D; ++i, ++c1, ++c2)
    checkCondition(*c1 == *c2,"Center coordinates do not match.");
}

template<int D,typename Sphere,typename FT>
void compare(const FT& tol,const Sphere& m1,const Sphere& m2,
	     const CGAL::Tag_false is_exact) {
  typedef typename Sphere::Coordinate_iterator CIt;

  // check radii:
  const FT r1 = m1.radius(), r2 = m2.radius();
  checkRelError(r1,r2,tol,"Radii do not match.");

  // check coordinates:
  CIt c1 = m1.center_begin();
  CIt c2 = m2.center_begin();
  for (int i=0; i<D; ++i, ++c1, ++c2)
    checkRelError(*c1,*c2,tol,"Center coordinates do not match.");
}

template<int D,
	 typename FT,
         typename Sqrt,
         typename Alg>
void test(const int n,const FT& tol) {
  typedef BallTraits<D,FT,Sqrt,Alg> Traits;
  typedef CGAL::Min_sphere_of_spheres_d<Traits> Min_sphere;
  typedef typename Traits::Sphere Sphere;
  using std::cout;
  using std::endl;

  // create random balls:
  std::vector<Sphere> S;
  FT coord[D];
  for (int i=0; i<n; ++i) {
    for (int j=0; j<D; ++j)
      coord[j] = static_cast<FT>(uniform());
    S.push_back(Sphere(coord,static_cast<FT>(uniform())));
  }
  
  cout << " constructor with balls:" << endl;
  Min_sphere  ms1(S.begin(),S.end());
  ms1.update();
  checkCondition(ms1.is_valid(),"Minsphere not valid.");

  cout << " default constructor and set():" << endl;
  Min_sphere  ms2;
  ms2.set(S.begin(),S.end());
  checkCondition(ms2.is_valid(),"Minsphere not valid.");
  compare<D,Min_sphere,FT>(tol,ms1,ms2,is_exact(tol));

  cout << " default constructor and insert():" << endl;
  Min_sphere  ms3;
  ms3.insert(S.begin(),S.end());
  ms3.update();
  checkCondition(ms3.is_valid(),"Minsphere not valid.");
  compare<D,Min_sphere,FT>(tol,ms1,ms3,is_exact(tol));

  cout << " default constructor and multiple insert()'s:" << endl;
  Min_sphere  ms4;
  for (int i=0; i<S.size(); ++i) {
    ms4.insert(S[i]);
    ms4.update();
  }
  checkCondition(ms4.is_valid(),"Minsphere not valid.");
  compare<D,Min_sphere,FT>(tol,ms1,ms4,is_exact(tol));

  cout << " clearing and multiple insert()'s:" << endl;
  ms2.clear();
  for (int i=0; i<S.size(); ++i) {
    ms2.insert(S[i]);
    ms2.update();
  }
  checkCondition(ms2.is_valid(),"Minsphere not valid.");
  compare<D,Min_sphere,FT>(tol,ms1,ms2,is_exact(tol));
}

int main () {
  using std::cout;
  using std::endl;
  const int N = 20;

  const FieldType T = 0;     // tolerance for exact computation
  const double t = 1.0e-10;  // tolerance for double computation
  const float tf = 5.0e-7f;  // tolerance for float computation

  cout << "test exact LP-algorithm with sqrts:" << endl;
  test<2,FieldType,CGAL::Tag_true,CGAL::LP_algorithm>(N,T);
  test<3,FieldType,CGAL::Tag_true,CGAL::LP_algorithm>(N,T);
  test<5,FieldType,CGAL::Tag_true,CGAL::LP_algorithm>(N,T);

  cout << "test exact LP-algorithm without sqrts:" << endl;
  test<2,FieldType,CGAL::Tag_false,CGAL::LP_algorithm>(N,T);
  test<3,FieldType,CGAL::Tag_false,CGAL::LP_algorithm>(N,T);
  test<5,FieldType,CGAL::Tag_false,CGAL::LP_algorithm>(N,T);

  cout << "test exact farthest-first heuristic with sqrts:" << endl;
  test<2,FieldType,CGAL::Tag_true,CGAL::Farthest_first_heuristic>(N,T);
  test<3,FieldType,CGAL::Tag_true,CGAL::Farthest_first_heuristic>(N,T);
  test<5,FieldType,CGAL::Tag_true,CGAL::Farthest_first_heuristic>(N,T);

  cout << "test exact farthest-first heuristic without sqrts:" << endl;
  test<2,FieldType,CGAL::Tag_false,CGAL::Farthest_first_heuristic>(N,T);
  test<3,FieldType,CGAL::Tag_false,CGAL::Farthest_first_heuristic>(N,T);
  test<5,FieldType,CGAL::Tag_false,CGAL::Farthest_first_heuristic>(N,T);

  cout << "test double farthest-first heuristic with sqrts:" << endl;
  test<2,double,CGAL::Tag_true,CGAL::Farthest_first_heuristic>(N,t);
  test<3,double,CGAL::Tag_true,CGAL::Farthest_first_heuristic>(N,t);
  test<5,double,CGAL::Tag_true,CGAL::Farthest_first_heuristic>(N,t);

  cout << "test double farthest-first heuristic without sqrts:" << endl;
  test<2,double,CGAL::Tag_false,CGAL::Farthest_first_heuristic>(N,t);
  test<3,double,CGAL::Tag_false,CGAL::Farthest_first_heuristic>(N,t);
  test<5,double,CGAL::Tag_false,CGAL::Farthest_first_heuristic>(N,t);
    
  cout << "test float farthest-first heuristic with sqrts:" << endl;
  test<2,float,CGAL::Tag_true,CGAL::Farthest_first_heuristic>(N,tf);
  test<3,float,CGAL::Tag_true,CGAL::Farthest_first_heuristic>(N,tf);
  test<5,float,CGAL::Tag_true,CGAL::Farthest_first_heuristic>(N,tf);

  cout << "test float farthest-first heuristic without sqrts:" << endl;
  test<2,float,CGAL::Tag_false,CGAL::Farthest_first_heuristic>(N,tf);
  test<3,float,CGAL::Tag_false,CGAL::Farthest_first_heuristic>(N,tf);
  test<5,float,CGAL::Tag_false,CGAL::Farthest_first_heuristic>(N,tf);

  cout << endl
       << "Notice: Because the LP-algorithm is not good when used with" << endl
       << "double/float arithmetic, we do not test these combinations." << endl
       << endl
       << "Done." << endl;
}
