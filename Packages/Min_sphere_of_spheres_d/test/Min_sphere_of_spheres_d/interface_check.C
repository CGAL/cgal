// Program testing the package extensively.
#include <CGAL/Cartesian_d.h>
#include <CGAL/Random.h>
#include <CGAL/Min_sphere_of_spheres_d.h>
#include <CGAL/Gmpq.h>
#include <iostream>
#include <vector>

template<typename FT>
CGAL::Tag_true is_exact(const FT&) {
  return CGAL::Tag_true();
}

CGAL::Tag_false is_exact(const double) {
  return CGAL::Tag_false();
}

void checkCondition(bool c,const char *msg) {
  if (!c) {
    std::cout << msg << endl;
    std::cout.flush();
    exit(-1);
  }
}

void checkRelError(const double exact,
		   const double approx,
		   const double tol,
		   const char *msg) {
  if (CGAL_NTS abs(exact) <= tol) 
    return;

  const double relerr = CGAL_NTS abs(exact-approx) / CGAL_NTS abs(exact);
  if (relerr > tol) {
    std::cout << "Relative error " << relerr << " too large." << std::endl;
    checkCondition(false,msg);
  }
}

template<typename Sphere,typename FT>
void compare(const FT& tol,const Sphere& m1,const Sphere& m2,
	     const CGAL::Tag_true is_exact) {
  typedef typename Sphere::Result Pair;
  typedef typename Sphere::Center_coordinate_iterator CIt;

  // check radii:
  const Pair r1 = m1.radius(), r2 = m2.radius();
  checkCondition(r1.first == r2.first && r1.second == r2.second,
		 "Radii do not match.");

  // check coordinates:
  CIt c1 = m1.center_begin(), c1end = m1.center_end();
  CIt c2 = m2.center_begin(), c2end = m2.center_end();
  for (; c1 != c1end; ++c1, ++c2) {
    checkCondition(c2 != c2end,"Center iterator miss-sized.");
    checkCondition((*c1).first == (*c2).first &&
		   (*c1).second == (*c2).second,
		   "Center coordinates do not match.");
  }
  checkCondition(c2 == c2end,"Center iterator miss-sized.");
}

template<typename Sphere,typename FT>
void compare(const FT& tol,const Sphere& m1,const Sphere& m2,
	     const CGAL::Tag_false is_exact) {
  typedef typename Sphere::Center_coordinate_iterator CIt;

  // check radii:
  const double r1 = m1.radius(), r2 = m2.radius();
  checkRelError(r1,r2,tol,"Radii do not match.");

  // check coordinates:
  CIt c1 = m1.center_begin(), c1end = m1.center_end();
  CIt c2 = m2.center_begin(), c2end = m2.center_end();
  for (; c1 != c1end; ++c1, ++c2) {
    checkCondition(c2 != c2end,"Center iterator miss-sized.");
    checkRelError(*c1,*c2,tol,"Center coordinates do not match.");
  }
  checkCondition(c2 == c2end,"Center iterator miss-sized.");
}

template<typename FT,
         typename Sqrt,
         typename Alg>
void test(const int d,const int n,const FT& tol) {
  typedef CGAL::Cartesian_d<FT> K;
  typedef CGAL::Min_sphere_of_spheres_d_traits_d<K,FT,Sqrt,Alg> Traits;
  typedef CGAL::Min_sphere_of_spheres_d<Traits> Min_sphere;
  typedef typename K::Point_d Point;
  typedef typename Traits::Sphere Sphere;
  using std::cout;
  using std::endl;

  // create random balls:
  std::vector<Sphere> S;
  std::vector<FT> coord;
  CGAL::Random  r;
  const int Low = 0, High = 10000;
  for (int i=0; i<n; ++i) {
    coord.clear();
    for (int j=0; j<d; ++j)
      coord.push_back(static_cast<FT>(r.get_double(Low,High)));
    Point p(d,coord.begin(),coord.end());
    S.push_back(Sphere(p,static_cast<FT>(r.get_double(Low,High))));
  }
  
  cout << " constructor with balls:" << endl;
  Min_sphere  ms1(S.begin(),S.end());
  checkCondition(ms1.is_valid(true,0),"Minsphere not valid.");

  cout << " default constructor and set():" << endl;
  Min_sphere  ms2;
  ms2.set(S.begin(),S.end());
  checkCondition(ms2.is_valid(true,0),"Minsphere not valid.");
  compare<Min_sphere,FT>(tol,ms1,ms2,is_exact(tol));

  cout << " default constructor and insert():" << endl;
  Min_sphere  ms3;
  ms3.insert(S.begin(),S.end());
  checkCondition(ms3.is_valid(true,0),"Minsphere not valid.");
  compare<Min_sphere,FT>(tol,ms1,ms3,is_exact(tol));
  
  cout << " default constructor and multiple insert()'s:" << endl;
  Min_sphere  ms4;
  for (int i=0; i<S.size(); ++i)
    ms4.insert(S[i]);
  checkCondition(ms4.is_valid(true,0),"Minsphere not valid.");
  compare<Min_sphere,FT>(tol,ms1,ms4,is_exact(tol));
  
  cout << " clearing and multiple insert()'s:" << endl;
  ms2.clear();
  for (int i=0; i<S.size(); ++i)
    ms2.insert(S[i]);
  checkCondition(ms2.is_valid(true,0),"Minsphere not valid.");
  compare<Min_sphere,FT>(tol,ms1,ms2,is_exact(tol));
}

int main () {
  using std::cout;
  using std::endl;
  const int N = 20;

  const CGAL::Gmpq T = 0;  // tolerance for exact computation
  const double t = 1e-10;  // tolerance for double computation

  cout << "test exact LP-algorithm with sqrts:" << endl;
  for (int d=2; d<6; ++d)
    test<CGAL::Gmpq,CGAL::Tag_true,CGAL::LP_algorithm>(d,N,T);

  cout << "test exact LP-algorithm without sqrts:" << endl;
  for (int d=2; d<6; ++d)
    test<CGAL::Gmpq,CGAL::Tag_false,CGAL::LP_algorithm>(d,N,T);

  cout << "test exact farthest-first heuristic with sqrts:" << endl;
  for (int d=2; d<6; ++d)
    test<CGAL::Gmpq,CGAL::Tag_true,CGAL::Farthest_first_heuristic>(d,N,T);

  cout << "test exact farthest-first heuristic without sqrts:" << endl;
  for (int d=2; d<6; ++d)
    test<CGAL::Gmpq,CGAL::Tag_false,CGAL::Farthest_first_heuristic>(d,N,T);

  cout << "test double farthest-first heuristic with sqrts:" << endl;
  for (int d=2; d<6; ++d)
    test<double,CGAL::Tag_true,CGAL::Farthest_first_heuristic>(d,N,t);

  cout << "test double farthest-first heuristic without sqrts:" << endl;
  for (int d=2; d<6; ++d)
    test<double,CGAL::Tag_false,CGAL::Farthest_first_heuristic>(d,N,t);

  cout << endl
       << "Notice: Because the LP-algorithm is not good when used with" << endl
       << "double arithmetic, we do not test this combination." << endl
       << endl
       << "Done." << endl;
}

