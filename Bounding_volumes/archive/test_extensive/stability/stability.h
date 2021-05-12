// (Refer to file stability.cpp or stability-exact.cpp for more information.)

#ifndef MINIBALL_STABILITY
#define MINIBALL_STABILITY

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <string>
#include <fstream>
#include <CGAL/Min_sphere_of_spheres_d.h>

// A functor to convert doubles to some other type T:
template<typename T>
struct fromDouble {
  typedef T Result;
  const T operator()(const double x) {
    return static_cast<T>(x);
  }
};

// A class to build an iterator from a functor and an iterator.
// Assumes that Functor::Result is the return type of
// Functor::operator().
template<typename Iterator,typename Functor>
class FunctorIterator {
private:
  Iterator it;
  Functor f;

public:
  typedef typename Functor::Result Result;

  FunctorIterator(Iterator it) : it(it) {}
  FunctorIterator(Iterator it,const Functor& f) : it(it), f(f) {}

  const Result operator*() {
    return f(*it);
  }

  FunctorIterator& operator++() {
    ++it;
    return *this;
  }

  FunctorIterator operator++(int) {
    FunctorIterator old(*this);
    ++(*this);
    return old;
  }

  bool operator!=(const FunctorIterator& i) {
    return it != i.it;
  }
};

// We start with the definition of a sphere class parametrized by
// the number type and the dimension:
template<typename Number,int d>
class Ball {
private: // representation:
  Number c[d];
  Number r;

public: // construction:
  // Construct a ball with radius r and center coordinates
  // from the range [first,first+d).
  template<typename InputIterator>
  Ball(InputIterator first,const Number& r) : r(r) {
    assert(r >= 0);
    std::copy(first,first+d,c);
  }

public: // access:
  // Returns the radius.
  const Number& radius() const {
    return r;
  }
  // Returns an iterator to the first coordinate of the center.
  const Number *center() const {
    return c;
  }

  // Returns the i-th coordinate for 0<=i<d.
  const Number& operator[](int i) const {
    assert(0<=i && i<d);
    return c[i];
  }
};

// The following is a traits class to make the Miniball package work
// with balls of the type Ball<d> defined above:
template<int d,
         typename StorageNumberType,
         typename NumberType,
         typename Sqrt>
struct BallTraits {
  // constants:
  static const int D = d;

  // types:
  typedef Ball<StorageNumberType,d> Sphere;
  typedef NumberType FT;
  typedef CGAL::Default_algorithm Algorithm;
  typedef Sqrt Use_square_roots;
  typedef FunctorIterator< const StorageNumberType *,
                           fromDouble<NumberType> > Cartesian_const_iterator;

  // routines:
  inline static Cartesian_const_iterator
    center_cartesian_begin(const Sphere& b) {
    return Cartesian_const_iterator(b.center());
  }
  inline static FT radius(const Sphere& b) {
    return static_cast<NumberType>(b.radius());
  }
};

// A portable random number generator class:
class Uniform {
private:
  int lastNo;

public:
  Uniform() : lastNo(230575L) {}

  double next() {
    const int a = 16807L, m = 2147483647L, q = 127773L, r = 2836L;
    int gamma = a * (lastNo % q) - r * (lastNo / q);
    if (gamma > 0)
      lastNo = gamma;
    else
      lastNo = gamma + m;
    return static_cast<double>(lastNo)/m;
  }
};

enum Perturbation {
  None=0, Tiny=1, Small=2, Medium=3, Large=4
};

// The following routine reads a set of balls in R^d from a certain
// file, embeds them into R^D and perturbs their coordinates by a
// certain amount of pseudo-random noise.  Finally, each ball is
// duplicated copies times.  The resulting set of balls is returned in
// the variable balls.
template<int d,int D,typename Number>
void readConfiguration(const std::string& name, // name of the configuration
                       Perturbation perturb,    // amount of perturbation
                       double magnitude,        // magnitude of perturbation
                       int copies,              // number of copies
                       std::vector< Ball<Number,D> >& balls) {
  // clear set and open data file:
  balls.clear();
  std::ifstream data(name.c_str());

  // determine amount of permutation:
  Uniform rng;
  const double amounts[5] = { 0, 1e-50, 1e-30, 1e-10, 1e-3 };
  double amount = magnitude * amounts[static_cast<int>(perturb)];

  // read number of input balls:
  int n;
  data >> n;

  for (int i=0; i<n; ++i) {
    // read coordinate and radius of the i-th ball:
    Number r,c[d];
    for (int j=0; j<d; ++j)
      data >> c[j];
    data >> r;

    // insert duplicates:
    for (int k=0; k<copies; ++k) {
      // embed and perturb:
      Number ce[D];
      for (int j=0; j<d; ++j)
        ce[j] = c[j] + (rng.next()-0.5)*amount;
      for (int j=d; j<D; ++j)
        ce[j] = (rng.next()-0.5)*amount;

      // add to set:
      balls.push_back(Ball<Number,D>(ce,r));
    }
  }
}

// Converts a string to some type T.
template <typename T>
T strTo(const std::string& str) {
  std::stringstream strm(str);
  T t;
  strm >> t;
  return t;
}

// Converts a value to a string.
template <typename T>
std::string toStr(const T& t) {
  std::stringstream strm;
  strm << t;
  return strm.str();
}

struct RunFlag {};
struct BuildFlag {};

template<int d,int D>
void atomicTest(const int,const std::string&,Perturbation,
                const std::vector< Ball<double,D> >& S,RunFlag);
template<int d,int D>
void atomicTest(const int,const std::string&,Perturbation,
                const std::vector< Ball<double,D> >& S,BuildFlag);

// (See routine basicTest() below.)
int from = 0, to = 1000, curr = -1;

template<int d,int D,typename Flag>
void basicTest(const std::string& name,double magnitude,
               Perturbation perturb,int copies) {
  // decide whether to do the test or not:
  ++curr;
  if (curr < from || curr > to)
    return;

  // construct the set of balls:
  std::vector< Ball<double,D> > S;
  readConfiguration<d,D>(name,perturb,magnitude,copies,S);

  // perform test:
  std::cout << std::endl << "    Perturbation: ";
  switch (perturb) {
  case None:
    std::cout << "none" << std::endl;
    break;
  case Tiny:
    std::cout << "10^-50" << std::endl;
    break;
  case Small:
    std::cout << "10^-30" << std::endl;
    break;
  case Medium:
    std::cout << "10^-10" << std::endl;
    break;
  case Large:
    std::cout << "10^-3" << std::endl;
    break;
  }
  atomicTest<d,D>(curr,name,perturb,S,Flag());
}

template <int d,typename Flag>
void test(const std::string& f,double magnitude) {
  using std::cout;
  using std::endl;

  cout << f << endl;

  cout << "  Set (d):      ";
  basicTest<d,d,Flag>(f,magnitude,None,1);
  basicTest<d,d,Flag>(f,magnitude,Tiny,1);
  basicTest<d,d,Flag>(f,magnitude,Small,1);
  basicTest<d,d,Flag>(f,magnitude,Medium,1);
  basicTest<d,d,Flag>(f,magnitude,Large,1);
  cout << endl;

  cout << "  Multiset (d): ";
  basicTest<d,d,Flag>(f,magnitude,None,10);
  basicTest<d,d,Flag>(f,magnitude,Tiny,10);
  basicTest<d,d,Flag>(f,magnitude,Small,10);
  basicTest<d,d,Flag>(f,magnitude,Medium,10);
  basicTest<d,d,Flag>(f,magnitude,Large,10);
  cout << endl;

  cout << "  Set (D):      ";
  basicTest<d,2*d,Flag>(f,magnitude,None,1);
  basicTest<d,2*d,Flag>(f,magnitude,Tiny,1);
  basicTest<d,2*d,Flag>(f,magnitude,Small,1);
  basicTest<d,2*d,Flag>(f,magnitude,Medium,1);
  basicTest<d,2*d,Flag>(f,magnitude,Large,1);
  cout << endl;

  cout << "  Multiset (D): ";
  basicTest<d,2*d,Flag>(f,magnitude,None,10);
  basicTest<d,2*d,Flag>(f,magnitude,Tiny,10);
  basicTest<d,2*d,Flag>(f,magnitude,Small,10);
  basicTest<d,2*d,Flag>(f,magnitude,Medium,10);
  basicTest<d,2*d,Flag>(f,magnitude,Large,10);
  cout << endl;

  cout << endl;
}

template<typename Flag>
void mainTest() {
  test<2,Flag>("data/cocircular_points_small_radius_2.data",1e8);
  test<2,Flag>("data/cocircular_points_large_radius_2.data",1e9);

  test<3,Flag>("data/almost_cospherical_points_3.data",1);
  test<10,Flag>("data/almost_cospherical_points_10.data",1);

  test<3,Flag>("data/longitude_latitude_model_3.data",1);

  test<3,Flag>("data/random_points_3.data",1);
  test<5,Flag>("data/random_points_5.data",1);
  test<10,Flag>("data/random_points_10.data",1);

  test<10,Flag>("data/simplex_10.data",1);
  test<15,Flag>("data/simplex_15.data",1);

  test<10,Flag>("data/cube_10.data",1);
  test<12,Flag>("data/cube_12.data",1);

  test<2,Flag>("data/balls_on_boundary_2.data",1);
  test<3,Flag>("data/balls_on_boundary_3.data",1);
}

#endif // MINIBALL_STABILITY
