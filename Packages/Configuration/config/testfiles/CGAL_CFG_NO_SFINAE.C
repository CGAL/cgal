// ======================================================================
//
// Copyright (c) 2003 The CGAL Consortium
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 
//
// file          : config/testfiles/CGAL_CFG_NO_SFINAE.C
// source        :
// author(s)     : Sylvain Pion
//
// coordinator   : Utrecht University
//
// ======================================================================

// ---------------------------------------------------------------------
// A short test program to evaluate a C++ compiler.
// This program is used by cgal_configure.
// The following documentation will be pasted in the generated configfile.
// ---------------------------------------------------------------------

//| This flag is set if the compiler doesn't support the SFINAE principle
//| (Substitution Failure Is Not An Error), which we eventually plan to use
//| for the next version of the kernel design.


class My_vector {   // 2D vector type
  int _x, _y;
public:

  My_vector(int x, int y) : _x(x), _y(y) {}

  int x() const { return _x; }
  int y() const { return _y; }
};

class My_point {   // 2D point type
  int _x, _y;
public:

  typedef My_vector Vector_2; // type returned by p-q;

  My_point(int x, int y) : _x(x), _y(y) {}

  int x() const { return _x; }
  int y() const { return _y; }
};

// Traits class determining if a type is a 2D vector.
template < typename > struct IsVector_2            { enum { value = false }; };
template <>           struct IsVector_2<My_vector> { enum { value = true  }; };

// Traits class determining if a type is a 2D point.
template < typename > struct IsPoint_2             { enum { value = false }; };
template <>           struct IsPoint_2 <My_point>  { enum { value = true  }; };


// Enable_if : tool for SFINAE
template < bool, typename >
struct Enable_if {};

template < typename T >
struct Enable_if <true, T> { typedef T type; };

// Workaround for making the following signature slightly different
// on the arguments (otherwise g++ 3.2 barfs on it).
template < typename T >
struct Same { typedef T type; };
#define SAME(T) typename Same<T>::type
// #define SAME(T) T


template < typename Vector_2 >
typename Enable_if< IsVector_2<Vector_2>::value, Vector_2 >::type
operator-(Vector_2 const &v, Vector_2 const &w) {
  return Vector_2(v.x() - w.x(), v.y() - w.y());
}

template < typename Vector_2 >
typename Enable_if< IsVector_2<Vector_2>::value, Vector_2 >::type
operator+(Vector_2 const &v, Vector_2 const &w) {
  return Vector_2(v.x() + w.x(), v.y() + w.y());
}


template < typename Point_2 >
typename Enable_if< IsPoint_2<Point_2>::value, Point_2 >::type::Vector_2
operator-(Point_2 const &p, SAME(Point_2) const &q) {
  typedef typename Point_2::Vector_2 Vector_2;
  return Vector_2(p.x() - q.x(), p.y() - q.y());
}

int main() {
  My_vector v(1,2), w(3,4);
  My_vector z = v - w;   // OK
  My_vector y = v + w;   // OK

  My_point p(1,2), q(3,4);
  My_vector r = p - q;   // OK
  // My_point s = p + q; // error :
                         // no match for `My_point& + My_point&' operator

  (void) z; (void) y; (void) r;
  return 0;
}
