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
template < typename T >
struct IsVector_2 {
  enum { value = false };  // By default, a type is not a 2D vector.
};

template <>
struct IsVector_2 <My_vector> {
  enum { value = true };   // But My_vector is a 2D vector.
};

// Traits class determining if a type is a 2D point.
template < typename T >
struct IsPoint_2 {
  enum { value = false };  // By default, a type is not a 2D point.
};

template <>
struct IsPoint_2 <My_point> {
  enum { value = true };   // But My_point is a 2D point.
};

// Enable_if : tool for SFINAE
template < typename T, typename, bool = T::value >
struct Enable_if;

template < typename T, typename R >
struct Enable_if <T, R, true> {
  typedef R type;
};

template < typename Vector_2 >
typename Enable_if< IsVector_2<Vector_2>, Vector_2 >::type
operator-(Vector_2 const &v, Vector_2 const &w) {
  return Vector_2(v.x() - w.x(), v.y() - w.y());
}

template < typename Point_2 >
typename Enable_if< IsPoint_2<Point_2>, typename Point_2::Vector_2 >::type
operator-(Point_2 const &p, Point_2 const &q) {
  typedef typename Point_2::Vector_2 Vector_2;
  return Vector_2(p.x() - q.x(), p.y() - q.y());
}

int main() {
  My_vector v(1,2), w(3,4);
  My_vector z = v - w;   // OK
  (void) z;

  My_point p(1,2), q(3,4);
  // My_point s = p + q; // error :
                         // no match for `My_point& + My_point&' operator
  My_vector r = p - q;   // OK
  (void) r;

  return 0;
}
