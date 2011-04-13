// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve.Bronnimann@sophia.inria.fr

#ifndef CGAL_POINT_D_H
#define CGAL_POINT_D_H

#ifndef CGAL_REP_CLASS_DEFINED
#error  no representation class defined
#endif  // CGAL_REP_CLASS_DEFINED

CGAL_BEGIN_NAMESPACE

template < class _R >
class Point_d : public _R::Point_d_base
{
public:
  typedef  _R   R;
  typedef typename R::RT              RT;
  typedef typename R::FT              FT;
  typedef CGAL::Point_d<R>            Self;
  typedef typename R::Point_d_base    RPoint_d;
  typedef const RT*                   const_iterator;

  Point_d()
    {}
  Point_d(const Self &p) : RPoint_d((const RPoint_d&)p)
    {}
  Point_d(const RPoint_d &p) : RPoint_d(p)
    {}
  Point_d(int dim, const Origin &o) : RPoint_d(dim, o)
    {}
  template <class InputIterator>
  Point_d(int dim, InputIterator first, InputIterator last)
        : RPoint_d (dim, first, last)
    {}
  Point_d(int dim, const RT &x, const RT &y, const RT &w)
    {
      CGAL_kernel_precondition( dim == 2 || dim == 3);
      RT e[3] = { x, y, w }; *this = RPoint_d(dim,e+0,e+3);
    }
  Point_d(int dim, const RT &x, const RT &y, const RT &z, const RT &w)
    {
      CGAL_kernel_precondition( dim == 3 || dim == 4);
      RT e[4] = { x, y, z, w }; *this = RPoint_d(dim,e+0,e+4);
    }
  Point_d(int dim,
          const RT &x, const RT &y, const RT &z, const RT &t, const RT &w)
    {
      CGAL_kernel_precondition( dim == 4);
      RT e[5] = { x, y, z, t, w }; *this = RPoint_d(dim,e+0,e+5);
    }
  Self& operator=(const Self& p)
    { RPoint_d::operator=(p); return *this; }
  bool operator==(const Self& p) const
    { return RPoint_d::operator==(p); }
  bool operator!=(const Self& p) const
    { return !(*this == p); }
  long id() const
    { return (long) PTR; }
  RT homogeneous(int i) const
    { return RPoint_d::homogeneous(i); }
  FT cartesian(int i) const
    { return RPoint_d::cartesian(i); }
  FT operator[](int i) const
    { return RPoint_d::operator[](i); }
  FT x() { return cartesian(0); }
  FT y() { return cartesian(1); }
  FT z() { return cartesian(2); }
  RT hx() { return homogeneous(0); }
  RT hy() { return homogeneous(1); }
  RT hz() { return homogeneous(2); }
  RT hw() { return homogeneous(3); }
  const_iterator begin() const
    { return RPoint_d::begin(); }
  const_iterator end() const
    { return RPoint_d::end(); }   
  int dimension() const
    { return RPoint_d::dimension(); }
};

#ifndef NO_OSTREAM_INSERT_POINT_D
template < class R >
std::ostream& 
operator<<(std::ostream& os, const Point_d<R>& p)
{
  typedef typename  R::Point_d_base    Point;
  return os << (const Point&)p;
}
#endif // NO_OSTREAM_INSERT_POINT_D

#ifndef NO_ISTREAM_EXTRACT_POINT_D
template < class R >
std::istream& 
operator>>(std::istream& is, Point_d<R> &p)
{
  typedef typename  R::Point_d_base    Point;
  return is >> (Point&)p;
}
#endif // NO_ISTREAM_EXTRACT_POINT_D
CGAL_END_NAMESPACE


#endif // CGAL_POINT_D_H

