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
  Point_d(int dim, const Origin &o) : RPoint_d(dim, o)
    {}
  Point_d(const Self &p) : RPoint_d((RPoint_d&)p)
    {}
  Point_d(const RPoint_d &p) : RPoint_d(p)
    {}
  template <class InputIterator>
  Point_d (int dim, InputIterator first, InputIterator last)
        : RPoint_d (dim, first, last)
    {}
  Self& operator=(const Self& p)
    { RPoint_d::operator=(p); return *this; }
  bool operator==(const Self& p) const
    { return RPoint_d::operator==(p); }
  bool operator!=(const Self& p) const
    { return !(*this == p); }
  int id() const
    { return (int) PTR; }
  RT homogeneous(int i) const
    { return RPoint_d::homogeneous(i); }
  FT cartesian(int i) const
    { return RPoint_d::cartesian(i); }
  FT operator[](int i) const
    { return RPoint_d::operator[](i); }
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

