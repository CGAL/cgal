// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve.Bronnimann@sophia.inria.fr

#ifndef CGAL_PLANE_D_H
#define CGAL_PLANE_D_H

#ifndef CGAL_REP_CLASS_DEFINED
#error  no representation class defined
#endif  // CGAL_REP_CLASS_DEFINED

CGAL_BEGIN_NAMESPACE

template <class _R>
class Plane_d : public _R::Plane_d_base
{
public:
  typedef _R                                    R;
  typedef typename R::FT                        FT;
  typedef typename R::RT                        RT;
  typedef const RT*                             const_iterator ;
  typedef RT*                                   iterator ;
  typedef CGAL::Plane_d<R>                      Self;
  typedef typename R::Plane_d_base              RPlane_d;
  typedef typename R::Point_d_base              Point_d;
  typedef typename R::Direction_d_base          Direction_d;
  typedef typename R::Vector_d_base             Vector_d;
  // typedef typename R::Line_d_base               Line_d;

  Plane_d(int d = 0) : RPlane_d(d)
    {}
  Plane_d(const Self &h) : RPlane_d((const RPlane_d&)h)
    {}
  Plane_d(const RPlane_d &h) : RPlane_d(h)
    {}
  Plane_d(const Point_d &p, const Direction_d &d) : RPlane_d(p,d)
    {}
  Plane_d(const Point_d &p, const Vector_d &v) : RPlane_d(p,v)
    {}
  Plane_d(int dim, const Point_d &p, const Point_d &q, const Point_d &r)
    {
      CGAL_kernel_precondition( dim == 3);
      const Point_d e[3] = { p, q, r }; *this = RPlane_d(e+0,e+3);
    }
  Plane_d(int dim, const Origin &p, const Point_d &q, const Point_d &r)
    {
      CGAL_kernel_precondition( dim == 3);
      const Point_d e[3] = { Point_d(3,p), q, r }; *this = RPlane_d(e+0,e+3);
    }
  Plane_d(int dim, const RT &a, const RT &b, const RT &c, const RT &d)
    {
      CGAL_kernel_precondition( dim == 3);
      const RT e[4] = { a, b, c, d }; *this = RPlane_d(dim,e+0,e+4);
    }
  template < class InputIterator >
  Plane_d(const int d, const InputIterator &first, const InputIterator &last)
        : RPlane_d(d,first,last)
    {}
  template < class InputIterator >
  Plane_d(const InputIterator &first, const InputIterator &last,
          const Self &o, Oriented_side side = POSITIVE)
	: RPlane_d(first,last,o,side)
    {}

  Self &operator=(const Self &h)
    { RPlane_d::operator=(h); return *this; }
  bool operator==(const Self &h) const
    { return RPlane_d::operator==(h); }
  bool operator!=(const Self &h) const
    { return RPlane_d::operator!=(h); }

  RT operator[](int i) const
    { return RPlane_d::operator[](i); }
  RT a() const
    { return RPlane_d::operator[](0); }
  RT b() const
    { return RPlane_d::operator[](1); }
  RT c() const
    { return RPlane_d::operator[](2); }
  RT d() const
    { return RPlane_d::operator[](3); }

  long id() const
    { return (long)PTR; }
  // Line_d perpendicular_line(const Point_d &p) const
  //   { return RPlane_d::perpendicular_line(p); }
  Self opposite() const
    { return RPlane_d::opposite(); }
  Point_d point() const
    { return RPlane_d::point(); }
  Point_d projection(const Point_d &p) const
    { return RPlane_d::projection(p); }
  Vector_d orthogonal_vector() const
    { return RPlane_d::orthogonal_vector(); }
  Direction_d orthogonal_direction() const
    { return RPlane_d::orthogonal_direction(); }
  Vector_d base(const int i) const
    { return RPlane_d::base(i); }

  Point_d to_plane_basis(const Point_d &p) const
    { return RPlane_d::to_plane_basis(p); }
  // Self transform(const Aff_transformation_d &t) const;
  //  { return RPlane_d::transform(); }
  Oriented_side  oriented_side(const Point_d &p) const
    { return RPlane_d::oriented_side(p); }
  bool           has_on_boundary(const Point_d &p) const
    { return RPlane_d::has_on_boundary(p); }
  // bool           has_on_boundary(const Line_d &l) const
  // { return RPlane_d::has_on_boundary(l); }
  bool           has_on_positive_side(const Point_d &p) const
    { return RPlane_d::has_on_positive_side(p); }
  bool           has_on_negative_side(const Point_d &p) const
    { return RPlane_d::has_on_negative_side(p); }
  bool           has_on(const Point_d &p) const
    { return RPlane_d::has_on(p); }
  bool           is_degenerate() const
    { return RPlane_d::is_degenerate(); }

  int            dimension() const { return ptr->d; }
  const_iterator begin()     const { return ptr()->e; }
  const_iterator end()       const { return ptr()->e+dimension(); }

protected:
  iterator       begin()           { return ptr()->e; }
  iterator       end()             { return ptr()->e+dimension(); }

private:
  const _d_tuple<FT>*ptr()   const { return (const _d_tuple<FT>*)PTR; }
  _d_tuple<FT>*  ptr()             { return (_d_tuple<FT>*)PTR; }
  void           new_rep(const int dim, const FT*h);
  void           new_rep(const int k, const Self*h);
};

#ifndef NO_OSTREAM_INSERT_PLANE_D
template < class R >
std::ostream& 
operator<<(std::ostream& os, const Plane_d<R>& p)
{
  typedef typename  R::Plane_d_base    Plane;
  return os << (const Plane&)p;
}
#endif // NO_OSTREAM_INSERT_PLANE_D

#ifndef NO_ISTREAM_EXTRACT_PLANE_D
template < class R >
std::istream& 
operator>>(std::istream& is, Plane_d<R> &p)
{
  typedef typename  R::Plane_d_base    Plane;
  return is >> (Plane&)p;
}
#endif // NO_ISTREAM_EXTRACT_PLANE_D

CGAL_END_NAMESPACE

#endif  // CGAL_PLANE_D_H
