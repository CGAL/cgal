#ifndef CGAL_CARTESIAN_REDEFINE_NAMES_2_H
#define CGAL_CTAG
#endif

#ifndef CGAL_CARTESIAN_LINE_2_C
#define CGAL_CARTESIAN_LINE_2_C

CGAL_BEGIN_NAMESPACE

template < class R >
inline
_Threetuple<typename R::FT>*
LineC2<R CGAL_CTAG>::ptr() const
{
  return (_Threetuple<FT>*)PTR;
}

template < class R >
CGAL_KERNEL_INLINE
void
LineC2<R CGAL_CTAG>::new_rep(const LineC2<R CGAL_CTAG>::Point_2 &p,
                             const LineC2<R CGAL_CTAG>::Point_2 &q)
{
  PTR = new _Threetuple<FT> (p.y() - q.y(),
                             q.x() - p.x(),
                             p.x()*q.y() - p.y()*q.x());
}

template < class R >
CGAL_KERNEL_INLINE
void
LineC2<R CGAL_CTAG>::new_rep(const typename R::FT &a,
                             const typename R::FT &b,
			     const typename R::FT &c)
{
  PTR = new _Threetuple<FT> (a, b, c);
}

template < class R >
inline
LineC2<R CGAL_CTAG>::LineC2()
{
  PTR = new _Threetuple<FT>;
}

template < class R >
inline
LineC2<R CGAL_CTAG>::LineC2(const LineC2<R CGAL_CTAG>  &l)
  : Handle((Handle&)l)
{}

template < class R >
inline
LineC2<R CGAL_CTAG>::LineC2(const LineC2<R CGAL_CTAG>::Point_2 &p,
                            const LineC2<R CGAL_CTAG>::Point_2 &q)
{
  new_rep(p, q);
}

template < class R >
inline
LineC2<R CGAL_CTAG>::LineC2(const typename R::FT &a,
                            const typename R::FT &b,
			    const typename R::FT &c)
{
  new_rep(a, b, c);
}

template < class R >
inline
LineC2<R CGAL_CTAG>::LineC2(const LineC2<R CGAL_CTAG>::Segment_2 &s)
{
  new_rep(s.start(), s.end());
}

template < class R >
inline
LineC2<R CGAL_CTAG>::LineC2(const LineC2<R CGAL_CTAG>::Ray_2 &r)
{
  new_rep(r.start(), r.second_point());
}

template < class R >
CGAL_KERNEL_INLINE
LineC2<R CGAL_CTAG>::LineC2(const LineC2<R CGAL_CTAG>::Point_2 &p,
                            const LineC2<R CGAL_CTAG>::Direction_2 &d)
{
  new_rep(-d.dy(), d.dx(), -d.dx()* p.y()  + d.dy() * p.x());
}


template < class R >
inline
LineC2<R CGAL_CTAG>::~LineC2()
{}


template < class R >
CGAL_KERNEL_INLINE
LineC2<R CGAL_CTAG> &
LineC2<R CGAL_CTAG>::operator=(const LineC2<R CGAL_CTAG> &l)
{

  Handle::operator=(l);
  return *this;
}

template < class R >
CGAL_KERNEL_MEDIUM_INLINE
bool LineC2<R CGAL_CTAG>::operator==(const LineC2<R CGAL_CTAG> &l) const
{
  if (  (a() * l.c() != l.a() * c())
      ||(b() * l.c() != l.b() * c()) )
      return false;
  int sc  = CGAL::sign(c());
  int slc = CGAL::sign(l.c());
  if ( sc == slc )
      return (sc == 0) ?  ( a()*l.b() ==  b()*l.a() )
                          && (CGAL::sign(a() ) == CGAL::sign( l.a() ))
                          && (CGAL::sign(b() ) == CGAL::sign( l.b() ))
                       : true;
  return false;
}

template < class R >
inline
bool LineC2<R CGAL_CTAG>::operator!=(const LineC2<R CGAL_CTAG> &l) const
{
  return !(*this == l);
}

template < class R >
inline
int LineC2<R CGAL_CTAG>::id() const
{
  return (int) PTR;
}

template < class R >
inline
typename R::FT LineC2<R CGAL_CTAG>::a() const
{
  return ptr()->e0;
}

template < class R >
inline
typename R::FT LineC2<R CGAL_CTAG>::b() const
{
  return ptr()->e1;
}

template < class R >
inline
typename R::FT LineC2<R CGAL_CTAG>::c() const
{
  return ptr()->e2;
}

template < class R >
CGAL_KERNEL_INLINE
typename R::FT
LineC2<R CGAL_CTAG>::x_at_y(const typename R::FT &y) const
{
  CGAL_kernel_precondition_msg( (a() != FT(0)),
  "Line::x_at_y(const typename R::FT &y) is undefined for horizontal line" );
  return ( -b()*y - c() ) / a();
}

template < class R >
CGAL_KERNEL_INLINE
typename R::FT
LineC2<R CGAL_CTAG>::y_at_x(const typename R::FT &x) const
{
  CGAL_kernel_precondition_msg( (b() != FT(0)),
  "Line::x_at_y(const typename R::FT &y) is undefined for vertical line");
  return ( -a()*x - c() ) / b();
}

template < class R >
inline
LineC2<R CGAL_CTAG>
LineC2<R CGAL_CTAG>::perpendicular(const LineC2<R CGAL_CTAG>::Point_2 &p) const
{
  return LineC2<R CGAL_CTAG>( -b() , a() , b() * p.x() - a() * p.y()  );
}


template < class R >
inline
LineC2<R CGAL_CTAG>
LineC2<R CGAL_CTAG>::opposite() const
{
  return LineC2<R CGAL_CTAG>( -a(), -b(), -c() );
}

template < class R >
CGAL_KERNEL_MEDIUM_INLINE
LineC2<R CGAL_CTAG>::Point_2
LineC2<R CGAL_CTAG>::point(int i) const
{
  if (i == 0)
    return is_vertical() ? Point_2( (-b()-c())/a(), FT(1) )
                         : Point_2( FT(1), -(a()+c())/b());
  if (i == 1)
    return is_vertical() ? Point_2( (-b()-c())/a() + b(), FT(1) - a() )
                         : Point_2( FT(1) + b(), -(a()+c())/b() - a() );
  // we add i times the direction
  if (is_vertical())
    return Point_2( (-b()-c())/a() + FT(i)*b(), FT(1) - FT(i)*a() );
  return Point_2( FT(1) + FT(i)*b(), -(a()+c())/b() - FT(i)*a() );
}

template < class R >
CGAL_KERNEL_INLINE
LineC2<R CGAL_CTAG>::Point_2
LineC2<R CGAL_CTAG>::point() const
{
  return is_vertical() ? Point_2( (-b()-c())/a(), FT(1) )
                       : Point_2( FT(1), -(a()+c())/b());
}

template < class R >
CGAL_KERNEL_MEDIUM_INLINE
LineC2<R CGAL_CTAG>::Point_2
LineC2<R CGAL_CTAG>::projection(const LineC2<R CGAL_CTAG>::Point_2 &p) const
{
  if (is_horizontal()) return Point_2(p.x(), -c()/b());
  if (is_vertical())   return Point_2( -c()/a(), p.y());

  FT ab = a()/b(), ba = b()/a(), ca = c()/a();
  FT y = ( -p.x() + ab*p.y() - ca ) / ( ba + ab );
  return Point_2(-ba * y - ca, y);
}

template < class R >
inline
LineC2<R CGAL_CTAG>::Direction_2
LineC2<R CGAL_CTAG>::direction() const
{
  return LineC2<R CGAL_CTAG>::Direction_2( b(), -a() );
}

template < class R >
CGAL_KERNEL_INLINE
Oriented_side
LineC2<R CGAL_CTAG>::
oriented_side(const LineC2<R CGAL_CTAG>::Point_2 &p) const
{
  return Oriented_side(CGAL::sign(a()*p.x() + b()*p.y() + c()));
}

template < class R >
inline
bool
LineC2<R CGAL_CTAG>::
has_on_boundary(const LineC2<R CGAL_CTAG>::Point_2 &p) const
{
  return (a()*p.x() + b()*p.y() + c()) == FT(0);
}

template < class R >
inline
bool
LineC2<R CGAL_CTAG>::
has_on_positive_side(const LineC2<R CGAL_CTAG>::Point_2 &p) const
{
  return (a()*p.x() + b()*p.y() + c()) >  FT(0);
}

template < class R >
CGAL_KERNEL_INLINE
bool
LineC2<R CGAL_CTAG>::
has_on_negative_side(const LineC2<R CGAL_CTAG>::Point_2 &p) const
{
  return (a()*p.x() + b()*p.y() + c()) <  FT(0);
}

template < class R >
inline
bool LineC2<R CGAL_CTAG>::is_horizontal() const
{
  return a() == LineC2<R CGAL_CTAG>::FT(0) ;
}

template < class R >
inline
bool LineC2<R CGAL_CTAG>::is_vertical() const
{
  return b() == FT(0) ;
}

template < class R >
inline
bool LineC2<R CGAL_CTAG>::is_degenerate() const
{
  return (a() == FT(0)) && (b() == FT(0)) ;
}

template < class R >
inline
LineC2<R CGAL_CTAG>
LineC2<R CGAL_CTAG>::
transform(const LineC2<R CGAL_CTAG>::Aff_transformation_2 &t) const
{
  return LineC2<R CGAL_CTAG>( t.transform(point(0) ), t.transform(direction() ));
}


#ifndef CGAL_NO_OSTREAM_INSERT_LINEC2
template < class R >
std::ostream &operator<<(std::ostream &os, const LineC2<R CGAL_CTAG> &l)
{

    switch(os.iword(IO::mode)) {
    case IO::ASCII :
        return os << l.a() << ' ' << l.b() << ' ' << l.c();
    case IO::BINARY :
        write(os, l.a());
        write(os, l.b());
        write(os, l.c());
        return os;
    default:
        return os << "LineC2(" << l.a() << ", " << l.b() << ", " << l.c() <<')';
    }
}
#endif // CGAL_NO_OSTREAM_INSERT_LINEC2

#ifndef CGAL_NO_ISTREAM_EXTRACT_LINEC2
template < class R >
std::istream &operator>>(std::istream &is, LineC2<R CGAL_CTAG> &p)
{
    typename R::FT a, b, c;
    switch(is.iword(IO::mode)) {
    case IO::ASCII :
        is >> a >> b >> c;
        break;
    case IO::BINARY :
        read(is, a);
        read(is, b);
        read(is, c);
        break;
    default:
        cerr << "" << std::endl;
        cerr << "Stream must be in ascii or binary mode" << std::endl;
        break;
    }
    p = LineC2<R CGAL_CTAG>(a, b, c);
    return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_LINEC2

CGAL_END_NAMESPACE

#endif
