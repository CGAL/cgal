#ifndef CGAL_CARTESIAN_REDEFINE_NAMES_2_H
#define CGAL_CTAG
#endif

#ifndef CGAL_CARTESIAN_CIRCLE_2_C
#define CGAL_CARTESIAN_CIRCLE_2_C

CGAL_BEGIN_NAMESPACE

template < class R >
inline
Circle_repC2<R> *CircleC2<R CGAL_CTAG>::ptr() const
{
  return (Circle_repC2<R>*)PTR;
}

template < class R >
CGAL_KERNEL_CTOR_INLINE
CircleC2<R CGAL_CTAG>::CircleC2()
{
  PTR = new Circle_repC2<R> ;
}

template < class R >
CGAL_KERNEL_CTOR_INLINE
CircleC2<R CGAL_CTAG>::CircleC2(const CircleC2<R CGAL_CTAG> &t)
  : Handle((Handle&)t)
{}

template < class R >
CGAL_KERNEL_CTOR_INLINE
CircleC2<R CGAL_CTAG>::CircleC2(const CircleC2<R CGAL_CTAG>::Point_2 &center,
                      const typename R::FT &squared_radius,
                      const Orientation &orient)
{
  CGAL_kernel_precondition( ( squared_radius >= FT(0) ) &&
                            ( orient    != COLLINEAR) );

  PTR = new Circle_repC2<R>(center, squared_radius, orient);
}

template < class R >
CGAL_KERNEL_CTOR_INLINE
CircleC2<R CGAL_CTAG>::CircleC2(const CircleC2<R CGAL_CTAG>::Point_2 &center,
                      const Orientation &orient)
{
  CGAL_kernel_precondition( ( orient    != COLLINEAR) );

  PTR = new Circle_repC2<R>(center, FT(0), orient);
}

template < class R >
CGAL_KERNEL_CTOR_MEDIUM_INLINE
CircleC2<R CGAL_CTAG>::CircleC2(const CircleC2<R CGAL_CTAG>::Point_2 &p,
                      const CircleC2<R CGAL_CTAG>::Point_2 &q,
                      const Orientation &orient)
{
  CGAL_kernel_precondition( orient != COLLINEAR);

  if ( p != q) {
    CircleC2<R CGAL_CTAG>::Point_2 center = midpoint(p,q);
    CircleC2<R CGAL_CTAG>::FT      squared_radius = squared_distance(p,center);

    PTR = new Circle_repC2<R>( center, squared_radius, orient);
  } else
    PTR = new Circle_repC2<R>( p, FT(0), orient);
}


template < class R >
CGAL_KERNEL_CTOR_MEDIUM_INLINE
CircleC2<R CGAL_CTAG>::CircleC2(const CircleC2<R CGAL_CTAG>::Point_2 &p,
                      const CircleC2<R CGAL_CTAG>::Point_2 &q,
                      const CircleC2<R CGAL_CTAG>::Point_2 &r)
{
  Orientation orient = CGAL::orientation(p,q,r);
  CGAL_kernel_precondition( orient != COLLINEAR);

  Point_2 center = circumcenter(p,q,r);
  FT      squared_radius = squared_distance(p,center);

  PTR = new Circle_repC2<R>(center, squared_radius, orient);
}


template < class R >
inline
CircleC2<R CGAL_CTAG>::~CircleC2()
{}


template < class R >
inline
CircleC2<R CGAL_CTAG> &CircleC2<R CGAL_CTAG>::operator=(const CircleC2<R CGAL_CTAG> &t)
{
  Handle::operator=(t);
  return *this;
}
template < class R >
CGAL_KERNEL_INLINE
bool CircleC2<R CGAL_CTAG>::operator==(const CircleC2<R CGAL_CTAG> &t) const
{
   return (center() == t.center()) &&
          (squared_radius() == t.squared_radius() &&
          orientation() == t.orientation());
}

template < class R >
inline
bool CircleC2<R CGAL_CTAG>::operator!=(const CircleC2<R CGAL_CTAG> &t) const
{
  return !(*this == t);
}

template < class R >
int CircleC2<R CGAL_CTAG>::id() const
{
  return (int)PTR;
}
template < class R >
inline
CircleC2<R CGAL_CTAG>::Point_2 CircleC2<R CGAL_CTAG>::center() const
{
 return ptr()->center;
}

template < class R >
inline
CircleC2<R CGAL_CTAG>::FT CircleC2<R CGAL_CTAG>::squared_radius() const
{
 return ptr()->squared_radius;
}

template < class R >
inline
Orientation CircleC2<R CGAL_CTAG>::orientation() const
{
 return ptr()->orient;
}


template < class R >
CGAL_KERNEL_MEDIUM_INLINE
Oriented_side CircleC2<R CGAL_CTAG>::oriented_side(const CircleC2<R CGAL_CTAG>::Point_2 &p) const
{
  return Oriented_side(bounded_side(p) * orientation());
}

template < class R >
CGAL_KERNEL_INLINE
Bounded_side CircleC2<R CGAL_CTAG>::bounded_side(const CircleC2<R CGAL_CTAG>::Point_2 &p) const
{
  return Bounded_side(CGAL::compare(squared_radius(),
                                    squared_distance(center(),p)));
}

template < class R >
inline
bool CircleC2<R CGAL_CTAG>::has_on_boundary(const CircleC2<R CGAL_CTAG>::Point_2 &p) const
{
  return squared_distance(center(),p) == squared_radius();
}

template < class R >
CGAL_KERNEL_INLINE
bool CircleC2<R CGAL_CTAG>::has_on_negative_side(const CircleC2<R CGAL_CTAG>::Point_2 &p) const
{
  if (orientation() == COUNTERCLOCKWISE) {
    return has_on_unbounded_side(p);
  }
  return has_on_bounded_side(p);
}

template < class R >
CGAL_KERNEL_INLINE
bool CircleC2<R CGAL_CTAG>::has_on_positive_side(const CircleC2<R CGAL_CTAG>::Point_2 &p) const
{
  if (orientation() == COUNTERCLOCKWISE) {
    return has_on_bounded_side(p);
  }
  return has_on_unbounded_side(p);
}

template < class R >
inline
bool CircleC2<R CGAL_CTAG>::has_on_bounded_side(const CircleC2<R CGAL_CTAG>::Point_2 &p) const
{
  return squared_distance(center(),p) < squared_radius();
}

template < class R >
inline
bool CircleC2<R CGAL_CTAG>::has_on_unbounded_side(const CircleC2<R CGAL_CTAG>::Point_2 &p) const
{
  return squared_distance(center(),p) > squared_radius();
}


template < class R >
inline
bool CircleC2<R CGAL_CTAG>::is_degenerate() const
{
  return is_zero(squared_radius());
}
template < class R >
inline
CircleC2<R CGAL_CTAG> CircleC2<R CGAL_CTAG>::opposite() const
{
  return CircleC2<R CGAL_CTAG>(center(),
                      squared_radius(),
                      CGAL::opposite(orientation()) );
}


template < class R >
CGAL_KERNEL_INLINE
Bbox_2 CircleC2<R CGAL_CTAG>::bbox() const
{
  double cx = CGAL::to_double(center().x());
  double cy = CGAL::to_double(center().y());
  double radius = sqrt(CGAL::to_double(squared_radius()));

  return Bbox_2(cx - radius, cy - radius, cx + radius, cy + radius);
}


template < class R >
CGAL_KERNEL_INLINE
CircleC2<R CGAL_CTAG>
CircleC2<R CGAL_CTAG>::orthogonal_transform(const CircleC2<R CGAL_CTAG>::Aff_transformation_2 &t) const
{
  Vector_2 vec(FT(1), FT(0) );   // unit vector
  vec = vec.transform(t);             // transformed
  FT  sq_scale = FT( vec*vec );       // squared scaling factor

  return CircleC2<R CGAL_CTAG>(t.transform(center()),
                      sq_scale * squared_radius(),
                      t.is_even() ? orientation()
                                  : CGAL::opposite(orientation()));
}


/*
template < class R >
inline
EllipseC2<CircleC2<R CGAL_CTAG>::FT> CircleC2<R CGAL_CTAG>::transform(
                                 const Aff_transformationC2<CircleC2<R CGAL_CTAG>::FT> &t) const
{
  return CircleC2<R CGAL_CTAG>(t.transform(center()),
                      squared_radius(),
                      orientation());
}
*/


#ifndef CGAL_NO_OSTREAM_INSERT_CIRCLEC2
template < class R >
CGAL_KERNEL_INLINE
std::ostream &operator<<(std::ostream &os, const CircleC2<R CGAL_CTAG> &c)
{
    switch(os.iword(IO::mode)) {
    case IO::ASCII :
        os << c.center() << ' ' << c.squared_radius() << ' '
           << (int)c.orientation();
        break;
    case IO::BINARY :
        os << c.center();
        write(os, c.squared_radius());
        write(os, (int)c.orientation());
        break;
    default:
        os << "CircleC2(" << c.center() <<  ", " << c.squared_radius() ;
        switch (c.orientation()) {
        case CLOCKWISE:
            os << ", clockwise)";
            break;
        case COUNTERCLOCKWISE:
            os << ", counterclockwise)";
            break;
        default:
            os << ", collinear)";
            break;
        }
        break;
    }
    return os;
}
#endif // CGAL_NO_OSTREAM_INSERT_CIRCLEC2

#ifndef CGAL_NO_ISTREAM_EXTRACT_CIRCLEC2
template < class R >
CGAL_KERNEL_INLINE
std::istream& operator>>(std::istream &is, CircleC2<R CGAL_CTAG> &c)
{
    CircleC2<R CGAL_CTAG>::Point_2 center;
    CircleC2<R CGAL_CTAG>::FT squared_radius;
    int o;
    switch(is.iword(IO::mode)) {
    case IO::ASCII :
        is >> center >> squared_radius >> o;
        break;
    case IO::BINARY :
        is >> center;
        read(is, squared_radius);
        is >> o;
        break;
    default:
        cerr << "" << std::endl;
        cerr << "Stream must be in ascii or binary mode" << endl;
        break;
    }
    c = CircleC2<R CGAL_CTAG>(center, squared_radius, (Orientation)o);
    return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_CIRCLEC2

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_CIRCLE_2_H

