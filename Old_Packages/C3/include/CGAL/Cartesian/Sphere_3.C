#ifndef CGAL_CARTESIAN_REDEFINE_NAMES_3_H
#define CGAL_CTAG
#endif

#ifndef CGAL_CARTESIAN_SPHERE_3_C
#define CGAL_CARTESIAN_SPHERE_3_C

CGAL_BEGIN_NAMESPACE

template < class R >
inline
Sphere_repC3<R> *SphereC3<R CGAL_CTAG>::ptr() const
{
  return (Sphere_repC3<R>*)PTR;
}

template < class R >
CGAL_KERNEL_CTOR_INLINE
SphereC3<R CGAL_CTAG>::SphereC3()
{
  PTR = new Sphere_repC3<R> ;
}

template < class R >
CGAL_KERNEL_CTOR_INLINE
SphereC3<R CGAL_CTAG>::
SphereC3(const SphereC3<R CGAL_CTAG> &t)
  : Handle((Handle&)t)
{}

template < class R >
CGAL_KERNEL_CTOR_INLINE
SphereC3<R CGAL_CTAG>::
SphereC3(const SphereC3<R CGAL_CTAG>::Point_3 &center,
                      const typename R::FT &squared_radius,
                      const Orientation &orient = COUNTERCLOCKWISE)
{
  CGAL_kernel_precondition( ( squared_radius >= FT(0) ) &&
                            ( orient    != COLLINEAR) );

  PTR = new Sphere_repC3<R>(center, squared_radius, orient);
}

template < class R >
CGAL_KERNEL_CTOR_INLINE
SphereC3<R CGAL_CTAG>::
SphereC3(const SphereC3<R CGAL_CTAG>::Point_3 &center,
         const Orientation &orient = COUNTERCLOCKWISE)
{
  CGAL_kernel_precondition( ( orient    != COLLINEAR) );

  PTR = new Sphere_repC3<R>(center, FT(0), orient);
}

template < class R >
CGAL_KERNEL_CTOR_MEDIUM_INLINE
SphereC3<R CGAL_CTAG>::
SphereC3(const SphereC3<R CGAL_CTAG>::Point_3 &p,
         const SphereC3<R CGAL_CTAG>::Point_3 &q,
         const Orientation &orient = COUNTERCLOCKWISE)
{
  CGAL_kernel_precondition( orient != COLLINEAR);

  SphereC3<R CGAL_CTAG>::Point_3 center = midpoint(p,q);
  SphereC3<R CGAL_CTAG>::FT      squared_radius = squared_distance(p,center);

  PTR = new Sphere_repC3<R>( center, squared_radius, orient);
}

template < class R >
CGAL_KERNEL_CTOR_MEDIUM_INLINE
SphereC3<R CGAL_CTAG>::
SphereC3(const SphereC3<R CGAL_CTAG>::Point_3 &p,
         const SphereC3<R CGAL_CTAG>::Point_3 &q,
         const SphereC3<R CGAL_CTAG>::Point_3 &r,
         const Orientation &orient = COUNTERCLOCKWISE)
{
  CGAL_kernel_precondition( orient != COLLINEAR);

  /****** ADD CIRCUMCENTER OF 3 POINTS IN 3D ******/
  Point_3 center = circumcenter(p,q,r);
  FT      squared_radius = squared_distance(p,center);

  PTR = new Sphere_repC3<R>(center, squared_radius, orient);
}

template < class R >
CGAL_KERNEL_CTOR_MEDIUM_INLINE
SphereC3<R CGAL_CTAG>::
SphereC3(const SphereC3<R CGAL_CTAG>::Point_3 &p,
         const SphereC3<R CGAL_CTAG>::Point_3 &q,
         const SphereC3<R CGAL_CTAG>::Point_3 &r,
         const SphereC3<R CGAL_CTAG>::Point_3 &s)
{
  Point_3 center = circumcenter(p,q,r,s);
  FT      squared_radius = squared_distance(p,center);

  PTR = new Sphere_repC3<R>(center, squared_radius, orient);
}

template < class R >
inline
SphereC3<R CGAL_CTAG>::~SphereC3()
{}

template < class R >
inline
SphereC3<R CGAL_CTAG> &SphereC3<R CGAL_CTAG>::
operator=(const SphereC3<R CGAL_CTAG> &t)
{
  Handle::operator=(t);
  return *this;
}

template < class R >
CGAL_KERNEL_INLINE
bool SphereC3<R CGAL_CTAG>::
operator==(const SphereC3<R CGAL_CTAG> &t) const
{
   return (center() == t.center()) &&
          (squared_radius() == t.squared_radius() &&
          orientation() == t.orientation());
}

template < class R >
inline
bool SphereC3<R CGAL_CTAG>::
operator!=(const SphereC3<R CGAL_CTAG> &t) const
{
  return !(*this == t);
}

template < class R >
int SphereC3<R CGAL_CTAG>::id() const
{
  return (int)PTR;
}

template < class R >
inline
SphereC3<R CGAL_CTAG>::Point_3
SphereC3<R CGAL_CTAG>::center() const
{
 return ptr()->center;
}

template < class R >
inline
SphereC3<R CGAL_CTAG>::FT
SphereC3<R CGAL_CTAG>::squared_radius() const
{
 return ptr()->squared_radius;
}

template < class R >
inline
Orientation SphereC3<R CGAL_CTAG>::orientation() const
{
 return ptr()->orient;
}

template < class R >
CGAL_KERNEL_MEDIUM_INLINE
Oriented_side
SphereC3<R CGAL_CTAG>::
oriented_side(const SphereC3<R CGAL_CTAG>::Point_3 &p) const
{
  return Oriented_side(bounded_side(p) * orientation());
}

template < class R >
CGAL_KERNEL_INLINE
Bounded_side
SphereC3<R CGAL_CTAG>::
bounded_side(const SphereC3<R CGAL_CTAG>::Point_3 &p) const
{
  return Bounded_side(CGAL::compare(squared_radius(),
                                    squared_distance(center(),p)));
}

template < class R >
inline
bool
SphereC3<R CGAL_CTAG>::
has_on_boundary(const SphereC3<R CGAL_CTAG>::Point_3 &p) const
{
  return squared_distance(center(),p) == squared_radius();
}

template < class R >
CGAL_KERNEL_INLINE
bool
SphereC3<R CGAL_CTAG>::
has_on_negative_side(const SphereC3<R CGAL_CTAG>::Point_3 &p) const
{
  if (orientation() == COUNTERCLOCKWISE) {
    return has_on_unbounded_side(p);
  }
  return has_on_bounded_side(p);
}

template < class R >
CGAL_KERNEL_INLINE
bool
SphereC3<R CGAL_CTAG>::
has_on_positive_side(const SphereC3<R CGAL_CTAG>::Point_3 &p) const
{
  if (orientation() == COUNTERCLOCKWISE) {
    return has_on_bounded_side(p);
  }
  return has_on_unbounded_side(p);
}

template < class R >
inline
bool
SphereC3<R CGAL_CTAG>::
has_on_bounded_side(const SphereC3<R CGAL_CTAG>::Point_3 &p) const
{
  return squared_distance(center(),p) < squared_radius();
}

template < class R >
inline
bool SphereC3<R CGAL_CTAG>::
has_on_unbounded_side(const SphereC3<R CGAL_CTAG>::Point_3 &p) const
{
  return squared_distance(center(),p) > squared_radius();
}

template < class R >
inline
bool SphereC3<R CGAL_CTAG>::
is_degenerate() const
{
  return is_zero(squared_radius());
}

template < class R >
inline
SphereC3<R CGAL_CTAG> SphereC3<R CGAL_CTAG>::
opposite() const
{
  return SphereC3<R CGAL_CTAG>(center(),
                      squared_radius(),
                      CGAL::opposite(orientation()) );
}

template < class R >
CGAL_KERNEL_INLINE
Bbox_3
SphereC3<R CGAL_CTAG>::bbox() const
{
  double cx = CGAL::to_double(center().x());
  double cy = CGAL::to_double(center().y());
  double cz = CGAL::to_double(center().z());
  double radius = sqrt(CGAL::to_double(squared_radius()));

  return Bbox_3(cx - radius, cy - radius, cz - radius,
                cx + radius, cy + radius, cz + radius);
}

template < class R >
CGAL_KERNEL_INLINE
SphereC3<R CGAL_CTAG>
SphereC3<R CGAL_CTAG>::
orthogonal_transform(const SphereC3<R CGAL_CTAG>::Aff_transformation_3 &t) const
{
  Vector_3 vec(FT(1), FT(0) );   // unit vector
  vec = vec.transform(t);             // transformed
  FT  sq_scale = FT( vec*vec );       // squared scaling factor

  return SphereC3<R CGAL_CTAG>(t.transform(center()),
                      sq_scale * squared_radius(),
                      t.is_even() ? orientation()
                                  : CGAL::opposite(orientation()));
}

/*
template < class R >
inline
EllipseC3<SphereC3<R CGAL_CTAG>::FT> SphereC3<R CGAL_CTAG>::transform(
                                 const Aff_transformationC3<SphereC3<R CGAL_CTAG>::FT> &t) const
{
  return SphereC3<R CGAL_CTAG>(t.transform(center()),
                      squared_radius(),
                      orientation());
}
*/


#ifndef CGAL_NO_OSTREAM_INSERT_SPHEREC3
template < class R >
CGAL_KERNEL_INLINE
std::ostream &operator<<(std::ostream &os, const SphereC3<R CGAL_CTAG> &c)
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
        os << "SphereC3(" << c.center() <<  ", " << c.squared_radius() ;
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
#endif // CGAL_NO_OSTREAM_INSERT_SPHEREC3

#ifndef CGAL_NO_ISTREAM_EXTRACT_SPHEREC3
template < class R >
CGAL_KERNEL_INLINE
std::istream& operator>>(std::istream &is, SphereC3<R CGAL_CTAG> &c)
{
    SphereC3<R CGAL_CTAG>::Point_3 center;
    SphereC3<R CGAL_CTAG>::FT squared_radius;
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
    c = SphereC3<R CGAL_CTAG>(center, squared_radius, (Orientation)o);
    return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_SPHEREC3

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_SPHERE_3_H

