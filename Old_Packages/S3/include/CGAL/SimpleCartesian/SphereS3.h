// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
// release       : 
// release_date  : 2000, October 15
//
// source        : webS3/S3.lw
// file          : include/CGAL/SimpleCartesian/SphereS3.h
// package       : S3 (1.7)
// maintainer    : Stefan Schirra <stschirr@mpi-sb.mpg.de>
// revision      : 1.7
// revision_date : 15 Oct 2000
// author(s)     : Stefan Schirra <Stefan.Schirra@@mpi-sb.mpg.de>
//                 based on code by
//                 Andreas Fabri and
//                 Herve Brönnimann
//
// coordinator   : MPI, Saarbrücken
// ======================================================================

#ifndef CGAL_SPHERES3_H
#define CGAL_SPHERES3_H

#include <CGAL/intersection_3.h>

namespace CGAL {

template <class FT>
class SphereS3
{
 public:

  SphereS3() {}

  SphereS3(const PointS3<FT>& p, const FT& s,
           const Orientation& o = COUNTERCLOCKWISE);

  SphereS3(const PointS3<FT>& p, const PointS3<FT>& q,
           const PointS3<FT>& r, const PointS3<FT>& u);

  SphereS3(const PointS3<FT>& p, const PointS3<FT>& q, const PointS3<FT>& r,
           const Orientation& o = COUNTERCLOCKWISE);

  SphereS3(const PointS3<FT>&  p, const PointS3<FT>&  q,
           const Orientation& o = COUNTERCLOCKWISE);

  SphereS3(const PointS3<FT>&  p,
           const Orientation& o = COUNTERCLOCKWISE);


  bool
  operator==(const SphereS3<FT>&) const;

  bool
  operator!=(const SphereS3<FT>&) const;

  PointS3<FT>  
  center() const;

  FT     
  squared_radius() const;

  Orientation 
  orientation() const;

  SphereS3<FT>   
  orthogonal_transform(const Aff_transformationS3<FT>& t) const;

  bool   
  is_degenerate() const;

  SphereS3<FT>   
  opposite() const;

  Oriented_side  
  oriented_side(const PointS3<FT>& p) const;

  bool   
  has_on_boundary(const PointS3<FT>& p) const;

  bool   
  has_on_positive_side(const PointS3<FT>& p) const;

  bool   
  has_on_negative_side(const PointS3<FT>& p) const;

  Bounded_side   
  bounded_side(const PointS3<FT>& p) const;

  bool   
  has_on_bounded_side(const PointS3<FT>& p) const;

  bool   
  has_on_unbounded_side(const PointS3<FT>& p) const;

  Bbox_3
  bbox() const;

 private:
  PointS3<FT>  center_;
  FT           squared_radius_;
  Orientation  orient;

};

template <class FT>
CGAL_KERNEL_CTOR_INLINE
SphereS3<FT>::
SphereS3(const PointS3<FT>& center,
         const FT& squared_radius,
         const Orientation& orient_)
 : center_(center), squared_radius_(squared_radius), orient(orient_)
{
  CGAL_kernel_precondition( ( squared_radius >= FT(0) ) &&
                            ( orient    != COLLINEAR) );
}

template <class FT>
CGAL_KERNEL_CTOR_INLINE
SphereS3<FT>::
SphereS3(const PointS3<FT>& center,
         const Orientation& orient_)
 : center_(center), squared_radius_(FT(0)), orient(orient_)
{
  CGAL_kernel_precondition( ( orient    != COLLINEAR) );
}

template <class FT>
CGAL_KERNEL_CTOR_MEDIUM_INLINE
SphereS3<FT>::
SphereS3(const PointS3<FT>& p, const PointS3<FT>& q,
         const Orientation& orient_)
{
  CGAL_kernel_precondition( orient_ != COLLINEAR);
  center_ = midpoint(p,q);
  squared_radius_ = squared_distance(p,center_);
  orient = orient_;
}

template <class FT>
CGAL_KERNEL_CTOR_MEDIUM_INLINE
SphereS3<FT>::
SphereS3(const PointS3<FT>& p, const PointS3<FT>& q, const PointS3<FT>& r,
         const Orientation& orient_)
{
  CGAL_kernel_precondition( orient_ != COLLINEAR);
  center_ = gp_linear_intersection( PlaneS3<FT>(p,q,r),
                                    bisector(p,q),
                                    bisector(p,r)); 
  squared_radius_ = squared_distance(p,center_);
  orient = orient_;
}

template <class FT>
CGAL_KERNEL_CTOR_MEDIUM_INLINE
SphereS3<FT>::
SphereS3(const PointS3<FT>& p, const PointS3<FT>& q,
         const PointS3<FT>& r, const PointS3<FT>& s)
{
  center_ = circumcenter(p,q,r,s);
  squared_radius_ = squared_distance(p,center_);
  orient = CGAL::orientation(p,q,r,s);
}

template <class FT>
CGAL_KERNEL_INLINE
bool
SphereS3<FT>::operator==(const SphereS3<FT>& t) const
{
  return center() == t.center() &&
         squared_radius() == t.squared_radius() &&
         orientation() == t.orientation();
}

template <class FT>
inline
bool
SphereS3<FT>::operator!=(const SphereS3<FT>& t) const
{ return !(*this == t); }

template <class FT>
inline
PointS3<FT>
SphereS3<FT>::center() const
{ return center_; }

template <class FT>
inline
FT
SphereS3<FT>::squared_radius() const
{ return squared_radius_; }

template <class FT>
inline
Orientation 
SphereS3<FT>::orientation() const
{ return orient; }

template <class FT>
CGAL_KERNEL_MEDIUM_INLINE
Oriented_side
SphereS3<FT>::oriented_side(const PointS3<FT>& p) const
{ return Oriented_side(bounded_side(p) * orientation()); }

template <class FT>
CGAL_KERNEL_INLINE
Bounded_side
SphereS3<FT>::bounded_side(const PointS3<FT>& p) const
{
  return Bounded_side(CGAL_NTS compare(squared_radius(),
                                       squared_distance(center(),p)));
}

template <class FT>
inline
bool
SphereS3<FT>::has_on_boundary(const PointS3<FT>& p) const
{ return squared_distance(center(),p) == squared_radius(); }

template <class FT>
CGAL_KERNEL_INLINE
bool
SphereS3<FT>::has_on_negative_side(const PointS3<FT>& p) const
{
  if (orientation() == COUNTERCLOCKWISE)
  { return has_on_unbounded_side(p); }
  return has_on_bounded_side(p);
}

template <class FT>
CGAL_KERNEL_INLINE
bool
SphereS3<FT>::has_on_positive_side(const PointS3<FT>& p) const
{
  if (orientation() == COUNTERCLOCKWISE)
  { return has_on_bounded_side(p); }
  return has_on_unbounded_side(p);
}

template <class FT>
inline
bool
SphereS3<FT>::has_on_bounded_side(const PointS3<FT>& p) const
{ return squared_distance(center(),p) < squared_radius(); }

template <class FT>
inline
bool 
SphereS3<FT>::has_on_unbounded_side(const PointS3<FT>& p) const
{ return squared_distance(center(),p) > squared_radius(); }

template <class FT>
inline
bool 
SphereS3<FT>::is_degenerate() const
{ return CGAL_NTS is_zero(squared_radius()); }

template <class FT>
inline
SphereS3<FT> 
SphereS3<FT>::opposite() const
{
  return SphereS3<FT>(center(), 
                      squared_radius(),
                      CGAL::opposite(orientation()) );
}

template <class FT>
CGAL_KERNEL_INLINE
Bbox_3
SphereS3<FT>::bbox() const
{
  // to be fixed !!!
  double cx = CGAL::to_double(center().x());
  double cy = CGAL::to_double(center().y());
  double cz = CGAL::to_double(center().z());
  double radius = CGAL::sqrt(CGAL::to_double(squared_radius()));

  return Bbox_3(cx - radius, cy - radius, cz - radius,
                cx + radius, cy + radius, cz + radius);
}

template <class FT>
CGAL_KERNEL_INLINE
SphereS3<FT>
SphereS3<FT>::orthogonal_transform(const Aff_transformationS3<FT>& t) const
{
  VectorS3<FT> vec(FT(1), FT(0) );   // unit vector
  vec = vec.transform(t);        // transformed
  FT  sq_scale = FT( vec*vec );  // squared scaling factor

  return SphereS3<FT>(t.transform(center()),
                      sq_scale * squared_radius(),
                      t.is_even() ? orientation()
                                  : CGAL::opposite(orientation()));
}

#ifndef CGAL_NO_OSTREAM_INSERT_SPHERES3
template <class FT>
CGAL_KERNEL_INLINE
std::ostream& 
operator<<(std::ostream& os, const SphereS3<FT>& c)
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
        os << "SphereS3(" << c.center() <<  ", " << c.squared_radius() ;
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
#endif // CGAL_NO_OSTREAM_INSERT_SPHERES3

#ifndef CGAL_NO_ISTREAM_EXTRACT_SPHERES3
template <class FT>
CGAL_KERNEL_INLINE
std::istream& 
operator>>(std::istream& is, SphereS3<FT>& c)
{
    PointS3<FT> center;
    FT squared_radius;
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
        std::cerr << "" << std::endl;
        std::cerr << "Stream must be in ascii or binary mode" << std::endl;
        break;
    }
    c = SphereS3<FT>(center, squared_radius, (Orientation)o);
    return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_SPHERES3

} // namespace CGAL

#endif // CGAL_SPHERES3_H
