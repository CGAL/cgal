// ==========================================================================
//
// Copyright (c) 1998 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// --------------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/Cartesian/Scaling_rep_3.h
// source        : include/CGAL/Cartesian/Scaling_rep_3.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve.Bronnimann@sophia.inria.fr
//
// coordinator   : INRIA Sophia-Antipolis (Herve.Bronnimann@sophia.inria.fr)
//
// ==========================================================================

#ifndef CGAL_CARTESIAN_REDEFINE_NAMES_3_H
#define CGAL_CTAG
#endif

#ifdef _MSC_VER
#define typename
#endif

#include <CGAL/Cartesian/distance_computations_3.h>

#ifndef CGAL_CARTESIAN_SEGMENT_3_C
#define CGAL_CARTESIAN_SEGMENT_3_C

CGAL_BEGIN_NAMESPACE

template < class R >
inline
_Twotuple< typename SegmentC3<R CGAL_CTAG>::Point_3 > *
SegmentC3<R CGAL_CTAG>::ptr() const
{
  return (_Twotuple< Point_3 >*)PTR;
}

t#ifndef CGAL_CARTESIAN_SCALING_REP_3_H
#define CGAL_CARTESIAN_SCALING_REP_3_H

#ifndef CGAL_CARTESIAN_REDEFINE_NAMES_3_H
#include <CGAL/Cartesian/redefine_names_3.C>
#endif

CGAL_BEGIN_NAMESPACE

template < class R >
class Scaling_repC3 : public Aff_transformation_rep_baseC3<R>
{
public:
  typedef typename R::FT                        FT;
  typedef typename R::RT                        RT;
  typedef Aff_transformation_rep_baseC3<R>      Aff_t_base_3;
  typedef typename Aff_t_base_3::Point_3                Point_3;
  typedef typename Aff_t_base_3::Vector_3               Vector_3;
  typedef typename Aff_t_base_3::Direction_3            Direction_3;
  typedef typename Aff_t_base_3::Plane_3                Plane_3;
  typedef typename Aff_t_base_3::Aff_transformation_3   Aff_transformation_3;

friend Aff_transformation_3 operator* CGAL_NULL_TMPL_ARGS
                              (const Aff_transformation_3 &a,
                               const Aff_transformation_3 &b);

  Scaling_repC3() {}
  Scaling_repC3(const FT &s) : scalefactor(s) {}
  ~Scaling_repC3() {}

  Point_3      transform(const Point_3 &p) const
  {
    return Point_3(scalefactor * p.x(), scalefactor * p.y(),
                   scalefactor * p.z());
  }

  Vector_3     transform(const Vector_3 &v) const
  {
    return Vector_3(scalefactor * v.x(), scalefactor * v.y(),
                    scalefactor * v.z());
  }

  Direction_3  transform(const Direction_3 &d) const
  {
    return d;
  }

  Aff_transformation_3 inverse() const
  {
    return Aff_transformation_3(SCALING, FT(1)/scalefactor);
  }

  Aff_transformation_3 general_form() const
  {
    FT ft0(0);

    return Aff_transformation_3(scalefactor, ft0, ft0,
                                ft0, scalefactor, ft0,
                                ft0, ft0, scalefactor);
  }

  Aff_transformation_3 transpose() const
  {
    return general_form();
  }

  virtual bool            is_even() const
  {
    return true;
  }

  virtual FT cartesian(int i, int j) const
  {
    if (i!=j) return FT(0);
    if (i==3) return FT(1);
    return scalefactor;
  }


  virtual std::ostream &print(std::ostream &os) const
  {
    FT ft0(0);
    os << "                   "<< scalefactor <<' '<< ft0 <<' '<< ft0 << ' '
       << ft0 << "\n";
    os << "                   "<< ft0 <<' '<< scalefactor <<' '<< ft0 << ' '
       << ft0 << "\n";
    os << "                   "<< ft0 <<' '<<  ft0 <<' '<< scalefactor <<' '
       << ft0 << ")\n";

    return os;
  }

private:
  FT   scalefactor;
};

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_SCALING_REP_3_H
