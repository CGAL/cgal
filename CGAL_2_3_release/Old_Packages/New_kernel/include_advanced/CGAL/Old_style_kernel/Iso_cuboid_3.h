// ============================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
// 
// release       :
// release_date  :
// 
// file          : include/CGAL/Old_style_kernel/Iso_cuboid_3.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 

#ifndef CGAL_OLD_STYLE_KERNEL_ISO_CUBOID_3_H
#define CGAL_OLD_STYLE_KERNEL_ISO_CUBOID_3_H

namespace CGAL {

template <class R_>
class Iso_cuboid_3< R_, Old_style_tag> : public R_::Iso_cuboid_3_base
{
 public:
  typedef  R_                        R;
  typedef typename R::RT             RT;
  typedef typename R::FT             FT;
  typedef typename R::Point_3_base   RPoint_3;
  typedef typename R::Iso_cuboid_3_base  RIso_cuboid_3;

  Iso_cuboid_3() : RIso_cuboid_3() {}
  Iso_cuboid_3(const RIso_cuboid_3& r) : RIso_cuboid_3(r) {}
  Iso_cuboid_3(const RPoint_3& p, const RPoint_3& q)
    : RIso_cuboid_3(p,q) {}
};

} // namespace CGAL
#endif // CGAL_OLD_STYLE_KERNEL_ISO_CUBOID_3_H
