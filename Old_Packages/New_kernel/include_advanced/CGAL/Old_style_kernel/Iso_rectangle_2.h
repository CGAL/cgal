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
// file          : include/CGAL/Old_style_kernel/Iso_rectangle_2.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 

#ifndef CGAL_OLD_STYLE_KERNEL_ISO_RECTANGLE_2_H
#define CGAL_OLD_STYLE_KERNEL_ISO_RECTANGLE_2_H

namespace CGAL {

template <class R_>
class Iso_rectangle_2< R_, Old_style_tag>
  : public R_::Iso_rectangle_2_base
{
 public:
  typedef  R_                        R;
  typedef typename R::RT             RT;
  typedef typename R::FT             FT;
  typedef typename R::Point_2_base   RPoint_2;
  typedef typename R::Iso_rectangle_2_base  RIso_rectangle_2;

  Iso_rectangle_2() : RIso_rectangle_2() {}
  Iso_rectangle_2(const RIso_rectangle_2& r) : RIso_rectangle_2(r) {}
  Iso_rectangle_2(const RPoint_2& p,
                  const RPoint_2& q)
    : RIso_rectangle_2(p,q) {}
};

} // namespace CGAL
#endif // CGAL_OLD_STYLE_KERNEL_ISO_RECTANGLE_2_H
