// ======================================================================
//
// Copyright (c) 2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
// release       : 
// release_date  : 2000, December 10
//
// source        : Traits_helpers.lw
// file          : include/CGAL/Kernel/traits_aids.h
// package       : Kernel_basic (3.17)
// maintainer    : Stefan Schirra <stschirr@mpi-sb.mpg.de>
// revision      : 0.95
// revision_date : 22 Feb 2000
// author(s)     : Stefan Schirra
//
// coordinator   : MPI, Saarbruecken
// ======================================================================

#ifndef CGAL_KERNEL_TRAITS_AIDS_H
#define CGAL_KERNEL_TRAITS_AIDS_H

namespace CGAL {

template <class Leftturn>
class Rightturn_by_leftturn
{
  public:
    Rightturn_by_leftturn(const Leftturn& l) : leftturn(l) {}

    template <class Point>
    bool
    operator()(const Point& p, const Point& q, const Point& r) const
    { return  leftturn( q, p, r); }
    
  private:
    Leftturn   leftturn;
};
template <class BinaryPredicate>
class Invert_binary
{
  public:
    Invert_binary(const BinaryPredicate& p) : pred(p) {}

    template <class Arg1, class Arg2>
    bool
    operator()(const Arg1& a1, const Arg2& a2)
    { return !pred(a1,a2); }

  protected:
    BinaryPredicate  pred;
};

} // namespace CGAL

#endif // CGAL_KERNEL_TRAITS_AIDS_H
