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
// release_date  : 
//
// file          : include/CGAL/Kernel/traits_aids.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
// coordinator   : MPI, Saarbruecken
// ======================================================================

#ifndef CGAL_KERNEL_TRAITS_AIDS_H
#define CGAL_KERNEL_TRAITS_AIDS_H

namespace CGAL {

template <class LeftTurn>
class Right_turn_by_left_turn
{
  public:
    Right_turn_by_left_turn(const LeftTurn& l) : left_turn(l) {}

    template <class Point>
    bool
    operator()(const Point& p, const Point& q, const Point& r) const
    { return  left_turn( q, p, r); }
    
  private:
    LeftTurn   left_turn;
};

#ifndef CGAL_NO_DEPRECATED_CODE
template <class Leftturn>
class Rightturn_by_leftturn
{
  public:
    Rightturn_by_leftturn(const Leftturn& l) : leftturn(l) {}

    template <class Point>
    bool
    operator()(const Point& p, const Point& q, const Point& r) const
    { return  left_turn( q, p, r); }
    
  private:
    Leftturn   leftturn;
};
#endif

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
