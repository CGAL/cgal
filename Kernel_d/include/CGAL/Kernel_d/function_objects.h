// Copyright (c) 1999  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Stefan Schirra

#ifndef CGAL_KERNEL_D_FUNCTION_OBJECTS_H
#define CGAL_KERNEL_D_FUNCTION_OBJECTS_H

// These functors come from the 2D-3D kernels.
// Since they have changed there, they now need to be copied here.

namespace CGAL {

namespace internal {

template <class ToBeConstructed>
class Construct
{
  public:
    typedef ToBeConstructed  result_type;

    ToBeConstructed
    operator()() const
    { return ToBeConstructed(); }

    template <class A1> 
    ToBeConstructed
    operator()( const A1& a1) const
    { return ToBeConstructed(a1); }

    template <class A1, class A2> 
    ToBeConstructed
    operator()( const A1& a1, const A2& a2) const
    { return ToBeConstructed(a1,a2); }

    template <class A1, class A2, class A3> 
    ToBeConstructed
    operator()( const A1& a1, const A2& a2, const A3& a3) const
    { return ToBeConstructed(a1,a2,a3); }

    template <class A1, class A2, class A3, class A4> 
    ToBeConstructed
    operator()( const A1& a1, const A2& a2, const A3& a3, const A4& a4) const
    { return ToBeConstructed(a1,a2,a3,a4); }

    template <class A1, class A2, class A3, class A4, class A5> 
    ToBeConstructed
    operator()( const A1& a1, const A2& a2, const A3& a3, const A4& a4,
	    const A5& a5) const
    { return ToBeConstructed(a1,a2,a3,a4,a5); }

    template <class A> 
    ToBeConstructed
    operator()( const A& a1, const A& a2, const A& a3,
                const A& a4, const A& a5, const A& a6 ) const
    { return ToBeConstructed(a1,a2,a3,a4,a5,a6); }

    template <class A> 
    ToBeConstructed
    operator()( const A& a1, const A& a2, const A& a3,
                const A& a4, const A& a5, const A& a6,
                const A& a7 ) const
    { return ToBeConstructed(a1,a2,a3,a4,a5,a6,a7); }

    template <class A> 
    ToBeConstructed
    operator()( const A& a1, const A& a2, const A& a3,
                const A& a4, const A& a5, const A& a6,
                const A& a7, const A& a8, const A& a9) const
    { return ToBeConstructed(a1,a2,a3,a4,a5,a6,a7,a8,a9); }

    template <class A> 
    ToBeConstructed
    operator()( const A& a1, const A& a2, const A& a3,
                const A& a4, const A& a5, const A& a6,
                const A& a7, const A& a8, const A& a9,
                const A& a10) const
    { return ToBeConstructed(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10); }

    template <class A> 
    ToBeConstructed
    operator()( const A& a1, const A& a2, const A& a3,
                const A& a4, const A& a5, const A& a6,
                const A& a7, const A& a8, const A& a9,
                const A& a10,const A& a11,const A& a12) const
    { return ToBeConstructed(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12); }

    template <class A> 
    ToBeConstructed
    operator()( const A& a1, const A& a2, const A& a3,
                const A& a4, const A& a5, const A& a6,
                const A& a7, const A& a8, const A& a9,
                const A& a10,const A& a11,const A& a12,
                const A& a13) const
    { return ToBeConstructed(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13); }

};

class Call_has_on_positive_side
{
  public:
    typedef bool           result_type;

    template <class Cls, class Arg>
    bool
    operator()( const Cls& c, const Arg& a) const
    { return c.has_on_positive_side(a); }
};

class Call_oriented_side
{
  public:
    typedef Oriented_side   result_type;

    template <class Cls, class Arg>
    Oriented_side
    operator()( const Cls& c, const Arg& a) const
    { return c.oriented_side(a); }
};

class Intersect
{
  public:
    typedef CGAL::Object   result_type;

    template <class T1, class T2>
    CGAL::Object
    operator()(const T1& t1, const T2& t2) const
    { return intersection( t1, t2); }
};

} // end namespace internal
} //namespace CGAL

#endif // CGAL_KERNEL_D_FUNCTION_OBJECTS_H
