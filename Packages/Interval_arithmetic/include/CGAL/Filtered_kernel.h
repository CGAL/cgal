// ============================================================================
//
// Copyright (c) 2001 The CGAL Consortium
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
// file          : include/CGAL/Filtered_kernel.h
// revision      : $Revision$
// revision_date : $Date$
// package       : ???
// author(s)     : Sylvain Pion
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

#ifndef CGAL_FILTERED_KERNEL_H
#define CGAL_FILTERED_KERNEL_H

// This file contains the definition of a generic kernel filter.
//
// TODO:
// - at the moment, it's restricted to IA filtering, but this should be easily
//   generalized to allow static filters...
// - at the moment, only the predicates are filtered.
//   Constructions will come later.
// - the kernel only works with traits only and as a pure traits only.

#include <CGAL/basic.h>
#include <CGAL/Filter_predicate.h>

CGAL_BEGIN_NAMESPACE

// This class is just used to encapsulate something.
template <class T>
class Blind_wrapper
{
public:
    typedef T  value_type;

    Blind_wrapper(const T& object) : obj(object) {}

    const T& get() const // operator() instead ?
    { return obj; }

private:
    T obj;
};

// Just forwards the construction of the CK kernel.
template < class T >
class Forward_construction
{
public:
    typedef Blind_wrapper<typename T::result_type>   result_type;

    Forward_construction() : CK_construction() {}

    result_type
    operator()() const
    { return CK_construction(); }

    template < class A1 >
    result_type
    operator()(const A1& a1) const
    { return CK_construction(a1.get()); }

    template < class A1, class A2 >
    result_type
    operator()(const A1& a1, const A2& a2) const
    { return CK_construction(a1.get(), a2.get()); }

    template < class A1, class A2, class A3 >
    result_type
    operator()(const A1& a1, const A2& a2, const A3& a3) const
    { return CK_construction(a1.get(), a2.get(), a3.get()); }

    template < class A1, class A2, class A3, class A4 >
    result_type
    operator()(const A1& a1, const A2& a2, const A3& a3, const A4& a4) const
    { return CK_construction(a1.get(), a2.get(), a3.get(), a4.get()); }

    template < class A1, class A2, class A3, class A4, class A5 >
    result_type
    operator()(const A1& a1, const A2& a2, const A3& a3, const A4& a4,
	       const A5& a5) const
    {
	return CK_construction(a1.get(), a2.get(), a3.get(), a4.get(),
		a5.get());
    }

    template < class A1, class A2, class A3, class A4, class A5, class A6 >
    result_type
    operator()(const A1& a1, const A2& a2, const A3& a3, const A4& a4,
	       const A5& a5, const A6& a6) const
    {
	return CK_construction(a1.get(), a2.get(), a3.get(), a4.get(),
		a5.get(), a6.get());
    }

    template < class A1, class A2, class A3, class A4, class A5, class A6,
               class A7 >
    result_type
    operator()(const A1& a1, const A2& a2, const A3& a3, const A4& a4,
	       const A5& a5, const A6& a6, const A7& a7) const
    {
	return CK_construction(a1.get(), a2.get(), a3.get(), a4.get(),
		a5.get(), a6.get(), a7.get());
    }

    template < class A1, class A2, class A3, class A4, class A5, class A6,
               class A7, class A8 >
    result_type
    operator()(const A1& a1, const A2& a2, const A3& a3, const A4& a4,
	       const A5& a5, const A6& a6, const A7& a7, const A8& a8) const
    {
	return CK_construction(a1.get(), a2.get(), a3.get(), a4.get(),
		a5.get(), a6.get(), a7.get(), a8.get());
    }

    template < class A1, class A2, class A3, class A4, class A5, class A6,
               class A7, class A8, class A9 >
    result_type
    operator()(const A1& a1, const A2& a2, const A3& a3, const A4& a4,
	       const A5& a5, const A6& a6, const A7& a7, const A8& a8,
	       const A9& a9) const
    {
	return CK_construction(a1.get(), a2.get(), a3.get(), a4.get(),
		a5.get(), a6.get(), a7.get(), a8.get(), a9.get());
    }

    template < class A1, class A2, class A3, class A4, class A5, class A6,
               class A7, class A8, class A9, class A10 >
    result_type
    operator()(const A1& a1, const A2& a2, const A3& a3, const A4& a4,
	       const A5& a5, const A6& a6, const A7& a7, const A8& a8,
	       const A9& a9, const A10& a10) const
    {
	return CK_construction(a1.get(), a2.get(), a3.get(), a4.get(),
		a5.get(), a6.get(), a7.get(), a8.get(), a9.get(), a10.get());
    }

    template < class A1, class A2, class A3, class A4, class A5, class A6,
               class A7, class A8, class A9, class A10, class A11 >
    result_type
    operator()(const A1& a1, const A2& a2, const A3& a3, const A4& a4,
	       const A5& a5, const A6& a6, const A7& a7, const A8& a8,
	       const A9& a9, const A10& a10, const A11& a11) const
    {
	return CK_construction(a1.get(), a2.get(), a3.get(), a4.get(),
		a5.get(), a6.get(), a7.get(), a8.get(), a9.get(), a10.get(),
		a11.get());
    }

    template < class A1, class A2, class A3, class A4, class A5, class A6,
               class A7, class A8, class A9, class A10, class A11, class A12 >
    result_type
    operator()(const A1& a1, const A2& a2, const A3& a3, const A4& a4,
	       const A5& a5, const A6& a6, const A7& a7, const A8& a8,
	       const A9& a9, const A10& a10, const A11& a11,
	       const A12& a12) const
    {
	return CK_construction(a1.get(), a2.get(), a3.get(), a4.get(),
		a5.get(), a6.get(), a7.get(), a8.get(), a9.get(), a10.get(),
		a11.get(), a12.get());
    }

    template < class A1, class A2, class A3, class A4, class A5, class A6,
               class A7, class A8, class A9, class A10, class A11, class A12,
	       class A13 >
    result_type
    operator()(const A1& a1, const A2& a2, const A3& a3, const A4& a4,
	       const A5& a5, const A6& a6, const A7& a7, const A8& a8,
	       const A9& a9, const A10& a10, const A11& a11,
	       const A12& a12, const A13& a13) const
    {
	return CK_construction(a1.get(), a2.get(), a3.get(), a4.get(),
		a5.get(), a6.get(), a7.get(), a8.get(), a9.get(), a10.get(),
		a11.get(), a12.get(), a13.get());
    }

private:
    T CK_construction;
};

// CK = construction kernel.
// EK = exact kernel called when needed by the filter.
// FK = filtering kernel
template <class CK,
          class EK = Simple_cartesian<MP_Float>,
	  class C2E_Converter = Default_converter<CK, EK>,
          class FK = Simple_cartesian<Interval_nt_advanced>,
	  class C2F_Converter = Default_converter<CK, FK> >
class Filtered_kernel
{
    typedef typename CK::Kernel_tag                       Kernel_tag; // ?

    // Define the object types :

#define CGAL_Blind_wrapper(X) typedef Blind_wrapper<typename CK::X> X;

    CGAL_Blind_wrapper(Point_2)
    CGAL_Blind_wrapper(Vector_2)
    CGAL_Blind_wrapper(Direction_2)
    CGAL_Blind_wrapper(Segment_2)
    CGAL_Blind_wrapper(Line_2)
    CGAL_Blind_wrapper(Ray_2)
    CGAL_Blind_wrapper(Triangle_2)
    CGAL_Blind_wrapper(Circle_2)
    CGAL_Blind_wrapper(Iso_rectangle_2)
    CGAL_Blind_wrapper(Aff_transformation_2)

    // CGAL_Blind_wrapper(Data_accessor_2) // ?
    // CGAL_Blind_wrapper(Conic_2)         // ?

    CGAL_Blind_wrapper(Point_3)
    CGAL_Blind_wrapper(Vector_3)
    CGAL_Blind_wrapper(Direction_3)
    CGAL_Blind_wrapper(Segment_3)
    CGAL_Blind_wrapper(Line_3)
    CGAL_Blind_wrapper(Plane_3)
    CGAL_Blind_wrapper(Ray_3)
    CGAL_Blind_wrapper(Triangle_3)
    CGAL_Blind_wrapper(Tetrahedron_3)
    CGAL_Blind_wrapper(Sphere_3)
    CGAL_Blind_wrapper(Iso_cuboid_3)
    CGAL_Blind_wrapper(Aff_transformation_3)

    // Define the predicate types and their accessor functions :

#define CGAL_Filter_pred(P, Pf) \
    typedef Filtered_predicate<EK::P, FK::P, C2E_Converter, C2F_Converter> P; \
    P Pf() const { return P(); }

    CGAL_Filter_pred(Orientation, orientation_2_object)
    ...

    // Define the construction types and their accessor functions :

#define CGAL_Filter_constr(C, Cf) \
    typedef Forward_construction<CK::C> C; \
    C Cf() const { return C(); }

    CGAL_Filter_constr(Construct_point_2, construct_point_2_object)
    ...

};

CGAL_END_NAMESPACE

#endif // CGAL_FILTERED_KERNEL_H
