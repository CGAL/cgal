// Copyright (c) 2010-2013 CNRS and LIRIS' Establishments (France).
// All rights reserved.
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
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
#ifndef CGAL_LINEAR_CELL_COMPLEX_COPY_FUNCTORS_H
#define CGAL_LINEAR_CELL_COMPLEX_COPY_FUNCTORS_H 1

#include <CGAL/Dimension.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/Cartesian_converter.h>

/* Definition of functors used internally to copy linear cell complex attributes
 * (we need functors as attributes are stored in tuple, thus all the access
 *  must be done at compiling time).
 */
namespace CGAL
{
template< typename LCC1, typename LCC2>
struct Attribute_converter_lcc_vertex_attributes;

template< typename Functor, typename LCC1, typename LCC2>
struct Modify_attribute_converter_lcc_vertex_attributes;
// ****************************************************************************
namespace internal
{
// ****************************************************************************
// Functor allowing to set the value of a point if point exist, have
// same dimension. For dim>3, if type of points are the same (because no converter).
template<typename Point1, typename Point2>
struct Set_point_d_if_same
{
  static void run(const Point1&, Point2&)
  {}
};

template<typename Point1>
struct Set_point_d_if_same<Point1, Point1>
{
  static void run(const Point1& p1, Point1& p2)
  {
    p2 = p1; // Copy of Point_d having same type
  }
};

template<typename Point1, typename Point2,
         typename T1=typename Ambient_dimension<Point1>::type,
         typename T2=typename Ambient_dimension<Point2>::type>
struct Set_point_if_possible
{
  static void run(const Point1&, Point2&)
  {}
};

template<typename Point1, typename Point2>
struct Set_point_if_possible<Point1, Point2,
    Dimension_tag<2>, Dimension_tag<2> >
{
  static void run(const Point1& p1, Point2& p2)
  {
    p2 = Cartesian_converter<typename Kernel_traits<Point1>::Kernel,
        typename Kernel_traits<Point2>::Kernel>(p1);
  }
};

template<typename Point1, typename Point2>
struct Set_point_if_possible<Point1, Point2,
    Dimension_tag<3>, Dimension_tag<3> >
{
  static void run(const Point1& p1, Point2& p2)
  {
    p2 = Cartesian_converter<typename Kernel_traits<Point1>::Kernel,
        typename Kernel_traits<Point2>::Kernel>()(p1);
  }
};

template<typename Point1, typename Point2>
struct Set_point_if_possible<Point1, Point2, Dynamic_dimension_tag,
                             Dynamic_dimension_tag>
{
  static void run(const Point1& p1, Point2& p2)
  {
    if ( p1.dimension()==p2.dimension() )
      Set_point_d_if_same<Point1, Point2>::run(p1, p2);
  }
};

// Set_point_if_exist if Attr1 has a point
template< typename Attr1, typename Attr2,
          typename Point1=typename Attr1::Point,
          typename Point2=typename Attr2::Point >
struct Set_point_if_exist
{
  static void run(const Attr1& a1, Attr2& a2)
  {
    Set_point_if_possible<Point1, Point2>::run(a1.point(), a2.point());
  }
};

template<typename Attr1, typename Attr2, typename Point1 >
struct Set_point_if_exist<Attr1, Attr2, Point1, CGAL::Void>
{
  static void run(const Attr1&, Attr2&)
  {}
};
// ****************************************************************************
#ifndef CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES
template <typename LCC1, typename LCC2, class Tuple,
          class Res=CGAL::cpp11::tuple<> >
struct Modify_tuple_of_converter_for_vertex_attribute;

// empty tuple and empty res: res is a tuple with only
// Point_converter_lcc_vertex_attributes
template <typename LCC1, typename LCC2>
struct Modify_tuple_of_converter_for_vertex_attribute
    <LCC1, LCC2, CGAL::cpp11::tuple<>, CGAL::cpp11::tuple<> >
{
  typedef CGAL::cpp11::tuple
  <Point_converter_lcc_vertex_attributes<LCC1, LCC2> > type;

  static type run(const CGAL::cpp11::tuple<>&)
  { return type(); }
};

// empty tuple but non empty res: we already have modified vertex attribute
template < typename LCC1, typename LCC2, class ... Res >
struct Modify_tuple_of_converter_for_vertex_attribute
    <LCC1, LCC2, CGAL::cpp11::tuple<>, CGAL::cpp11::tuple<Res...> >
{
  typedef CGAL::cpp11::tuple<Res...> type;

  static type run(const CGAL::cpp11::tuple<>&)
  { CGAL_assertion(false); return type(); }
};

// empty res but non empty empty tuple, Firsttype is vertex attribute.
template < typename LCC1, typename LCC2, class Firsttype, class ... Tuple >
struct Modify_tuple_of_converter_for_vertex_attribute
    <LCC1, LCC2, CGAL::cpp11::tuple<Firstype, Tuple...>, CGAL::cpp11::tuple<> >
{
  typedef typename Modify_tuple_of_converter_for_vertex_attribute
  <LCC1, LCC2, CGAL::cpp11::tuple<Tuple...>,
  CGAL::cpp11::tuple
  <Modify_attribute_converter_lcc_vertex_attributes<Firsttype, LCC1, LCC2> >::
  type type;

  static type run(const CGAL::cpp11::tuple<Firsttype, Tuple...>& t)
  {
    return type
        (Modify_attribute_converter_lcc_vertex_attributes<Firsttype, LCC1, LCC2>
         (t.get_head()),
         t.get_tail());
  }
};

// non empty res, non empty tuple, we copy
template < typename LCC1, typename LCC2, class Firsttype,
           class ... Tuple, class ... Res >
struct Modify_tuple_of_converter_for_vertex_attribute
    <LCC1, LCC2, CGAL::cpp11::tuple<Firsttype,Tuple...>,
    CGAL::cpp11::tuple<Res...> >
{
  typedef typename Modify_tuple_of_converter_for_vertex_attribute
  <LCC1, LCC2, CGAL::cpp11::tuple<Tuple...>,
  CGAL::cpp11::tuple<Res...,Functor> >::type type;

  static type run(const CGAL::cpp11::tuple<>&)
  { CGAL_assertion(false); return type; }
};
#else // CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES
template <class LCC1, class LCC2, class Tuple>
struct Modify_tuple_of_converter_for_vertex_attribute;

template <class LCC1, class LCC2>
struct Modify_tuple_of_converter_for_vertex_attribute
    <LCC1,LCC2,CGAL::cpp11::tuple<> > {
  typedef CGAL::cpp11::tuple
  <Attribute_converter_lcc_vertex_attributes<LCC1, LCC2> > type;

  static type run(const CGAL::cpp11::tuple<>&)
  { return type(); }
};

template <class LCC1, class LCC2, class T1>
struct Modify_tuple_of_converter_for_vertex_attribute
    <LCC1,LCC2,CGAL::cpp11::tuple<T1> > {
  typedef CGAL::cpp11::tuple
  <Modify_attribute_converter_lcc_vertex_attributes<T1, LCC1, LCC2> > type;

  static type run(const CGAL::cpp11::tuple<T1>& t)
  {
    return type
        (Modify_attribute_converter_lcc_vertex_attributes<T1, LCC1, LCC2>
         (t.get_head()),
         t.get_tail());
  }
};

template <class LCC1, class LCC2, class T1, class T2>
struct Modify_tuple_of_converter_for_vertex_attribute
    <LCC1,LCC2,CGAL::cpp11::tuple<T1,T2> > {
  typedef CGAL::cpp11::tuple
  <Modify_attribute_converter_lcc_vertex_attributes<T1, LCC1, LCC2>,T2> type;

  static type run(const CGAL::cpp11::tuple<T1,T2>& t)
  {
    return type
        (Modify_attribute_converter_lcc_vertex_attributes<T1, LCC1, LCC2>
         (t.get_head()),
         t.get_tail());
  }
};

template <class LCC1, class LCC2, class T1, class T2, class T3>
struct Modify_tuple_of_converter_for_vertex_attribute
    <LCC1,LCC2,CGAL::cpp11::tuple<T1,T2,T3> > {
  typedef CGAL::cpp11::tuple
  <Modify_attribute_converter_lcc_vertex_attributes<T1, LCC1, LCC2>,T2,T3>
  type;

  static type run(const CGAL::cpp11::tuple<T1,T2,T3>& t)
  {
    return type
        (Modify_attribute_converter_lcc_vertex_attributes<T1, LCC1, LCC2>
         (t.get_head()),
         t.get_tail());
  }
};

template <class LCC1, class LCC2, class T1, class T2, class T3, class T4>
struct Modify_tuple_of_converter_for_vertex_attribute
    <LCC1,LCC2,CGAL::cpp11::tuple<T1,T2,T3,T4> > {
  typedef CGAL::cpp11::tuple
  <Modify_attribute_converter_lcc_vertex_attributes<T1, LCC1, LCC2>,T2,T3,T4>
  type;

  static type run(const CGAL::cpp11::tuple<T1,T2,T3,T4>& t)
  {
    return type
        (Modify_attribute_converter_lcc_vertex_attributes<T1, LCC1, LCC2>
         (t.get_head()),
         t.get_tail());
  }
};

template <class LCC1, class LCC2, class T1, class T2, class T3, class T4,
          class T5>
struct Modify_tuple_of_converter_for_vertex_attribute
    <LCC1,LCC2,CGAL::cpp11::tuple<T1,T2,T3,T4,T5> > {
  typedef CGAL::cpp11::tuple
  <Modify_attribute_converter_lcc_vertex_attributes<T1, LCC1, LCC2>,T2,T3,T4,
  T5> type;

  static type run(const CGAL::cpp11::tuple<T1,T2,T3,T4,T5>& t)
  {
    return type
        (Modify_attribute_converter_lcc_vertex_attributes<T1, LCC1, LCC2>
         (t.get_head()),
         t.get_tail());
  }
};

template <class LCC1, class LCC2, class T1, class T2, class T3, class T4,
          class T5, class T6>
struct Modify_tuple_of_converter_for_vertex_attribute
    <LCC1,LCC2,CGAL::cpp11::tuple<T1,T2,T3,T4,T5,T6> > {
  typedef CGAL::cpp11::tuple
  <Modify_attribute_converter_lcc_vertex_attributes<T1, LCC1, LCC2>,T2,T3,T4,
  T5,T6> type;

  static type run(const CGAL::cpp11::tuple<T1,T2,T3,T4,T5,T6>& t)
  {
    return type
        (Modify_attribute_converter_lcc_vertex_attributes<T1, LCC1, LCC2>
         (t.get_head()),
         t.get_tail());
  }
};

template <class LCC1, class LCC2, class T1, class T2, class T3, class T4,
          class T5, class T6, class T7>
struct Modify_tuple_of_converter_for_vertex_attribute
    <LCC1,LCC2,CGAL::cpp11::tuple<T1,T2,T3,T4,T5,T6,T7> > {
  typedef CGAL::cpp11::tuple
  <Modify_attribute_converter_lcc_vertex_attributes<T1, LCC1, LCC2>,T2,T3,T4,
  T5,T6,T7> type;

  static type run(const CGAL::cpp11::tuple<T1,T2,T3,T4,T5,T6,T7>& t)
  {
    return type
        (Modify_attribute_converter_lcc_vertex_attributes<T1, LCC1, LCC2>
         (t.get_head()),
         t.get_tail());
  }
};

template <class LCC1, class LCC2, class T1, class T2, class T3, class T4,
          class T5, class T6, class T7, class T8>
struct Modify_tuple_of_converter_for_vertex_attribute
    <LCC1,LCC2,CGAL::cpp11::tuple<T1,T2,T3,T4,T5,T6,T7,T8> > {
  typedef CGAL::cpp11::tuple
  <Modify_attribute_converter_lcc_vertex_attributes<T1, LCC1, LCC2>,T2,T3,T4,
  T5,T6,T7,T8> type;

  static type run(const CGAL::cpp11::tuple<T1,T2,T3,T4,T5,T6,T7,T8>& t)
  {
    return type
        (Modify_attribute_converter_lcc_vertex_attributes<T1, LCC1, LCC2>
         (t.get_head()),
         t.get_tail());
  }
};

template <class LCC1, class LCC2, class T1, class T2, class T3, class T4,
          class T5, class T6, class T7, class T8, class T9>
struct Modify_tuple_of_converter_for_vertex_attribute
    <LCC1,LCC2,CGAL::cpp11::tuple<T1,T2,T3,T4,T5,T6,T7,T8,T9> > {
  typedef CGAL::cpp11::tuple
  <Modify_attribute_converter_lcc_vertex_attributes<T1, LCC1, LCC2>,T2,T3,T4,
  T5,T6,T7,T8,T9> type;

  static type run(const CGAL::cpp11::tuple<T1,T2,T3,T4,T5,T6,T7,T8,T9>& t)
  {
    return type
        (Modify_attribute_converter_lcc_vertex_attributes<T1, LCC1, LCC2>
         (t.get_head()),
         t.get_tail());
  }
};
#endif
// ****************************************************************************
} // namespace internal
// LCC1 is the existing map, to convert into map2.
// ****************************************************************************
template< typename LCC1, typename LCC2>
struct Attribute_converter_lcc_vertex_attributes
{
  typename LCC2::Vertex_attribute_handle operator()
  (const LCC1& map1, LCC2& map2, typename LCC1::Dart_const_handle dh1,
   typename LCC2::Dart_handle dh2) const
  {
    map2.set_vertex_attribute(dh2, map2.create_vertex_attribute());
    internal::Set_point_if_exist<typename LCC1::Vertex_attribute,
                       typename LCC2::Vertex_attribute>::
        run( dh1->vertex_attribute(), dh2->vertex_attribute() );
    return dh2->vertex_attribute();
  }
};
// ****************************************************************************
template< typename Functor, typename LCC1, typename LCC2>
struct Modify_attribute_converter_lcc_vertex_attributes
{
  Modify_attribute_converter_lcc_vertex_attributes(const Functor& f): myf(f)
  {}

  typename LCC2::Vertex_attribute_handle operator()
  (const LCC1& map1, LCC2& map2, typename LCC1::Dart_const_handle dh1,
   typename LCC2::Dart_handle dh2) const
  {
    myf(map1, map2, dh1, dh2);
    if ( dh2->vertex_attribute()==NULL )
      map2.set_vertex_attribute(dh2, map2.create_vertex_attribute());
    internal::Set_point_if_exist<typename LCC1::Vertex_attribute,
                       typename LCC2::Vertex_attribute>::
        run( dh1->vertex_attribute(), dh2->vertex_attribute() );
    return dh2->vertex_attribute();
  }
  private:
  Functor myf;
};
// Case where the two i-attributes are non void.
template< typename LCC1, typename LCC2, unsigned int i,
          typename Info1, typename Info2,
          typename Point2 >
struct Default_converter_two_non_void_attributes_lcc
{ // Here Info1!=Info2 but Point2!=CGAL::Void (thus Linear_cell_complex)
  static typename LCC2::template Attribute_handle<i>::type
  run(LCC2& map2, typename LCC1::template Attribute_const_handle<i>::type ah)
  {
    typename LCC2::template Attribute_handle<i>::type
        res=map2.template create_attribute<i>();
    if ( ah!=NULL )
    {
      // Copy the point of ah if it exists and have same dimension
      CGAL::internal::Set_point_if_exist<typename LCC1::template Attribute_type<i>::type,
          typename LCC2::template Attribute_type<i>::type>::run( *ah, *res );
    }
    return res;
  }
};

template< typename LCC1, typename LCC2, unsigned int i, typename Info, typename Point2 >
struct Default_converter_two_non_void_attributes_lcc<LCC1, LCC2, i, Info, Info, Point2>
{ // Here Info1==Info2 but Point2!=CGAL::Void (thus Linear_cell_complex)
  static typename LCC2::template Attribute_handle<i>::type
  run(LCC2& map2, typename LCC1::template Attribute_const_handle<i>::type ah)
  {
    typename LCC2::template Attribute_handle<i>::type
        res=map2.template create_attribute<i>();
    if ( ah!=NULL )
    {
      res->info()=ah->info();
      // Copy the point of ah if it exists and have same dimension
      CGAL::internal::Set_point_if_exist<typename LCC1::template Attribute_type<i>::type,
          typename LCC2::template Attribute_type<i>::type>::run( *ah, *res );
    }
    return res;
  }
};

template< typename LCC1, typename LCC2, unsigned int i, typename Point2 >
struct Default_converter_two_non_void_attributes_lcc<LCC1, LCC2, i, void, void, Point2>
{ // Here Info1==Info2==void but Point2!=CGAL::Void (thus Linear_cell_complex)
  static typename LCC2::template Attribute_handle<i>::type
  run(LCC2& map2, typename LCC1::template Attribute_const_handle<i>::type ah)
  {
    typename LCC2::template Attribute_handle<i>::type res=
      map2.template create_attribute<i>();
    if ( ah!=NULL )
    {
      // Copy the point of ah if it exists and have same dimension
      CGAL::internal::Set_point_if_exist<typename LCC1::template Attribute_type<i>::type,
          typename LCC2::template Attribute_type<i>::type>::
          run( *ah, *res );
    }
    return res;
  }
};
// ****************************************************************************
} // namespace CGAL

#endif // CGAL_LINEAR_CELL_COMPLEX_COPY_FUNCTORS_H
