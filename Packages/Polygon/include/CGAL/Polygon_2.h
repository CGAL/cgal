// ============================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-0.9-I-06 $
// release_date  : $CGAL_Date: 1998/03/11 $
//
// file          : include/CGAL/Polygon_2.h
// source        : 
// revision      : 1.8a
// revision_date : 13 Mar 1998
// author(s)     : Wieger Wesselink <wieger@cs.ruu.nl>
//
// coordinator   : Utrecht University
//
// ============================================================================

#ifndef CGAL_POLYGON_2_H
#define CGAL_POLYGON_2_H

#ifndef CGAL_BASIC_H
#include <CGAL/basic.h>
#endif // CGAL_BASIC_H

#ifdef CGAL_CFG_NO_MEMBER_TEMPLATES
#include <vector.h>
#include <list.h>
#endif

#include <iterator.h>

#ifndef CGAL_CIRCULATOR_H
#include <CGAL/circulator.h>
#endif // CGAL_CIRCULATOR_H
#ifndef CGAL_ENUM_H
#include <CGAL/enum.h>
#endif // CGAL_ENUM_H

#ifdef CGAL_REP_CLASS_DEFINED
#ifndef CGAL_POLYGON_TRAITS_2_H
#include <CGAL/Polygon_traits_2.h>
#endif // CGAL_POLYGON_TRAITS_2_H
#ifndef CGAL_AFF_TRANSFORMATION_2_H
#include <CGAL/Aff_transformation_2.h>
#endif // CGAL_AFF_TRANSFORMATION_2_H
#endif // CGAL_REP_CLASS_DEFINED

#ifndef CGAL_POLYGON_2_ALGORITHMS_H
#include <CGAL/Polygon_2_algorithms.h>
#endif // CGAL_POLYGON_2_ALGORITHMS_H
#ifndef CGAL_POLYGON_2_EDGE_ITERATOR_H
#include <CGAL/Polygon_2_edge_iterator.h>
#endif // CGAL_POLYGON_2_EDGE_ITERATOR_H
#ifndef CGAL_POLYGON_2_EDGE_CIRCULATOR_H
#include <CGAL/Polygon_2_edge_circulator.h>
#endif // CGAL_POLYGON_2_EDGE_CIRCULATOR_H

//-----------------------------------------------------------------------//
//                          CGAL_Polygon_2
//-----------------------------------------------------------------------//

template <class _Traits, class _Container>
class CGAL_Polygon_2 {
  private:
    _Container d_container;

  public:
    //--------------------------------------------------------
    //             Types
    //--------------------------------------------------------

    typedef _Traits Traits;
    typedef _Container Container;

    typedef typename _Traits::FT FT;
    typedef typename _Traits::Point_2 Point_2;
    typedef typename _Traits::Segment_2 Segment_2;

    typedef typename _Container::difference_type difference_type;
    typedef typename _Container::value_type value_type;

    //-------------------------------------------------------//
    // this intermediary step is required by Sun C++ 4.1
    typedef typename _Container::iterator iterator;
    typedef typename _Container::const_iterator const_iterator;
    //-------------------------------------------------------//

    typedef iterator Vertex_iterator;
    typedef const_iterator Vertex_const_iterator;

    typedef CGAL_Bidirectional_circulator_from_container<_Container>
            Vertex_circulator;

    typedef CGAL_Bidirectional_const_circulator_from_container<_Container>
            Vertex_const_circulator;

    typedef CGAL_Polygon_2_edge_iterator<_Traits,_Container>
            Edge_const_iterator;

    typedef CGAL_Polygon_2_const_edge_circulator<_Traits,_Container>
            Edge_const_circulator;

    //--------------------------------------------------------
    //             Creation
    //--------------------------------------------------------

    CGAL_Polygon_2()
      { }

    CGAL_Polygon_2(const CGAL_Polygon_2<_Traits,_Container>& polygon)
      : d_container(polygon.d_container) { }

    CGAL_Polygon_2<_Traits,_Container>&
    operator=(const CGAL_Polygon_2<_Traits,_Container>& polygon)
    {
      d_container = polygon.d_container;
      return *this;
    }

    ~CGAL_Polygon_2()
      { }

#ifdef CGAL_CFG_NO_MEMBER_TEMPLATES
    // the following typedefs are required for Sun C++ 4.2
    typedef typename CGAL_STD::vector<Point_2>::const_iterator         v_ci;
    typedef typename CGAL_STD::vector<Point_2>::const_reverse_iterator v_cri;
    typedef typename CGAL_STD::vector<Point_2>::iterator               v_i;
    typedef typename CGAL_STD::vector<Point_2>::reverse_iterator       v_ri;
    typedef typename CGAL_STD::list<Point_2>::const_iterator           l_ci;
    typedef typename CGAL_STD::list<Point_2>::const_reverse_iterator   l_cri;
    typedef typename CGAL_STD::list<Point_2>::iterator                 l_i;
    typedef typename CGAL_STD::list<Point_2>::reverse_iterator         l_ri;

    CGAL_Polygon_2(v_ci first, v_ci last)
      { copy(first, last, back_inserter(d_container)); }
    CGAL_Polygon_2(v_cri first, v_cri last)
      { copy(first, last, back_inserter(d_container)); }
    CGAL_Polygon_2(v_i first, v_i last)
      { copy(first, last, back_inserter(d_container)); }
    CGAL_Polygon_2(v_ri first, v_ri last)
      { copy(first, last, back_inserter(d_container)); }
    CGAL_Polygon_2(l_ci first, l_ci last)
      { copy(first, last, back_inserter(d_container)); }
    CGAL_Polygon_2(l_cri first, l_cri last)
      { copy(first, last, back_inserter(d_container)); }
    CGAL_Polygon_2(l_i first, l_i last)
      { copy(first, last, back_inserter(d_container)); }
    CGAL_Polygon_2(l_ri first, l_ri last)
      { copy(first, last, back_inserter(d_container)); }
#else
    template <class InputIterator>
    CGAL_Polygon_2(InputIterator first, InputIterator last)
      { copy(first, last, back_inserter(d_container)); }
#endif

    //--------------------------------------------------------
    //             Operations
    //--------------------------------------------------------

    Vertex_iterator insert(Vertex_iterator position, const Point_2& x)
      { return d_container.insert(position,x); }

#ifndef CGAL_CFG_NO_MEMBER_TEMPLATES
    template <class InputIterator>
    void insert(Vertex_iterator position,
                InputIterator first,
                InputIterator last)
      { d_container.insert(position, first, last); }
#endif

    void push_back(const Point_2& x)
      { d_container.insert(d_container.end(), x); }

    void erase(Vertex_iterator position)
      { d_container.erase(position); }

    void erase(Vertex_iterator first, Vertex_iterator last)
      { d_container.erase(first,last); }

    void reverse_orientation()
    {
      if (size() <= 1)
        return;

      typename _Container::iterator i = d_container.begin();
      reverse(++i, d_container.end());
    }

    //--------------------------------------------------------
    //             Traversal of a polygon
    //--------------------------------------------------------

    Vertex_iterator vertices_begin()
      { return d_container.begin(); }

    Vertex_iterator vertices_end()
      { return d_container.end(); }

    Vertex_const_iterator vertices_begin() const
      { return d_container.begin(); }

    Vertex_const_iterator vertices_end() const
      { return d_container.end(); }

    Vertex_circulator vertices_circulator()
      { return Vertex_circulator(&d_container, d_container.begin()); }

    Vertex_const_circulator vertices_circulator() const
      { return Vertex_const_circulator(&d_container, d_container.begin()); }

    Edge_const_iterator edges_begin() const
      { return Edge_const_iterator(&d_container, d_container.begin()); }

    Edge_const_iterator edges_end() const
      { return Edge_const_iterator(&d_container, d_container.end()); }

    Edge_const_circulator edges_circulator() const
      { return Edge_const_circulator(vertices_circulator()); }

    //--------------------------------------------------------
    //             Predicates
    //--------------------------------------------------------

    bool is_simple() const
      { return CGAL_is_simple_2(d_container.begin(),
                                d_container.end(),
                                Traits()); }

    bool is_convex() const
    {
      return CGAL_is_convex_2(d_container.begin(),
                              d_container.end(),
                              Traits());
    }

    CGAL_Orientation orientation() const
    {
      CGAL_polygon_precondition(is_simple());
      return CGAL_orientation_2(d_container.begin(),
                                d_container.end(),
                                Traits());
    }

    CGAL_Oriented_side oriented_side(const Point_2& value) const
    {
      CGAL_polygon_precondition(is_simple());
      return CGAL_oriented_side_2(d_container.begin(),
                                  d_container.end(),
                                  value,
                                  Traits());
    }

    CGAL_Bounded_side bounded_side(const Point_2& value) const
    {
      CGAL_polygon_precondition(is_simple());
      return CGAL_bounded_side_2(d_container.begin(),
                                 d_container.end(),
                                 value,
                                 Traits());
    }

    CGAL_Bbox_2 bbox() const
      {  return CGAL_bbox_2(d_container.begin(), d_container.end()); }

    FT area() const
    {
      FT area(0);
      CGAL_area_2(d_container.begin(), d_container.end(), area, Traits());
      return area;
    }

    Vertex_const_iterator left_vertex() const
    {
      return CGAL_left_vertex_2(d_container.begin(),
                                d_container.end(),
                                Traits());
    }

    Vertex_iterator left_vertex()
    {
      return CGAL_left_vertex_2(d_container.begin(),
                                d_container.end(),
                                Traits());
    }

    Vertex_const_iterator right_vertex() const
    {
      return CGAL_right_vertex_2(d_container.begin(),
                                 d_container.end(),
                                 Traits());
    }

    Vertex_iterator right_vertex()
    {
      return CGAL_right_vertex_2(d_container.begin(),
                                 d_container.end(),
                                 Traits());
    }

    Vertex_const_iterator top_vertex() const
    {
      return CGAL_top_vertex_2(d_container.begin(),
                               d_container.end(),
                               Traits());
    }

    Vertex_iterator top_vertex()
    {
      return CGAL_top_vertex_2(d_container.begin(),
                               d_container.end(),
                               Traits());
    }

    Vertex_const_iterator bottom_vertex() const
    { 
      return CGAL_bottom_vertex_2(d_container.begin(),
                                  d_container.end(),
                                  Traits());
    }

    Vertex_iterator bottom_vertex()
    {
      return CGAL_bottom_vertex_2(d_container.begin(),
                                  d_container.end(),
                                  Traits());
    }

    bool is_counterclockwise_oriented() const
      { return orientation() == CGAL_COUNTERCLOCKWISE; }

    bool is_clockwise_oriented() const
      { return orientation() == CGAL_CLOCKWISE; }

    bool is_collinear_oriented() const
      { return orientation() == CGAL_COLLINEAR; }

    bool has_on_positive_side(const Point_2& q) const
      { return oriented_side(q) == CGAL_ON_POSITIVE_SIDE; }

    bool has_on_negative_side(const Point_2& q) const
      { return oriented_side(q) == CGAL_ON_NEGATIVE_SIDE; }

    bool has_on_boundary(const Point_2& q) const
      { return bounded_side(q) == CGAL_ON_BOUNDARY; }

    bool has_on_bounded_side(const Point_2& q) const
      { return bounded_side(q) == CGAL_ON_BOUNDED_SIDE; }

    bool has_on_unbounded_side(const Point_2& q) const
      { return bounded_side(q) == CGAL_ON_UNBOUNDED_SIDE; }

    //--------------------------------------------------------
    //             Random access methods
    //--------------------------------------------------------

#ifndef CGAL_CFG_NO_LAZY_INSTANTIATION
    const Point_2& vertex(int i) const
      { return *(vertices_begin() + i); }

    Point_2& vertex(int i)
      { return *(vertices_begin() + i); }

    const Point_2& operator[](int i) const
      { return vertex(i); }

    Point_2& operator[](int i)
      { return vertex(i); }

    Segment_2 edge(int i) const
      { return *(edges_begin() + i); }
#endif

    //--------------------------------------------------------
    //             Miscellaneous
    //--------------------------------------------------------

    int size() const
      { return d_container.size(); }

    bool is_empty() const
      { return d_container.empty(); }

    const _Container& container() const
      { return d_container; }

    bool identical(const CGAL_Polygon_2<_Traits,_Container> &q) const
      { return this == &q; }
};

//-----------------------------------------------------------------------//
//               Globally defined operators
//-----------------------------------------------------------------------//

template <class _Traits, class _Container1, class _Container2>
bool operator==( const CGAL_Polygon_2<_Traits,_Container1> &x,
                 const CGAL_Polygon_2<_Traits,_Container2> &y );

template <class _Traits, class _Container1, class _Container2>
inline
bool
operator!=(const CGAL_Polygon_2<_Traits,_Container1> &x,
           const CGAL_Polygon_2<_Traits,_Container2> &y);

#ifdef CGAL_REP_CLASS_DEFINED
#ifndef CGAL_POLYGON_TRAITS_2_H
#include <CGAL/Polygon_traits_2.h>
#endif // CGAL_POLYGON_TRAITS_2_H
template <class Transformation, class _Traits, class _Container>
CGAL_Polygon_2<_Traits,_Container>
CGAL_transform(const Transformation& t, const CGAL_Polygon_2<_Traits,_Container>& p);
#endif // CGAL_REP_CLASS_DEFINED

//-----------------------------------------------------------------------//
//               I/O
//-----------------------------------------------------------------------//

template <class _Traits, class _Container>
istream &operator>>(istream &is, CGAL_Polygon_2<_Traits,_Container>& p);

template <class _Traits, class _Container>
ostream &operator<<(ostream &os, const CGAL_Polygon_2<_Traits,_Container>& p);

//-----------------------------------------------------------------------//
//                         implementation
//-----------------------------------------------------------------------//

#ifdef CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION
#include <CGAL/Polygon_2.C>
#endif

template <class _Traits, class _Container1, class _Container2>
inline
bool
operator!=(const CGAL_Polygon_2<_Traits,_Container1> &x,
           const CGAL_Polygon_2<_Traits,_Container2> &y)
{
  return !(x==y);
}

#endif

