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

#include <CGAL/basic.h>

#ifdef CGAL_CFG_NO_MEMBER_TEMPLATES
#include <vector>
#include <list>
#endif

#include <iterator>

#include <CGAL/circulator.h>
#include <CGAL/enum.h>

#ifdef CGAL_REP_CLASS_DEFINED
#include <CGAL/Polygon_traits_2.h>
#include <CGAL/Aff_transformation_2.h>
#endif // CGAL_REP_CLASS_DEFINED

#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/Polygon_2_edge_iterator.h>
#include <CGAL/Polygon_2_edge_circulator.h>

CGAL_BEGIN_NAMESPACE

//-----------------------------------------------------------------------//
//                          Polygon_2
//-----------------------------------------------------------------------//


template <class Traits_P, class Container_P
        = std::vector<CGAL_TYPENAME_MSVC_NULL Traits_P::Point_2> >
class Polygon_2 {

  private:
    Container_P d_container;
    Traits_P traits;

  public:
    //--------------------------------------------------------
    //             Types
    //--------------------------------------------------------

    typedef Traits_P Traits;
    typedef Container_P Container;

    typedef typename Traits_P::FT FT;
    typedef typename Traits_P::Point_2 Point_2;
    typedef typename Traits_P::Segment_2 Segment_2;

    typedef typename Container_P::difference_type difference_type;
    typedef typename Container_P::value_type value_type;
    typedef typename Container_P::pointer pointer;
    typedef typename Container_P::reference reference;
    typedef typename Container_P::const_reference const_reference;


    //-------------------------------------------------------//
    // this intermediary step is required by Sun C++ 4.1
    typedef typename Container_P::iterator iterator;
    typedef typename Container_P::const_iterator const_iterator;
    //-------------------------------------------------------//

    typedef iterator Vertex_iterator;
    typedef const_iterator Vertex_const_iterator;

    typedef Bidirectional_circulator_from_container<Container_P>
            Vertex_circulator;

    typedef Bidirectional_const_circulator_from_container<Container_P>
            Vertex_const_circulator;

    typedef Polygon_2_edge_iterator<Traits_P,Container_P>
            Edge_const_iterator;

    typedef Polygon_2_const_edge_circulator<Traits_P,Container_P>
            Edge_const_circulator;

    //--------------------------------------------------------
    //             Creation
    //--------------------------------------------------------

    Polygon_2(Traits p_traits = Traits()) : traits(p_traits)
      { }

    Polygon_2(const Polygon_2<Traits_P,Container_P>& polygon)
      : d_container(polygon.d_container), traits(polygon.traits) { }

    Polygon_2<Traits_P,Container_P>&
    operator=(const Polygon_2<Traits_P,Container_P>& polygon)
    {
      d_container = polygon.d_container;
      traits = polygon.traits;
      return *this;
    }

    ~Polygon_2()
      { }

#ifdef CGAL_CFG_NO_MEMBER_TEMPLATES
    // the following typedefs are required for Sun C++ 4.2
    typedef typename std::vector<Point_2>::const_iterator         v_ci;
    typedef typename std::vector<Point_2>::const_reverse_iterator v_cri;
    typedef typename std::vector<Point_2>::iterator               v_i;
    typedef typename std::vector<Point_2>::reverse_iterator       v_ri;
    typedef typename std::list<Point_2>::const_iterator           l_ci;
    typedef typename std::list<Point_2>::const_reverse_iterator   l_cri;
    typedef typename std::list<Point_2>::iterator                 l_i;
    typedef typename std::list<Point_2>::reverse_iterator         l_ri;

    Polygon_2(v_ci first, v_ci last)
      { copy(first, last, back_inserter(d_container)); }
    Polygon_2(v_cri first, v_cri last)
      { copy(first, last, back_inserter(d_container)); }
    Polygon_2(v_i first, v_i last)
      { copy(first, last, back_inserter(d_container)); }
    Polygon_2(v_ri first, v_ri last)
      { copy(first, last, back_inserter(d_container)); }
    Polygon_2(l_ci first, l_ci last)
      { copy(first, last, back_inserter(d_container)); }
    Polygon_2(l_cri first, l_cri last)
      { copy(first, last, back_inserter(d_container)); }
    Polygon_2(l_i first, l_i last)
      { copy(first, last, back_inserter(d_container)); }
    Polygon_2(l_ri first, l_ri last)
      { copy(first, last, back_inserter(d_container)); }
#else
   template <class InputIterator>
    Polygon_2(InputIterator first, InputIterator last,
            Traits p_traits = Traits())
        : d_container(first,last), traits(p_traits) {}
//      { std::copy(first, last, std::back_inserter(d_container)); }

/*
    template <class Circulator>
    Polygon_2(Circulator start)
	{
	    if (start != NULL) {
		Circulator cur = start;
		do {
		    d_container.push_back(*cur); ++cur;
		} while (cur != start);
	    }
	}
*/
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

      typename Container_P::iterator i = d_container.begin();
      std::reverse(++i, d_container.end());
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
      { return is_simple_2(d_container.begin(), d_container.end(), traits); }

    bool is_convex() const
    {
      return is_convex_2(d_container.begin(), d_container.end(), traits);
    }

    Orientation orientation() const
    {
      CGAL_polygon_precondition(is_simple());
      return orientation_2(d_container.begin(), d_container.end(), traits);
    }

    Oriented_side oriented_side(const Point_2& value) const
    {
      CGAL_polygon_precondition(is_simple());
      return oriented_side_2(d_container.begin(), d_container.end(),
                                  value, traits);
    }

    Bounded_side bounded_side(const Point_2& value) const
    {
      CGAL_polygon_precondition(is_simple());
      return bounded_side_2(d_container.begin(), d_container.end(),
                                 value, traits);
    }

    Bbox_2 bbox() const
      {  return bbox_2(d_container.begin(), d_container.end()); }

    FT area() const
    {
      FT area(0);
      area_2(d_container.begin(), d_container.end(), area, traits);
      return area;
    }

    Vertex_const_iterator left_vertex() const
    {
      return left_vertex_2(d_container.begin(), d_container.end(), traits);
    }

    Vertex_iterator left_vertex()
    {
      return left_vertex_2(d_container.begin(), d_container.end(), traits);
    }

    Vertex_const_iterator right_vertex() const
    {
      return right_vertex_2(d_container.begin(), d_container.end(), traits);
    }

    Vertex_iterator right_vertex()
    {
      return right_vertex_2(d_container.begin(), d_container.end(), traits);
    }

    Vertex_const_iterator top_vertex() const
    {
      return top_vertex_2(d_container.begin(), d_container.end(), traits);
    }

    Vertex_iterator top_vertex()
    {
      return top_vertex_2(d_container.begin(), d_container.end(), traits);
    }

    Vertex_const_iterator bottom_vertex() const
    { 
      return bottom_vertex_2(d_container.begin(), d_container.end(), traits);
    }

    Vertex_iterator bottom_vertex()
    {
      return bottom_vertex_2(d_container.begin(), d_container.end(), traits);
    }

    bool is_counterclockwise_oriented() const
      { return orientation() == COUNTERCLOCKWISE; }

    bool is_clockwise_oriented() const
      { return orientation() == CLOCKWISE; }

    bool is_collinear_oriented() const
      { return orientation() == COLLINEAR; }

    bool has_on_positive_side(const Point_2& q) const
      { return oriented_side(q) == ON_POSITIVE_SIDE; }

    bool has_on_negative_side(const Point_2& q) const
      { return oriented_side(q) == ON_NEGATIVE_SIDE; }

    bool has_on_boundary(const Point_2& q) const
      { return bounded_side(q) == ON_BOUNDARY; }

    bool has_on_bounded_side(const Point_2& q) const
      { return bounded_side(q) == ON_BOUNDED_SIDE; }

    bool has_on_unbounded_side(const Point_2& q) const
      { return bounded_side(q) == ON_UNBOUNDED_SIDE; }

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

    const Container_P& container() const
      { return d_container; }

    bool identical(const Polygon_2<Traits_P,Container_P> &q) const
      { return this == &q; }


    Traits_P const &traits_member() const { return traits;}
};

//-----------------------------------------------------------------------//
//               Globally defined operators
//-----------------------------------------------------------------------//

template <class Traits_P, class Container1_P, class Container2_P>
bool operator==( const Polygon_2<Traits_P,Container1_P> &x,
                 const Polygon_2<Traits_P,Container2_P> &y );

template <class Traits_P, class Container1_P, class Container2_P>
inline
bool
operator!=(const Polygon_2<Traits_P,Container1_P> &x,
           const Polygon_2<Traits_P,Container2_P> &y);

#ifdef CGAL_REP_CLASS_DEFINED

CGAL_END_NAMESPACE

#    include <CGAL/Polygon_traits_2.h>

CGAL_BEGIN_NAMESPACE

template <class Transformation, class Traits_P, class Container_P>
Polygon_2<Traits_P,Container_P>
transform(const Transformation& t, const Polygon_2<Traits_P,Container_P>& p);

#endif // CGAL_REP_CLASS_DEFINED

//-----------------------------------------------------------------------//
//               I/O
//-----------------------------------------------------------------------//

template <class Traits_P, class Container_P>
std::istream &operator>>(std::istream &is, Polygon_2<Traits_P,Container_P>& p);

template <class Traits_P, class Container_P>
std::ostream &operator<<(std::ostream &os, const Polygon_2<Traits_P,Container_P>& p);

//-----------------------------------------------------------------------//
//                         implementation
//-----------------------------------------------------------------------//

CGAL_END_NAMESPACE

#ifdef CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION
#include <CGAL/Polygon_2.C>
#endif

CGAL_BEGIN_NAMESPACE

template <class Traits_P, class Container1_P, class Container2_P>
inline
bool
operator!=(const Polygon_2<Traits_P,Container1_P> &x,
           const Polygon_2<Traits_P,Container2_P> &y)
{
  return !(x==y);
}

CGAL_END_NAMESPACE

#endif

