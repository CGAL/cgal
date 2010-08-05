// Copyright (c) 1997  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
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
// Author(s)     : Geert-Jan Giezeman <geert@cs.uu.nl>
//                 Wieger Wesselink

#ifndef CGAL_POLYGON_2_H
#define CGAL_POLYGON_2_H

#include <CGAL/basic.h>
#include <vector>
#include <list>
#include <iterator>

#include <CGAL/circulator.h>
#include <CGAL/enum.h>

#include <CGAL/Aff_transformation_2.h>

#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/Polygon_2/Polygon_2_vertex_circulator.h>
#include <CGAL/Polygon_2/Polygon_2_edge_iterator.h>
#include <CGAL/Polygon_2/Polygon_2_edge_circulator.h>

namespace CGAL {

template <class Traits_P, class Container_P
        = std::vector<typename Traits_P::Point_2> >
class Polygon_2 {

  public:

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

    typedef typename Container::iterator       Vertex_iterator;
    typedef typename Container::iterator       Vertex_const_iterator;
    //typedef typename Container::const_iterator Vertex_const_iterator; ??

    typedef Polygon_circulator<Container_P>    Vertex_const_circulator;
    typedef Vertex_const_circulator            Vertex_circulator;

    typedef Polygon_2_edge_iterator<Traits_P,Container_P>
            Edge_const_iterator;

    typedef Polygon_2_const_edge_circulator<Traits_P,Container_P>
            Edge_const_circulator;

    //--------------------------------------------------------
    //             Creation
    //--------------------------------------------------------

    Polygon_2(const Traits & p_traits = Traits()) : traits(p_traits) {}

    Polygon_2(const Polygon_2<Traits_P,Container_P>& polygon)
      : d_container(polygon.d_container), traits(polygon.traits) {}

    template <class InputIterator>
    Polygon_2(InputIterator first, InputIterator last,
              Traits p_traits = Traits())
        : d_container(), traits(p_traits)
    {
      // Sun STL switches off member templates for binary backward compat.
      std::copy(first, last, std::back_inserter(d_container));
    }

    //--------------------------------------------------------
    //             Operations
    //--------------------------------------------------------

    void set(Vertex_iterator pos, const Point_2& x)
     { *pos = x; }

    void set(Polygon_circulator<Container>const &pos, const Point_2& x)
     {
       *pos.mod_iterator() = x;
     }

    Vertex_iterator insert(Vertex_iterator pos, const Point_2& x)
      {
        return d_container.insert(pos,x);
      }

    Vertex_iterator insert(Vertex_circulator pos, const Point_2& x)
      {
        return d_container.insert(pos.mod_iterator(),x);
      }

    template <class InputIterator>
    void insert(Vertex_iterator pos,
                InputIterator first,
                InputIterator last)
      { d_container.insert(pos, first, last); }

    template <class InputIterator>
    void insert(Vertex_circulator pos,
                InputIterator first,
                InputIterator last)
      { d_container.insert(pos.mod_iterator(), first, last); }

    void push_back(const Point_2& x)
      { d_container.insert(d_container.end(), x); }

    void erase(Vertex_iterator pos)
      { d_container.erase(pos); }

    void erase(Vertex_circulator pos)
      { d_container.erase(pos.mod_iterator()); }

    void erase(Vertex_iterator first, Vertex_iterator last)
      {
        d_container.erase(first, last);
      }

    void clear()
    {
      d_container.clear();
    }

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

//    Vertex_iterator vertices_begin()
//      { return d_container.begin(); }

//    Vertex_iterator vertices_end()
//      { return d_container.end(); }

    Vertex_const_iterator vertices_begin() const
      { return const_cast<Polygon_2&>(*this).d_container.begin(); }

    Vertex_const_iterator vertices_end() const
      { return const_cast<Polygon_2&>(*this).d_container.end(); }

//    Vertex_const_circulator vertices_circulator() const
//      { return Vertex_const_circulator(&d_container, d_container.begin()); }

    Vertex_const_circulator vertices_circulator() const
      { 
        Polygon_2& self = const_cast<Polygon_2&>(*this);
        return Vertex_const_circulator(&self.d_container,
               self.d_container.begin());
      }

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
    {
      return is_simple_2(d_container.begin(),d_container.end(), traits);
    }

    bool is_convex() const
    {
      return is_convex_2(d_container.begin(),d_container.end(), traits);
    }

    Orientation orientation() const
    {
      return orientation_2(d_container.begin(), d_container.end(), traits);
    }

    Oriented_side oriented_side(const Point_2& value) const
    {
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
    {
      return bbox_2(d_container.begin(), d_container.end()); 
    }

    FT area() const
    {
      return polygon_area_2(d_container.begin(), d_container.end(), traits);
    }

    Vertex_const_iterator left_vertex() const
    {
       Polygon_2 &self = const_cast<Polygon_2&>(*this);
       return left_vertex_2(self.d_container.begin(),
                            self.d_container.end(), traits);
    }

    //Vertex_iterator left_vertex()
    //{
      //return left_vertex_2(d_container.begin(), d_container.end(), traits);
    //}

    Vertex_const_iterator right_vertex() const
    {
       Polygon_2 &self = const_cast<Polygon_2&>(*this);
       return right_vertex_2(self.d_container.begin(),
                             self.d_container.end(), traits);
    }

//    Vertex_iterator right_vertex()
//    {
//      return right_vertex_2(d_container.begin(), d_container.end(), traits);
//    }

    Vertex_const_iterator top_vertex() const
    {
       Polygon_2 &self = const_cast<Polygon_2&>(*this);
       return top_vertex_2(self.d_container.begin(),
                           self.d_container.end(), traits);
    }

//    Vertex_iterator top_vertex()
//    {
//      return top_vertex_2(d_container.begin(), d_container.end(), traits);
//    }

    Vertex_const_iterator bottom_vertex() const
    {
       Polygon_2 &self = const_cast<Polygon_2&>(*this);
       return bottom_vertex_2(self.d_container.begin(),
                              self.d_container.end(), traits);
    }

//    Vertex_iterator bottom_vertex()
//    {
//      return bottom_vertex_2(d_container.begin(), d_container.end(), traits);
//    }

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

  const Point_2& vertex(std::size_t i) const
      { return *(d_container.begin() + i); }

//    Point_2& vertex(std::size_t i)
//      { return *(d_container.begin() + i); }

  const Point_2& operator[](std::size_t i) const
      { return vertex(i); }

//    Point_2& operator[](std::size_t i)
//      { return vertex(i); }

  Segment_2 edge(std::size_t i) const
      { return *(edges_begin() + i); }

    //--------------------------------------------------------
    //             Miscellaneous
    //--------------------------------------------------------

  std::size_t size() const
      { return d_container.size(); }

    bool is_empty() const
      { return d_container.empty(); }

    const Container_P& container() const
      { return d_container; }

    bool identical(const Polygon_2<Traits_P,Container_P> &q) const
      { return this == &q; }


    Traits_P const &traits_member() const { return traits;}

  private:
    Container_P d_container;
    Traits_P traits;
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

template <class Transformation, class Traits_P, class Container_P>
Polygon_2<Traits_P,Container_P>
transform(const Transformation& t, const Polygon_2<Traits_P,Container_P>& p);

//-----------------------------------------------------------------------//
//               I/O
//-----------------------------------------------------------------------//

template <class Traits_P, class Container_P>
std::istream &operator>>(std::istream &is, Polygon_2<Traits_P,Container_P>& p);

template <class Traits_P, class Container_P>
std::ostream
&operator<<(std::ostream &os, const Polygon_2<Traits_P,Container_P>& p);

//-----------------------------------------------------------------------//
//                         implementation
//-----------------------------------------------------------------------//

} //namespace CGAL

#include <CGAL/Polygon_2/Polygon_2_impl.h>

namespace CGAL {

template <class Traits_P, class Container1_P, class Container2_P>
inline
bool
operator!=(const Polygon_2<Traits_P,Container1_P> &x,
           const Polygon_2<Traits_P,Container2_P> &y)
{
  return !(x==y);
}

} //namespace CGAL

#endif
