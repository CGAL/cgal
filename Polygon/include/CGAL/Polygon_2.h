// Copyright (c) 1997  
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
// Author(s)     : Geert-Jan Giezeman <geert@cs.uu.nl>
//                 Wieger Wesselink

/*!
  \file Polygon_2.h
 */

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

/// \ingroup PkgPolygon2
/// The class Polygon_2 implements polygons. The Polygon_2 is
/// parameterized by a traits class and a container class.  The latter
/// can be any class that fulfills the requirements for an STL
/// container. It defaults to the std::vector class.
///
/// \cgalHeading{Implementation}
///
/// The methods `is_simple()`, `is_convex()`, `orientation()`,
/// `oriented_side()`, `bounded_side()`, `bbox()`, `area()`, `left_vertex()`,
/// `right_vertex()`, `top_vertex()` and `bottom_vertex()` are all
/// implemented using the algorithms on sequences of 2D points.  See
/// the corresponding global functions for information about which
/// algorithms were used and what complexity they have.
///

template <class Traits_P, class Container_P
        = std::vector<typename Traits_P::Point_2> >
class Polygon_2 {

  public:

    /// \name Types
    /// @{

    /// The traits type.
    typedef Traits_P Traits;
    /// The container type.
    typedef Container_P Container;

    /// The number type of the coordinates of the points of the polygon.
    typedef typename Traits_P::FT FT;
    /// The point type of the polygon.
    typedef typename Traits_P::Point_2 Point_2;
    /// The type of a segment between two points of the polygon.
    typedef typename Traits_P::Segment_2 Segment_2;

    /// @}

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
    typedef typename Container::iterator       Vertex_const_iterator;
    typedef Polygon_circulator<Container_P>    Vertex_const_circulator;

    /// \name Iterators
    ///
    /// The following types denote iterators that allow to traverse
    /// the vertices and edges of a polygon.  Since 
    /// a polygon can be viewed as a circular as well as a
    /// linear data structure both circulators and iterators are
    /// defined.  
    ///
    /// \note At least conceptually, the circulators and iterators are
    /// non-mutable.  The enforcement depends on preprocessor flags. 
    ///
    /// \note The iterator category is in all cases bidirectional, except
    /// for Vertex_iterator, which has the same iterator category as
    /// `Container::iterator`. In fact all of them should have
    /// the same iterator category as `Container::iterator`. However,
    /// due to compiler problems this is currently not possible.
    ///
    /// @{

    /// 
    typedef typename Container::iterator       Vertex_iterator;


    //typedef typename Container::const_iterator Vertex_const_iterator; ??

#ifdef DOXYGEN_RUNNING
  typedef unspecified_type Vertex_circulator;
    typedef unspecified_type Edge_const_iterator;

    typedef unspecified_type Edge_const_circulator;
#else 
    typedef Vertex_const_circulator            Vertex_circulator;
    /// 
    typedef Polygon_2_edge_iterator<Traits_P,Container_P>
            Edge_const_iterator;

    /// 
    typedef Polygon_2_const_edge_circulator<Traits_P,Container_P>
            Edge_const_circulator;
#endif // DOXYGEN_RUNNING    
    /// @}

    /// \name Creation
    /// @{

    /// Creates an empty polygon.
    Polygon_2(const Traits & p_traits = Traits()) : traits(p_traits) {}

    /// Copy constructor.
    Polygon_2(const Polygon_2<Traits_P,Container_P>& polygon)
      : d_container(polygon.d_container), traits(polygon.traits) {}

    /// Introduces a polygon with vertices from the sequence
    /// defined by the range \c [first,last).
    /// The value type of \c InputIterator must be \c Point_2.
    template <class InputIterator>
    Polygon_2(InputIterator first, InputIterator last,
              Traits p_traits = Traits())
        : d_container(), traits(p_traits)
    {
      // Sun STL switches off member templates for binary backward compat.
      std::copy(first, last, std::back_inserter(d_container));
    }

    /// @}

    /// \name Modifiers
    /// @{

    /// Acts as `*i = q`, except that that would be illegal because
    /// the iterator is not mutable.
    void set(Vertex_iterator i, const Point_2& q)
     { *i = q; }

    void set(Polygon_circulator<Container>const &i, const Point_2& q)
     {
       *i.mod_iterator() = q;
     }

    /// Inserts the vertex `q` before `i`. The return value points to
    /// the inserted vertex.
    Vertex_iterator insert(Vertex_iterator i, const Point_2& q)
      {
        return d_container.insert(i,q);
      }

    Vertex_iterator insert(Vertex_circulator i, const Point_2& q)
      {
        return d_container.insert(i.mod_iterator(),q);
      }

    /// Inserts the vertices in the range `[first, last)`
    /// before `i`.  The value type of points in the range
    /// `[first,last)} must be \ccStyle{Point_2`.
    template <class InputIterator>
    void insert(Vertex_iterator i,
                InputIterator first,
                InputIterator last)
      { d_container.insert(i, first, last); }

    template <class InputIterator>
    void insert(Vertex_circulator i,
                InputIterator first,
                InputIterator last)
      { d_container.insert(i.mod_iterator(), first, last); }

    /// Has the same semantics as `p.insert(p.vertices_end(), q)`.
    void push_back(const Point_2& x)
      { d_container.insert(d_container.end(), x); }

    /// Erases the vertex pointed to by `i`.
    Vertex_iterator erase(Vertex_iterator i)
      {
        return d_container.erase(i);
      }

    Vertex_circulator erase(Vertex_circulator i)
      {
        return Vertex_circulator(&d_container,
                                 d_container.erase(i.mod_iterator()));
      }

    /// Erases the vertices in the range `[first, last)`.
    Vertex_iterator erase(Vertex_iterator first, Vertex_iterator last)
      {
        return d_container.erase(first, last);
      }

    /// Erases the vertices in the range `[first, last)`.
    void clear()
    {
      d_container.clear();
    }

    /// Reverses the orientation of the polygon. The vertex pointed to
    ///  by `p.vertices_begin()` remains the same.
    void reverse_orientation()
    {
      if (size() <= 1)
        return;
      typename Container_P::iterator i = d_container.begin();
      std::reverse(++i, d_container.end());
    }

    /// @}

    /// \name Access Functions 
    /// The following methods of the class Polygon_2
    /// return circulators and iterators that allow to traverse the
    /// vertices and edges.
    /// @{

    /// Returns a constant iterator that allows to traverse the
    /// vertices of the polygon.
    Vertex_const_iterator vertices_begin() const
      { return const_cast<Polygon_2&>(*this).d_container.begin(); }

    /// Returns the corresponding past-the-end iterator.
    Vertex_const_iterator vertices_end() const
      { return const_cast<Polygon_2&>(*this).d_container.end(); }

//    Vertex_const_circulator vertices_circulator() const
//      { return Vertex_const_circulator(&d_container, d_container.begin()); }

    /// Returns a mutable circulator that allows to traverse the
    /// vertices of the polygon.
    Vertex_const_circulator vertices_circulator() const
      { 
        Polygon_2& self = const_cast<Polygon_2&>(*this);
        return Vertex_const_circulator(&self.d_container,
               self.d_container.begin());
      }

    /// Returns a non-mutable iterator that allows to traverse the
    /// edges of the polygon.
    Edge_const_iterator edges_begin() const
      { return Edge_const_iterator(&d_container, d_container.begin()); }

    /// Returns the corresponding past-the-end iterator.
    Edge_const_iterator edges_end() const
      { return Edge_const_iterator(&d_container, d_container.end()); }

    /// Returns a non-mutable circulator that allows to traverse the
    /// edges of the polygon.
    Edge_const_circulator edges_circulator() const
      { return Edge_const_circulator(vertices_circulator()); }

    /// @}

    /// \name Predicates
    /// @{

    /// Returns whether this is a simple polygon.
    bool is_simple() const
    {
      return is_simple_2(d_container.begin(),d_container.end(), traits);
    }

    /// Returns whether this is convex.
    bool is_convex() const
    {
      return is_convex_2(d_container.begin(),d_container.end(), traits);
    }

    /// Returns the orientation. If the number of vertices
    /// `p.size() < 3` then \c COLLINEAR is returned.
    /// \pre `p.is_simple()`.
    Orientation orientation() const
    {
      return orientation_2(d_container.begin(), d_container.end(), traits);
    }

    /// Returns `POSITIVE_SIDE`, or `NEGATIVE_SIDE`,
    /// or `ON_ORIENTED_BOUNDARY`, depending on where point
    /// `q` is. 
    /// \pre `p.is_simple()`.
    Oriented_side oriented_side(const Point_2& value) const
    {
      return oriented_side_2(d_container.begin(), d_container.end(),
                                  value, traits);
    }

    /// Returns the symbolic constant `ON_BOUNDED_SIDE`,
    /// `ON_BOUNDARY` or `ON_UNBOUNDED_SIDE`,
    /// depending on where point `q` is. \pre
    /// `p.is_simple()`.
    Bounded_side bounded_side(const Point_2& value) const
    {
      CGAL_polygon_precondition(is_simple());
      return bounded_side_2(d_container.begin(), d_container.end(),
                                 value, traits);
    }

    /// Returns the smallest bounding box containing this polygon.
    Bbox_2 bbox() const
    {
      return bbox_2(d_container.begin(), d_container.end()); 
    }

    /// Returns the signed area of the polygon. This means that the
    /// area is positive for counter clockwise polygons and negative
    /// for clockwise polygons.
    FT area() const
    {
      return polygon_area_2(d_container.begin(), d_container.end(), traits);
    }

    /// Returns the leftmost vertex of the polygon with the smallest
    /// `x`-coordinate.
    Vertex_const_iterator left_vertex() const
    {
       Polygon_2 &self = const_cast<Polygon_2&>(*this);
       return left_vertex_2(self.d_container.begin(),
                            self.d_container.end(), traits);
    }

    /// Returns the rightmost vertex of the polygon with the largest
    /// `x`-coordinate.
    Vertex_const_iterator right_vertex() const
    {
       Polygon_2 &self = const_cast<Polygon_2&>(*this);
       return right_vertex_2(self.d_container.begin(),
                             self.d_container.end(), traits);
    }

    /// Returns topmost vertex of the polygon with the largest
    /// `y`-coordinate.
    Vertex_const_iterator top_vertex() const
    {
       Polygon_2 &self = const_cast<Polygon_2&>(*this);
       return top_vertex_2(self.d_container.begin(),
                           self.d_container.end(), traits);
    }

    /// Returns the bottommost vertex of the polygon with the
    /// smallest `y`-coordinate.
    Vertex_const_iterator bottom_vertex() const
    {
       Polygon_2 &self = const_cast<Polygon_2&>(*this);
       return bottom_vertex_2(self.d_container.begin(),
                              self.d_container.end(), traits);
    }

    /// @}


    /// \name 
    /// For convenience we provide the following Boolean functions:
    /// @{

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

    /// @}


    /// \name Random Access Methods
    /// @{

    /// Returns a (const) reference to the `i`-th vertex.
    const Point_2& vertex(std::size_t i) const
      { return *(d_container.begin() + i); }


    /// Returns a (const) reference to the `i`-th vertex.
    const Point_2& operator[](std::size_t i) const
      { return vertex(i); }

    /// Returns the `i`-th edge.
    Segment_2 edge(std::size_t i) const
      { return *(edges_begin() + i); }

    /// @}

    /// \name Miscellaneous
    /// @{

    /// Returns the number of vertices of the polygon.
    std::size_t size() const
      { return d_container.size(); }

    /// Returns `size() == 0`.
    bool is_empty() const
      { return d_container.empty(); }

    /// Returns a const reference to the sequence of vertices of the polygon.
    const Container_P& container() const
      { return d_container; }

    /// @}

    bool identical(const Polygon_2<Traits_P,Container_P> &q) const
      { return this == &q; }

    Traits_P const &traits_member() const { return traits;}

  private:
    Container_P d_container;
    Traits_P traits;
};

/// @} // polygon_2

/// \name Global Operators
/// @{

/// Test for equality: two polygons are equal iff there exists a
/// cyclic permutation of the vertices of `p2` such that they are
/// equal to the vertices of `p1`. Note that the template argument
/// `Container` of `p1` and `p2` may be different.
/// \memberof Polygon_2
template <class Traits_P, class Container1_P, class Container2_P>
bool operator==( const Polygon_2<Traits_P,Container1_P> &p1,
                 const Polygon_2<Traits_P,Container2_P> &p2 );

/// Test for inequality. 
/// \memberof Polygon_2
template <class Traits_P, class Container1_P, class Container2_P>
inline
bool
operator!=(const Polygon_2<Traits_P,Container1_P> &p1,
           const Polygon_2<Traits_P,Container2_P> &p2);

/// Returns the image of the polygon \c p under the transformation \c t.
/// \memberof Polygon_2
template <class Transformation, class Traits_P, class Container_P>
Polygon_2<Traits_P,Container_P>
transform(const Transformation& t, const Polygon_2<Traits_P,Container_P>& p);

/// @} // global operators

/// \name I/O
/// The information output in the `std::iostream` is the number of points
/// followed by the output of the coordinates of the vertices.
/// @{

/// Inserts the polygon `p` into the stream `os`. \pre The insert
/// operator must be defined for `Point_2`.
/// \memberof Polygon_2
template <class Traits_P, class Container_P>
std::istream &operator>>(std::istream &is, Polygon_2<Traits_P,Container_P>& p);

/// Reads a polygon from stream `is` and assigns it
/// to `p`. \pre The extract operator must be defined for `Point_2`.
/// \memberof Polygon_2
template <class Traits_P, class Container_P>
std::ostream
&operator<<(std::ostream &os, const Polygon_2<Traits_P,Container_P>& p);

/// @} // IO

} //namespace CGAL

//-----------------------------------------------------------------------//
//                         implementation
//-----------------------------------------------------------------------//

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
