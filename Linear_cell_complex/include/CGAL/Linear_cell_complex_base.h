// Copyright (c) 2011 CNRS and LIRIS' Establishments (France).
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
// SPDX-License-Identifier: LGPL-3.0+
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
#ifndef CGAL_LINEAR_CELL_COMPLEX_BASE_H
#define CGAL_LINEAR_CELL_COMPLEX_BASE_H 1

#include <CGAL/Combinatorial_map_functors.h>
#include <CGAL/internal/Combinatorial_map_internal_functors.h>
#include <CGAL/Linear_cell_complex_operations.h>
#include <CGAL/Origin.h>

#include<map>

namespace CGAL {

  // Tags allowing to make some template specialization for CMap and GMap,
  // if necessary.
  struct Combinatorial_map_tag;
  struct Generalized_map_tag;


  /** @file Linear_cell_complex_base.h
   * Definition of a linear cell complex, i.e. a combinatorial data structure
   * (cmap or gmap) with points associated to all vertices.
   */

  /**  Linear_cell_complex_base class.
   * The Linear_cell_complex a nD object with linear geometry, ie
   * an nD map with point associated to each vertex.
   */
  template < unsigned int d_, unsigned int ambient_dim,
             class Traits_,
             class Items_,
             class Alloc_,
             template<unsigned int,class,class,class,class>
             class Map,
             class Refs,
             class Storage_>
  class Linear_cell_complex_base:
    public Map<d_, Refs, Items_, Alloc_, Storage_>
  {
  public:
    typedef Linear_cell_complex_base<d_, ambient_dim,
                                Traits_, Items_, Alloc_, Map,
                                Refs, Storage_>  Self;
    typedef Map<d_, Refs, Items_, Alloc_, Storage_> Base;

    typedef Traits_ Traits;
    typedef Items_  Items;
    typedef Alloc_  Alloc;
    typedef Storage_ Storage;

    static const unsigned int ambient_dimension = ambient_dim;
    static const unsigned int dimension = Base::dimension;

    typedef typename Storage::Dart_handle       Dart_handle;
    typedef typename Storage::Dart_const_handle Dart_const_handle;
    typedef typename Storage::Helper            Helper;

    typedef typename Storage::Point  Point;
    typedef typename Storage::Vector Vector;
    typedef typename Storage::FT     FT;

    typedef typename Base::Dart_range Dart_range;

    /// Typedef for attributes
    template<int i>
    struct Attribute_type: public Base::template Attribute_type<i>
    {};
    template<int i>
    struct Attribute_handle: public Base::template Attribute_handle<i>
    {};
    template<int i>
    struct Attribute_const_handle:
      public Base::template Attribute_const_handle<i>
    {};
    template<int i>
    struct Attribute_range: public Base::template Attribute_range<i>
    {};
    template<int i>
    struct Attribute_const_range:
      public Base::template Attribute_const_range<i>
    {};

    typedef typename Base::template Attribute_type<0>::type Vertex_attribute;
    typedef typename Base::template Attribute_handle<0>::type
    Vertex_attribute_handle;
    typedef typename Base::template Attribute_const_handle<0>::type
    Vertex_attribute_const_handle;

    typedef typename Base::template Attribute_range<0>::type
    Vertex_attribute_range;
    typedef typename Base::template Attribute_const_range<0>::type
    Vertex_attribute_const_range;

    typedef typename Base::size_type size_type;
    typedef typename Base::Use_index Use_index;
    typedef typename Base::Exception_no_more_available_mark
    Exception_no_more_available_mark;

    /// To use previous definition of create_dart methods.
    using Base::create_dart;
    using Base::attribute;
    using Base::null_handle;
    using Base::point_of_vertex_attribute;
    using Base::other_extremity;
    using Base::darts;
    
    using Base::are_attributes_automatically_managed;
    using Base::mark;
    using Base::is_marked;
    using Base::unmark;
    using Base::free_mark;
    using Base::get_new_mark;
    using Base::next;
    using Base::previous;
    using Base::opposite;
    using Base::is_next_exist;
    using Base::is_previous_exist;
    
    using Base::make_segment;
    using Base::make_triangle;
    using Base::make_quadrangle;

    using Base::insert_cell_0_in_cell_1;
    using Base::insert_cell_0_in_cell_2;
    using Base::insert_dangling_cell_1_in_cell_2;
    using Base::insert_cell_1_in_cell_2;
    using Base::insert_cell_2_in_cell_3;

    Linear_cell_complex_base() : Base()
    {}

    /** Copy the given linear cell complex into *this.
     *  Note that both LCC can have different dimensions and/or non void attributes.
     *  @param alcc the linear cell complex to copy.
     *  @post *this is valid.
     */
    Linear_cell_complex_base(const Self& alcc) : Base(alcc)
    {}

    template < class LCC2 >
    Linear_cell_complex_base(const LCC2& alcc) : Base(alcc)
    {}

    template < class LCC2, typename Converters >
    Linear_cell_complex_base(const LCC2& alcc, Converters& converters) :
      Base(alcc, converters)
    {}

    template < class LCC2, typename Converters, typename DartInfoConverter >
    Linear_cell_complex_base(const LCC2& alcc, Converters& converters,
                             const DartInfoConverter& dartinfoconverter) :
      Base(alcc, converters, dartinfoconverter)
    {}

    template < class LCC2, typename Converters, typename DartInfoConverter,
               typename Pointconverter >
    Linear_cell_complex_base(const LCC2& alcc, Converters& converters,
                             const DartInfoConverter& dartinfoconverter,
                             const Pointconverter& pointconverter) :
      Base(alcc, converters, dartinfoconverter, pointconverter)
    {}

    /** Affectation operation. Copies one map to the other.
     * @param amap a lcc.
     * @return A copy of that lcc.
     */
    Self & operator= (const Self & alcc)
    {
      if (this!=&alcc)
      {
        Self tmp(alcc);
        this->swap(tmp);
      }
      return *this;
    }

    /** Create a vertex attribute.
     * @return an handle on the new attribute.
     */
#ifndef CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES
    template<typename ...Args>
    Vertex_attribute_handle create_vertex_attribute(const Args&... args)
    { return Base::template create_attribute<0>(args...); }
#else
    Vertex_attribute_handle create_vertex_attribute()
    { return Base::template create_attribute<0>(); }

    template<typename T1>
    Vertex_attribute_handle create_vertex_attribute(const T1& t1)
    { return Base::template create_attribute<0>(t1); }

    template<typename T1, typename T2>
    Vertex_attribute_handle create_vertex_attribute
    (const T1& t1, const T2 &t2)
    { return Base::template create_attribute<0>(t1, t2); }

    template<typename T1, typename T2, typename T3>
    Vertex_attribute_handle create_vertex_attribute
    (const T1& t1, const T2 &t2, const T3 &t3)
    { return Base::template create_attribute<0>(t1, t2, t3); }

    template<typename T1, typename T2, typename T3, typename T4>
    Vertex_attribute_handle create_vertex_attribute
    (const T1& t1, const T2 &t2, const T3 &t3, const T4 &t4)
    { return Base::template create_attribute<0>(t1, t2, t3, t4); }

    template<typename T1, typename T2, typename T3, typename T4, typename T5>
    Vertex_attribute_handle create_vertex_attribute
    (const T1& t1, const T2 &t2, const T3 &t3, const T4 &t4,
     const T5 &t5)
    { return Base::template create_attribute<0>(t1, t2, t3, t4, t5); }

    template<typename T1, typename T2, typename T3, typename T4, typename T5,
             typename T6>
    Vertex_attribute_handle create_vertex_attribute
    (const T1& t1, const T2 &t2, const T3 &t3, const T4 &t4,
     const T5 &t5, const T6 &t6)
    { return Base::template create_attribute<0>(t1, t2, t3, t4, t5, t6); }

    template<typename T1, typename T2, typename T3, typename T4, typename T5,
             typename T6, typename T7>
    Vertex_attribute_handle create_vertex_attribute
    (const T1& t1, const T2 &t2, const T3 &t3, const T4 &t4,
     const T5 &t5, const T6 &t6, const T7 &t7)
    { return Base::template create_attribute<0>(t1, t2, t3, t4, t5, t6, t7); }

    template<typename T1, typename T2, typename T3, typename T4, typename T5,
             typename T6, typename T7, typename T8>
    Vertex_attribute_handle create_vertex_attribute
    (const T1& t1, const T2 &t2, const T3 &t3, const T4 &t4,
     const T5 &t5, const T6 &t6, const T7 &t7, const T8 &t8)
    { return Base::template create_attribute<0>(t1, t2, t3, t4, t5, t6, t7,
                                                t8); }

    template<typename T1, typename T2, typename T3, typename T4, typename T5,
             typename T6, typename T7, typename T8, typename T9>
    Vertex_attribute_handle create_vertex_attribute
    (const T1& t1, const T2 &t2, const T3 &t3, const T4 &t4,
     const T5 &t5, const T6 &t6, const T7 &t7, const T8 &t8, const T9 &t9)
    { return Base::template create_attribute<0>(t1, t2, t3, t4, t5, t6, t7,
                                                t8, t9); }
#endif // CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES

    /**
     * Create a new dart associated with an handle through an attribute.
     * @param ahandle the point handle to associated with the dart.
     * @return a Dart_handle on the new dart.
     */
    Dart_handle create_dart(Vertex_attribute_handle ahandle)
    {
      Dart_handle res = create_dart();
      set_vertex_attribute_of_dart(res,ahandle);
      return res;
    }

    /** Create a new dart associated with a point.
     * @param apoint the point to associated with the dart.
     * @return a Dart_handle on the new dart.
     */
    Dart_handle create_dart(const Point& apoint)
    { return create_dart(create_vertex_attribute(apoint)); }

    /** Erase a given vertex attribute.
     * @param ahandle the handle to the vertex attribute to erase.
     */
    void erase_vertex_attribute(Vertex_attribute_handle ahandle)
    { Base::template erase_attribute<0>(ahandle); }

    /** Set the vertex attribute of the given dart.
     * @param adart a dart.
     * @param ah the attribute to set.
     */
    void set_vertex_attribute_of_dart(Dart_handle adart,
                                      Vertex_attribute_handle ah)
    {
      return CGAL::internal::Set_i_attribute_of_dart_functor<Self, 0>::
          run(*this, adart, ah);
    }

    /** Set the vertex attribute of all the darts of the vertex.
     * @param adart a dart of the vertex.
     * @param ah the attribute to set.
     */
    void set_vertex_attribute(Dart_handle adart,
                              Vertex_attribute_handle ah)
    { return CGAL::Set_i_attribute_functor<Self, 0>::run(*this, adart, ah); }

    /// @return the Vertex_attribute_range for all vertex_attributes.
    Vertex_attribute_range& vertex_attributes()
    { return this->template attributes<0>(); }

    /// @return the Vertex_attribute_const_range for all vertex_attributes.
    Vertex_attribute_const_range& vertex_attributes() const
    { return this->template attributes<0>(); }

    /// @return the size of the vertex_attribute container.
    typename Base::size_type number_of_vertex_attributes() const
    { return Base::template number_of_attributes<0>(); }

    /// Get the vertex_attribute associated with a dart.
    /// @param a dart
    /// @return the vertex_attribute.
    Vertex_attribute_handle vertex_attribute(Dart_handle adart)
    { return this->template attribute<0>(adart); }

    /// Get the vertex_attribute associated with a const dart.
    /// @param a dart
    /// @return the vertex_const_attribute.
    Vertex_attribute_const_handle
    vertex_attribute(Dart_const_handle adart) const
    { return this->template attribute<0>(adart); }

    /// Get the point associated with a dart.
    /// @param a dart
    /// @return the point.
    Point& point(Dart_handle adart)
    {
      CGAL_assertion(this->template attribute<0>(adart)!=null_handle );
      return point_of_vertex_attribute(this->template attribute<0>(adart));
    }

    /// Get the point associated with a const dart.
    /// @param a dart
    /// @return the point.
    const Point& point(Dart_const_handle adart) const
    {
      CGAL_assertion(this->template attribute<0>(adart)!=null_handle );
      return point_of_vertex_attribute(this->template attribute<0>(adart));
    }

    /** Test if the lcc is valid.
     * A Linear_cell_complex is valid if it is a valid Combinatorial_map with
     * an attribute associated to each dart.
     * @return true iff the map is valid.
     */
    bool is_valid() const
    {
      bool valid = Base::is_valid();
      for (typename Dart_range::const_iterator it(darts().begin()),
             itend(darts().end()); valid && it != itend; ++it)
      {
        if ( vertex_attribute(it)==null_handle )
        {
          std::cerr << "Map not valid: dart "<<&(*it)
                    <<" does not have a vertex."<< std::endl;
          valid = false;
        }
      }
      return valid;
    }

    /** validate the lcc
     */
    void correct_invalid_attributes()
    {
      // Copy of the code in Map::correct_invalid_attributes() to avoid
      // 2 iterations through the darts of the map.

      std::vector<size_type> marks(dimension+1);
      for ( unsigned int i=0; i<=dimension; ++i)
        marks[i] = Base::INVALID_MARK;

      Helper::template
        Foreach_enabled_attributes<Reserve_mark_functor<Self> >::
          run(*this, marks);

      for ( typename Dart_range::iterator it(darts().begin()),
             itend(darts().end()); it!=itend; ++it)
      {
        Helper::template Foreach_enabled_attributes
          <internal::Correct_invalid_attributes_functor<Self> >::
          run(*this, it, marks);

        if ( vertex_attribute(it)==null_handle )
        {
          // If a dart don't have a 0-attribute, we create a Point at the origin
          set_vertex_attribute(it, create_vertex_attribute(CGAL::ORIGIN));
        }
      }

      for ( unsigned int i=0; i<=dimension; ++i)
        if ( marks[i]!=Base::INVALID_MARK )
        {
          CGAL_assertion( this->is_whole_map_marked(marks[i]) );
          free_mark(marks[i]);
        }

      Helper::template
        Foreach_enabled_attributes<internal::Cleanup_useless_attributes<Self> >::
          run(*this);
    }

    /** test if the two given facets have the same geometry
     * @return true iff the two facets have the same geometry.
     */
    bool are_facets_same_geometry(Dart_const_handle d1,
                                  Dart_const_handle d2) const
    {
      Dart_const_handle s1=d1;
      Dart_const_handle s2=d2;
      while (is_previous_exist(d1) && previous(s1)!=d1)
      {
        s1=previous(s1);
        if (!is_previous_exist(d2)) return false;
        s2=previous(s2);
      }

      d1=s1;
      d2=s2;
      do
      {
        if (is_next_exist(d1)!=is_next_exist(d2))
          return false;

        if (point(d1)!=point(d2))
          return false;

        d1=next(d1);
        d2=next(d2);
      }
      while(is_next_exist(d1) && d1!=s1);

      if (is_next_exist(d1)!=is_next_exist(d2))
          return false;

      if (d1==s1 && d2!=s2) return false;

      return true;
    }

    /** test if the two given facets have the same geometry but with
     *  opposite orientations.
     * @return true iff the two facets have the same geometry with opposite
     *         orientation.
     */
    bool are_facets_opposite_and_same_geometry(Dart_const_handle d1,
                                               Dart_const_handle d2) const
    {
      Dart_const_handle s1=d1;
      Dart_const_handle s2=d2;
      while (is_previous_exist(d1) && previous(s1)!=d1)
      {
        s1=previous(s1);
        if (!is_next_exist(d2)) return false;
        s2=next(s2);
      }

      d1=s1;
      d2=s2;
      do
      {
        if (is_next_exist(d1)!=is_previous_exist(d2))
          return false;

        if (other_extremity(d2)!=null_handle &&
             point(d1)!=point(other_extremity(d2)))
          return false;

        // The only case where d2 could have no other_extremity
        // is the end of an open path. In this case d1 is the
        // beginning of an open path and we do not compare points
        // but this is the correct thing to do.

        d1=next(d1);
        d2=previous(d2);
      }
      while(is_next_exist(d1) && d1!=s1);

      if (is_next_exist(d1)!=is_previous_exist(d2))
          return false;

      if (d1==s1 && d2!=s2) return false;

      return true;
    }

    /// Sew3 the marked facets having same geometry
    /// (a facet is considered marked if one of its dart is marked).
    unsigned int sew3_same_facets(size_type AMark)
    {
      unsigned int res = 0;

      std::map<Point, std::vector<Dart_handle> > one_dart_per_facet;
      size_type mymark = get_new_mark();

      // First we fill the std::map by one dart per facet, and by using
      // the minimal point as index.
      for (typename Dart_range::iterator it(darts().begin()),
             itend(darts().end()); it!=itend; ++it )
      {
        if ( !is_marked(it, mymark) && is_marked(it, AMark) )
        {
          Point min_point=point(it);
          Dart_handle min_dart = it;
          mark(it, mymark);
          typename Base::template
            Dart_of_cell_range<2>::iterator it2(*this,it);
          ++it2;
          for ( ; it2.cont(); ++it2 )
          {
            Point cur_point=point(it2);
            this->mark(it2, mymark);
            if ( cur_point < min_point )
            {
              min_point = cur_point;
              min_dart = it2;
            }
          }
          one_dart_per_facet[min_point].push_back(min_dart);
        }
        else
          this->mark(it, mymark);
      }

      // Second we run through the map: candidates for sew3 have necessary the
      // same minimal point.
      typename std::map<Point, std::vector<Dart_handle> >::iterator
        itmap=one_dart_per_facet.begin(),
        itmapend=one_dart_per_facet.end();

      for (; itmap!=itmapend; ++itmap)
      {
        for (typename std::vector<Dart_handle>::iterator
               it1=(itmap->second).begin(),
               it1end=(itmap->second).end(); it1!=it1end; ++it1)
        {
          typename std::vector<Dart_handle>::iterator it2=it1;
          for (++it2; it2!=it1end; ++it2)
          {
            if (*it1!=*it2 &&
                !this->template is_opposite_exist<3>(*it1) &&
                !this->template is_opposite_exist<3>(*it2) &&
                are_facets_opposite_and_same_geometry
                (*it1, this->previous(*it2)))
            {
              ++res;
              this->template sew<3>(*it1,
                                    this->other_orientation(previous(*it2)));
            }
          }
        }
      }

      CGAL_assertion( this->is_whole_map_marked(mymark) );
      this->free_mark(mymark);
      return res;
    }

    /// Sew3 the facets having same geometry
    /// (all the facets of the map are considered)
    unsigned int sew3_same_facets()
    {
      size_type mark = this->get_new_mark();
      this->negate_mark(mark);
      unsigned int res=sew3_same_facets(mark);
      this->free_mark(mark);
      return res;
    }

    /** Create a segment given 2 points.
     * @param p0 the first point.
     * @param p1 the second point.
     * if closed==true, the edge has no 2-free dart.
     * @return the dart of the new segment incident to p0.
     */
    Dart_handle make_segment(const Point& p0,const Point& p1,
                             bool closed=false)
    {
      return make_segment(create_vertex_attribute(p0),
                          create_vertex_attribute(p1),
                          closed);
    }

    /** Create a triangle given 3 points.
     * @param p0 the first point.
     * @param p1 the second point.
     * @param p2 the third point.
     * @return the dart of the new triangle incident to p0 and edge p0p1.
     */
    Dart_handle make_triangle(const Point& p0,
                              const Point& p1,
                              const Point& p2)
    {
      return make_triangle(create_vertex_attribute(p0),
                           create_vertex_attribute(p1),
                           create_vertex_attribute(p2));
    }

    /** Create a quadrangle given 4 points.
     * @param p0 the first point.
     * @param p1 the second point.
     * @param p2 the third point.
     * @param p3 the fourth point.
     * @return the dart of the new quadrangle incident to p0 and edge p0p1.
     */
    Dart_handle make_quadrangle(const Point& p0,
                                const Point& p1,
                                const Point& p2,
                                const Point& p3)
    {
      return make_quadrangle(create_vertex_attribute(p0),
                             create_vertex_attribute(p1),
                             create_vertex_attribute(p2),
                             create_vertex_attribute(p3));
    }


    /** Create a tetrahedron given 4 Vertex_attribute_handle.
     * @param h0 the first vertex handle.
     * @param h1 the second vertex handle.
     * @param h2 the third vertex handle.
     * @param h3 the fourth vertex handle.
     * @return the dart of the new tetrahedron incident to h0, to edge
     *         h0,h1 and to facet h0,h1,h2.
     */
    Dart_handle make_tetrahedron(Vertex_attribute_handle h0,
                                 Vertex_attribute_handle h1,
                                 Vertex_attribute_handle h2,
                                 Vertex_attribute_handle h3)
    {
      Dart_handle d1 = make_triangle(h0, h1, h2);
      Dart_handle d2 = make_triangle(h1, h0, h3);
      Dart_handle d3 = make_triangle(h1, h3, h2);
      Dart_handle d4 = make_triangle(h3, h0, h2);

      return this->make_combinatorial_tetrahedron(d1, d2, d3, d4);
    }

    /** Create a tetrahedron given 4 points.
     * @param p0 the first point.
     * @param p1 the second point.
     * @param p2 the third point.
     * @param p3 the fourth point.
     * @return the dart of the new tetrahedron incident to p0, to edge
     *         p0,p1 and to facet p0,p1,p2.
     */
    Dart_handle make_tetrahedron(const Point& p0,
                                 const Point& p1,
                                 const Point& p2,
                                 const Point& p3)
    {
      return make_tetrahedron(create_vertex_attribute(p0),
                              create_vertex_attribute(p1),
                              create_vertex_attribute(p2),
                              create_vertex_attribute(p3));
    }

    /** Create an hexahedron given 8 Vertex_attribute_handle.
     *    (8 vertices, 12 edges and 6 facets)
     * \verbatim
     *       4----7
     *      /|   /|
     *     5----6 |
     *     | 3--|-2
     *     |/   |/
     *     0----1
     * \endverbatim
     * @param h0 the first vertex handle.
     * @param h1 the second vertex handle.
     * @param h2 the third vertex handle.
     * @param h3 the fourth vertex handle.
     * @param h4 the fifth vertex handle.
     * @param h5 the sixth vertex handle.
     * @param h6 the seventh vertex handle.
     * @param h7 the height vertex handle.
     * @return the dart of the new hexahedron incident to h0, to edge
     *         h0,h5 and to the facet (h0,h5,h6,h1).
     */
    Dart_handle make_hexahedron(Vertex_attribute_handle h0,
                                Vertex_attribute_handle h1,
                                Vertex_attribute_handle h2,
                                Vertex_attribute_handle h3,
                                Vertex_attribute_handle h4,
                                Vertex_attribute_handle h5,
                                Vertex_attribute_handle h6,
                                Vertex_attribute_handle h7)
    {
      Dart_handle d1 = make_quadrangle(h0, h5, h6, h1);
      Dart_handle d2 = make_quadrangle(h1, h6, h7, h2);
      Dart_handle d3 = make_quadrangle(h2, h7, h4, h3);
      Dart_handle d4 = make_quadrangle(h3, h4, h5, h0);
      Dart_handle d5 = make_quadrangle(h0, h1, h2, h3);
      Dart_handle d6 = make_quadrangle(h5, h4, h7, h6);

      return this->make_combinatorial_hexahedron(d1, d2, d3, d4, d5, d6);
    }

    /** Create an hexahedron given 8 points.
     * \verbatim
     *       4----7
     *      /|   /|
     *     5----6 |
     *     | 3--|-2
     *     |/   |/
     *     0----1
     * \endverbatim
     * @param p0 the first point.
     * @param p1 the second point.
     * @param p2 the third point.
     * @param p3 the fourth point.
     * @param p4 the fifth point.
     * @param p5 the sixth point.
     * @param p6 the seventh point.
     * @param p7 the height point.
     * @return the dart of the new hexahedron incident to p0, to edge
     *         p0,p5 and to the facet (p0,p5,p6,p1).
     */
    Dart_handle make_hexahedron(const Point& p0,
                                const Point& p1,
                                const Point& p2,
                                const Point& p3,
                                const Point& p4,
                                const Point& p5,
                                const Point& p6,
                                const Point& p7)
    {
      return make_hexahedron(create_vertex_attribute(p0),
                             create_vertex_attribute(p1),
                             create_vertex_attribute(p2),
                             create_vertex_attribute(p3),
                             create_vertex_attribute(p4),
                             create_vertex_attribute(p5),
                             create_vertex_attribute(p6),
                             create_vertex_attribute(p7));
    }

    /** Compute the barycenter of a given cell.
     * @param adart a dart incident to the cell.
     * @param adim the dimension of the cell.
     * @return the barycenter of the cell.
     */
    template<unsigned int i>
    Point barycenter(Dart_const_handle adart) const
    {
      return CGAL::Barycenter_functor<Self, i>::run(*this, adart);
    }

    /** Insert a point in a given 1-cell.
     * @param dh a dart handle to the 1-cell
     * @param p the point to insert
     * @param update_attributes a boolean to update the enabled attributes
     * @return a dart handle to the new vertex containing p.
     */
    Dart_handle insert_point_in_cell_1(Dart_handle dh, const Point& p,
                                       bool update_attributes)
    {
      return this->insert_cell_0_in_cell_1(dh,
                                           create_vertex_attribute(p),
                                           update_attributes);
    }

    /** Insert a point in a given 2-cell.
     * @param dh a dart handle to the 2-cell
     * @param p the point to insert
     * @param update_attributes a boolean to update the enabled attributes
     * @return a dart handle to the new vertex containing p.
     */
    Dart_handle insert_point_in_cell_2(Dart_handle dh, const Point& p,
                                       bool update_attributes)
    {
      Vertex_attribute_handle v = create_vertex_attribute(p);

      Dart_handle first = this->insert_cell_0_in_cell_2(dh, v, update_attributes);

      if ( first==null_handle ) // If the triangulated facet was made of one dart
        erase_vertex_attribute(v);

#ifdef CGAL_CMAP_TEST_VALID_INSERTIONS
      CGAL_assertion( is_valid() );
#endif

      return first;
    }

    /** Insert a point in a given i-cell.
     * @param dh a dart handle to the i-cell
     * @param p the point to insert
     * @param update_attributes a boolean to update the enabled attributes
     * @return a dart handle to the new vertex containing p.
     */
    template <unsigned int i>
    Dart_handle insert_point_in_cell(Dart_handle dh, const Point& p,
                                     bool update_attributes = true)
    {
      CGAL_static_assertion(1<=i && i<=2);
      if (i==1) return insert_point_in_cell_1(dh, p, update_attributes);
      return insert_point_in_cell_2(dh, p, update_attributes);
    }

    /** Insert a dangling edge in a given facet.
     * @param dh a dart of the facet (!=NULL).
     * @param p the coordinates of the new vertex.
     * @param update_attributes a boolean to update the enabled attributes
     * @return a dart of the new edge, incident to the new vertex.
     */
    Dart_handle insert_dangling_cell_1_in_cell_2(Dart_handle dh,
                                                 const Point& p,
                                                 bool update_attributes = true)
    {
      return this->insert_dangling_cell_1_in_cell_2
          (dh, create_vertex_attribute(p), update_attributes);
    }

    /** Insert a point in a given i-cell.
     * @param dh a dart handle to the i-cell
     * @param p the point to insert
     * @param update_attributes a boolean to update the enabled attributes
     * @return a dart handle to the new vertex containing p.
     */
    template <unsigned int i>
    Dart_handle insert_barycenter_in_cell(Dart_handle dh, bool update_attributes = true)
    { return insert_point_in_cell<i>(dh, barycenter<i>(dh), update_attributes); }

    /** Compute the dual of a Linear_cell_complex.
     * @param alcc the lcc in which we build the dual of this lcc.
     * @param adart a dart of the initial lcc, NULL by default.
     * @return adart of the dual lcc, the dual of adart if adart!=NULL,
     *         any dart otherwise.
     * As soon as we don't modify this lcc and alcc lcc, we can iterate
     * simultaneously through all the darts of the two lcc and we have
     * each time of the iteration two "dual" darts.
     */
    Dart_handle dual_points_at_barycenter(Self & alcc, Dart_handle adart=null_handle)
    {
      Dart_handle res = Base::dual(alcc, adart);

      // Now the lcc alcc is topologically correct, we just need to add
      // its geometry to each vertex (the barycenter of the corresponding
      // dim-cell in the initial map).
      typename Dart_range::iterator it2 = alcc.darts().begin();
      for (typename Dart_range::iterator it(darts().begin());
           it!=darts().end(); ++it, ++it2)
      {
        if (vertex_attribute(it2)==null_handle)
        {
          alcc.set_vertex_attribute(it2, alcc.create_vertex_attribute
                                    (barycenter<dimension>(it)));
        }
      }

      return res;
    }

    /** Set the status of the managment of the attributes of the Map
     */
    void set_update_attributes(bool newval)
    {
      if (this->automatic_attributes_management == false && newval == true)
      {
        // We need to recode this function because correct_invalid_attributes
        // is not a virtual function.
        correct_invalid_attributes();
      }

      this->automatic_attributes_management = newval;
    }
  };

} // namespace CGAL

#endif // CGAL_LINEAR_CELL_COMPLEX_BASE_H //
// EOF //
