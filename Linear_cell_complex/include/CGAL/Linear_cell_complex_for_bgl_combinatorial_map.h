// Copyright (c) 2017 CNRS and LIRIS' Establishments (France).
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
#ifndef CGAL_LINEAR_CELL_COMPLEX_FOR_BGL_COMBINATORIAL_MAP_H
#define CGAL_LINEAR_CELL_COMPLEX_FOR_BGL_COMBINATORIAL_MAP_H 1

#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex_traits.h>
#include <CGAL/CMap_linear_cell_complex_storages.h>
#include <CGAL/Cell_attribute_with_point.h>

namespace CGAL {

  /** @file Linear_cell_complex_for_bgl_combinatorial_map.h
   * Definition of a linear cell complex based on combinatorial map, to use with BGL.
   * This cmap has points associated to all vertices and faces enables. Moreover
   * these cells have id.
   */

  struct Linear_cell_complex_bgl_min_items
  {
    /// Dart_wrapper defines the type of darts used.
    template <class LCC>
    struct Dart_wrapper
    {
      typedef CGAL::Tag_true Darts_with_id;
      typedef CGAL::Cell_attribute_with_point_and_id<LCC> Vertex_attribute;
      typedef CGAL::Cell_attribute_with_id<LCC> Face_attribute;
      typedef CGAL::cpp11::tuple<Vertex_attribute, void, Face_attribute> Attributes;
    };
  };
  
  // Linear_cell_complex_for_bgl_combinatorial_map class.
  template < unsigned int d_, unsigned int ambient_dim = d_,
             class Traits_ = Linear_cell_complex_traits<ambient_dim>,
             class Items_ = Linear_cell_complex_bgl_min_items,
             class Alloc_ = CGAL_ALLOCATOR(int),
             template<unsigned int,class,class,class,class>
             class CMap = Combinatorial_map_base,
             class Storage_ = CMap_linear_cell_complex_storage_1<d_, ambient_dim,
                                                                 Traits_, Items_,
                                                                 Alloc_> >
    class Linear_cell_complex_for_bgl_combinatorial_map:
        public Linear_cell_complex_for_combinatorial_map
                                       <d_, ambient_dim, Traits_,
                                        Items_, Alloc_, CMap, Storage_>
    {
    public:
      typedef Linear_cell_complex_for_bgl_combinatorial_map<d_, ambient_dim,
                          Traits_, Items_, Alloc_, CMap, Storage_>  Self;

      typedef Linear_cell_complex_for_combinatorial_map<d_, ambient_dim,
                          Traits_, Items_, Alloc_, CMap, Storage_> Base;

      typedef Traits_ Traits;
      typedef Items_  Items;
      typedef Alloc_  Alloc;

      static const unsigned int ambient_dimension = Base::ambient_dimension;
      static const unsigned int dimension = Base::dimension;

      typedef typename Base::Dart_handle       Dart_handle;
      typedef typename Base::Dart_const_handle Dart_const_handle;
      typedef typename Base::Helper            Helper;

      typedef typename Base::Point  Point;
      typedef typename Base::Vector Vector;
      typedef typename Base::FT     FT;

      typedef typename Base::Dart_range Dart_range;

      typedef typename Base::template Attribute_type<0>::type Vertex_attribute;
      typedef typename Base::template Attribute_handle<0>::type
      Vertex_attribute_handle;
      typedef typename Base::template Attribute_const_handle<0>::type
      Vertex_attribute_const_handle;

      typedef typename Base::template Attribute_range<0>::type
      Vertex_attribute_range;
      typedef typename Base::template Attribute_const_range<0>::type
      Vertex_attribute_const_range;

      typedef typename Base::template Attribute_type<2>::type Face_attribute;
      typedef typename Base::template Attribute_handle<2>::type
      Face_attribute_handle;
      typedef typename Base::template Attribute_const_handle<2>::type
      Face_attribute_const_handle;

      typedef typename Base::template Attribute_range<2>::type
      Face_attribute_range;
      typedef typename Base::template Attribute_const_range<2>::type
      Face_attribute_const_range;

      typedef typename Base::size_type size_type;

      typedef typename Base::Use_index Use_index;
      typedef typename Base::Storage Storage;
      typedef typename Base::Exception_no_more_available_mark
      Exception_no_more_available_mark;

      Linear_cell_complex_for_bgl_combinatorial_map() : Base()
      {}

      /** Copy the given linear cell complex into *this.
       *  Note that both LCC can have different dimensions and/or non void attributes.
       *  @param alcc the linear cell complex to copy.
       *  @post *this is valid.
       */
      Linear_cell_complex_for_bgl_combinatorial_map(const Self& alcc) : Base(alcc)
      {}

      template < class LCC2 >
      Linear_cell_complex_for_bgl_combinatorial_map(const LCC2& alcc) : Base(alcc)
      {}      

      template < class LCC2, typename Converters >
      Linear_cell_complex_for_bgl_combinatorial_map(const LCC2& alcc,
                                                    Converters& converters) :
        Base(alcc, converters)
      {}

      template < class LCC2, typename Converters, typename DartInfoConverter >
      Linear_cell_complex_for_bgl_combinatorial_map(const LCC2& alcc,
                                                    Converters& converters,
                                                    const DartInfoConverter&
                                                    dartinfoconverter) :
        Base(alcc, converters, dartinfoconverter)
      {}

      template < class LCC2, typename Converters, typename DartInfoConverter,
                 typename PointConverter >
      Linear_cell_complex_for_bgl_combinatorial_map(const LCC2& alcc,
                                                    Converters& converters,
                                                    const DartInfoConverter&
                                                    dartinfoconverter,
                                                    const PointConverter&
                                                    pointconverter) :
        Base(alcc, converters, dartinfoconverter, pointconverter)
      {}

    };

} // namespace CGAL

#endif // CGAL_LINEAR_CELL_COMPLEX_FOR_BGL_COMBINATORIAL_MAP_H //
// EOF //
