// Copyright (c) 2011 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
#ifndef CGAL_LINEAR_CELL_COMPLEX_FOR_GENERALIZED_MAP_H
#define CGAL_LINEAR_CELL_COMPLEX_FOR_GENERALIZED_MAP_H 1

#include <CGAL/Linear_cell_complex_fwd.h>
#include <CGAL/Linear_cell_complex_base.h>
#include <CGAL/Linear_cell_complex_traits.h>
#include <CGAL/Linear_cell_complex_min_items.h>
#include <CGAL/Generalized_map.h>
#include <CGAL/GMap_linear_cell_complex_storages.h>

namespace CGAL {

  /** @file Linear_cell_complex_for_generalized_map.h
   * Definition of a linear cell complex based on generalized map, having
   * points associated to all vertices.
   */

  // Linear_cell_complex_for_generalized_map class.
  // No difference with class Linear_cell_complex_base except the default
  // template parameters for Refs class which is a generalized map.
  template < unsigned int d_, unsigned int ambient_dim,
             class Traits_, class Items_, class Alloc_,
             template<unsigned int,class,class,class,class> class CMap,
             class Storage_ >
    class Linear_cell_complex_for_generalized_map:
        public Linear_cell_complex_base<d_, ambient_dim, Traits_,
                                        Items_, Alloc_, CMap,
                                        Linear_cell_complex_for_generalized_map
                                        <d_, ambient_dim,
                                         Traits_, Items_,
                                         Alloc_, CMap, Storage_>,
                                        Storage_>
    {
    public:
      typedef Linear_cell_complex_for_generalized_map<d_, ambient_dim,
                          Traits_, Items_, Alloc_, CMap, Storage_>  Self;

      typedef Linear_cell_complex_base<d_, ambient_dim,
                          Traits_, Items_, Alloc_, CMap, Self, Storage_> Base;

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

      typedef typename Base::size_type size_type;

      typedef typename Base::Use_index Use_index;
      typedef typename Base::Storage Storage;
      typedef typename Base::Exception_no_more_available_mark
      Exception_no_more_available_mark;

      Linear_cell_complex_for_generalized_map() : Base()
      {}

      /** Copy the given linear cell complex into *this.
       *  Note that both LCC can have different dimensions and/or non void attributes.
       *  @param alcc the linear cell complex to copy.
       *  @post *this is valid.
       */
      Linear_cell_complex_for_generalized_map(const Self & alcc) : Base(alcc)
      {}

      template <unsigned int d2,  unsigned int ambient_dim2, class Traits2,
                class Items2, class Alloc2,
                template<unsigned int,class,class,class,class> class CMap2,
                class Storage2>
      Linear_cell_complex_for_generalized_map
      (const Linear_cell_complex_for_generalized_map<d2, ambient_dim2,
       Traits2, Items2, Alloc2, CMap2, Storage2>& alcc) : Base(alcc)
      {}

      template <unsigned int d2,  unsigned int ambient_dim2, class Traits2,
                class Items2, class Alloc2,
                template<unsigned int,class,class,class,class> class CMap2,
                class Storage2, typename Converters>
      Linear_cell_complex_for_generalized_map
      (const Linear_cell_complex_for_generalized_map<d2, ambient_dim2,
       Traits2, Items2, Alloc2, CMap2, Storage2>& alcc,
       const Converters& converters) : Base(alcc, converters)
      {}

      template <unsigned int d2,  unsigned int ambient_dim2, class Traits2,
                class Items2, class Alloc2,
                template<unsigned int,class,class,class,class> class CMap2,
                class Storage2, typename Converters, typename DartInfoConverter>
      Linear_cell_complex_for_generalized_map
      (const Linear_cell_complex_for_generalized_map<d2, ambient_dim2,
       Traits2, Items2, Alloc2, CMap2, Storage2>& alcc,
       const Converters& converters,
       const DartInfoConverter& dartinfoconverter) :
        Base(alcc, converters, dartinfoconverter)
      {}

      template <unsigned int d2,  unsigned int ambient_dim2, class Traits2,
                class Items2, class Alloc2,
                template<unsigned int,class,class,class,class> class CMap2,
                class Storage2, typename Converters, typename DartInfoConverter,
                typename PointConverter>
      Linear_cell_complex_for_generalized_map
      (const Linear_cell_complex_for_generalized_map<d2, ambient_dim2,
       Traits2, Items2, Alloc2, CMap2, Storage2>& alcc,
       const Converters& converters,
       const DartInfoConverter& dartinfoconverter,
       const PointConverter& pointconverter) :
        Base(alcc, converters, dartinfoconverter, pointconverter)
      {}

      Self & operator= (const Self & alcc)
      {
        Base::operator=(alcc);
        return *this;
      }
    };

} // namespace CGAL

#endif // CGAL_LINEAR_CELL_COMPLEX_FOR_GENERALIZED_MAP_H //
// EOF //
