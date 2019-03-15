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
#ifndef CGAL_LINEAR_CELL_COMPLEX_FOR_COMBINATORIAL_MAP_H
#define CGAL_LINEAR_CELL_COMPLEX_FOR_COMBINATORIAL_MAP_H 1

#include <CGAL/Linear_cell_complex_base.h>
#include <CGAL/Linear_cell_complex_traits.h>
#include <CGAL/Linear_cell_complex_min_items.h>
#include <CGAL/Combinatorial_map.h>
#include <CGAL/CMap_linear_cell_complex_storages.h>
#include <CGAL/boost/graph/properties.h>

namespace CGAL {

  /** @file Linear_cell_complex_for_combinatorial_map.h
   * Definition of a linear cell complex based on combinatorial map, having
   * points associated to all vertices.
   */

  // Linear_cell_complex_for_combinatorial_map class.
  // No difference with class Linear_cell_complex_base except the default
  // template parameters for Refs class which is a combinatorial map.
  template < unsigned int d_, unsigned int ambient_dim = d_,
             class Traits_ = Linear_cell_complex_traits<ambient_dim>,
#if defined(CGAL_CMAP_DART_DEPRECATED) && !defined(CGAL_NO_DEPRECATED_CODE)
             class Items_ = Linear_cell_complex_min_items<d_>,
#else
             class Items_ = Linear_cell_complex_min_items,
#endif
             class Alloc_ = CGAL_ALLOCATOR(int),
             template<unsigned int,class,class,class,class>
             class CMap = Combinatorial_map_base,
             class Storage_ = CMap_linear_cell_complex_storage_1<d_, ambient_dim,
                                                                 Traits_, Items_,
                                                                 Alloc_> >
    class Linear_cell_complex_for_combinatorial_map:
        public Linear_cell_complex_base<d_, ambient_dim, Traits_,
                                        Items_, Alloc_, CMap,
                                        Linear_cell_complex_for_combinatorial_map
                                        <d_, ambient_dim,
                                         Traits_, Items_,
                                         Alloc_, CMap, Storage_>,
                                        Storage_>
    {
    public:
      typedef Linear_cell_complex_for_combinatorial_map<d_, ambient_dim,
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

      Linear_cell_complex_for_combinatorial_map() : Base()
      {}

      /** Copy the given linear cell complex into *this.
       *  Note that both LCC can have different dimensions and/or non void attributes.
       *  @param alcc the linear cell complex to copy.
       *  @post *this is valid.
       */
      Linear_cell_complex_for_combinatorial_map(const Self& alcc) : Base(alcc)
      {}

      template <unsigned int d2,  unsigned int ambient_dim2, class Traits2,
                class Items2, class Alloc2,
                template<unsigned int,class,class,class,class> class CMap2,
                class Storage2>
      Linear_cell_complex_for_combinatorial_map
      (const Linear_cell_complex_for_combinatorial_map<d2, ambient_dim2,
       Traits2, Items2, Alloc2, CMap2, Storage2>& alcc) : Base(alcc)
      {}      

      template <unsigned int d2,  unsigned int ambient_dim2, class Traits2,
                class Items2, class Alloc2,
                template<unsigned int,class,class,class,class> class CMap2,
                class Storage2, typename Converters>
      Linear_cell_complex_for_combinatorial_map
      (const Linear_cell_complex_for_combinatorial_map<d2, ambient_dim2,
       Traits2, Items2, Alloc2, CMap2, Storage2>& alcc,
       const Converters& converters) : Base(alcc, converters)
      {}

      template <unsigned int d2,  unsigned int ambient_dim2, class Traits2,
                class Items2, class Alloc2,
                template<unsigned int,class,class,class,class> class CMap2,
                class Storage2, typename Converters, typename DartInfoConverter>
      Linear_cell_complex_for_combinatorial_map
      (const Linear_cell_complex_for_combinatorial_map<d2, ambient_dim2,
       Traits2, Items2, Alloc2, CMap2, Storage2>& alcc,
       const Converters& converters,
       const DartInfoConverter& dartinfoconverter) :
        Base(alcc, converters, dartinfoconverter)
      {}

      template <unsigned int d2,  unsigned int ambient_dim2, class Traits2,
                class Items2, class Alloc2,
                template<unsigned int,class,class,class,class> class CMap2,
                class Storage2, typename Converters,
                typename DartInfoConverter, typename PointConverter>
      Linear_cell_complex_for_combinatorial_map
      (const Linear_cell_complex_for_combinatorial_map<d2, ambient_dim2,
       Traits2, Items2, Alloc2, CMap2, Storage2>& alcc,
       const Converters& converters, const DartInfoConverter& dartinfoconverter,
       const PointConverter& pointconverter) :
        Base(alcc, converters, dartinfoconverter, pointconverter)
      {}

      /** Import the given hds which should be a model of an halfedge graph. */
      template<class HEG, class PointConverter>
      void import_from_halfedge_graph(const HEG& heg              ,
                                      const PointConverter& pointconverter,
                                      boost::unordered_map
                                      <typename boost::graph_traits<HEG>::halfedge_descriptor,
                                      Dart_handle>* origin_to_copy=NULL,
                                      boost::unordered_map
                                      <Dart_handle,
                                      typename boost::graph_traits<HEG>::halfedge_descriptor>*
                                      copy_to_origin=NULL)

      {
        boost::unordered_map
            <typename boost::graph_traits<HEG>::halfedge_descriptor,
            Dart_handle> local_dartmap;
        if (origin_to_copy==NULL) // Used local_dartmap if user does not provides its own unordered_map
        { origin_to_copy=&local_dartmap; }

        Base::import_from_halfedge_graph(heg, origin_to_copy, copy_to_origin);

        typedef typename boost::property_map<HEG,vertex_point_t>::const_type
            Point_property_map;
        Point_property_map ppmap = get(CGAL::vertex_point, heg);

        typename boost::unordered_map
          <typename boost::graph_traits<HEG>::halfedge_descriptor,
           Dart_handle>::iterator dartmap_iter, dartmap_iter_end=origin_to_copy->end();
        for (dartmap_iter=origin_to_copy->begin(); dartmap_iter!=dartmap_iter_end;
             ++dartmap_iter)
        {
          if (this->vertex_attribute(dartmap_iter->second)==NULL)
          {
            this->set_vertex_attribute(dartmap_iter->second,
                                 this->create_vertex_attribute());
            pointconverter.run(ppmap[source(dartmap_iter->first, heg)],
                this->point(dartmap_iter->second));
          }
        }
      }

      /** Import the given hds which should be a model of an halfedge graph. */
      template<class HEG>
      void import_from_halfedge_graph(const HEG& heg,
                                      boost::unordered_map
                                      <typename boost::graph_traits<HEG>::halfedge_descriptor,
                                      Dart_handle>* origin_to_copy=NULL,
                                      boost::unordered_map
                                      <Dart_handle,
                                      typename boost::graph_traits<HEG>::halfedge_descriptor>*
                                      copy_to_origin=NULL)
      {
         typedef typename boost::property_traits<typename boost::property_map
            <HEG, vertex_point_t>::type>::value_type HEG_point;

        CGAL::internal::Set_point_if_possible_cmap<HEG_point, Point> default_point_converter;
        import_from_halfedge_graph(heg, default_point_converter,
                                   origin_to_copy, copy_to_origin);
      }

    };

} // namespace CGAL

#endif // CGAL_LINEAR_CELL_COMPLEX_FOR_COMBINATORIAL_MAP_H //
// EOF //
