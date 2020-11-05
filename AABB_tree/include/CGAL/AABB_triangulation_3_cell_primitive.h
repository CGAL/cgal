// Copyright (c) 2017 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author        : Jane Tournois
//

#ifndef CGAL_AABB_TRIANGULATION_3_CELL_PRIMITIVE_H_
#define CGAL_AABB_TRIANGULATION_3_CELL_PRIMITIVE_H_

#include <CGAL/license/AABB_tree.h>


#include <CGAL/AABB_primitive.h>
#include <CGAL/result_of.h>
#include <iterator>

namespace CGAL
{
  namespace internal
  {
    template <class GeomTraits, class Iterator>
    struct Point_from_cell_iterator_proprety_map
    {
      //classical typedefs
      typedef Iterator key_type;
      typedef typename GeomTraits::Point_3 value_type;
      typedef typename cpp11::result_of<
        typename GeomTraits::Construct_vertex_3(typename GeomTraits::Tetrahedron_3, int)
      >::type reference;
      typedef boost::readable_property_map_tag category;

      inline friend
        typename Point_from_cell_iterator_proprety_map<GeomTraits, Iterator>::reference
        get(Point_from_cell_iterator_proprety_map<GeomTraits, Iterator>, Iterator it)
      {
        typename GeomTraits::Construct_point_3 point;
        return point(it->vertex(1)->point());
      }
    };

    template <class GeomTraits, class Iterator>
    struct Tet_from_cell_iterator_proprety_map
    {
      //classical typedefs
      typedef Iterator                           key_type;
      typedef typename GeomTraits::Tetrahedron_3 value_type;
      typedef value_type                         reference;
      typedef boost::readable_property_map_tag category;

      inline friend
        reference
        get(Tet_from_cell_iterator_proprety_map<GeomTraits, Iterator>, key_type it)
      {
        typename GeomTraits::Construct_point_3 point;
        return value_type(point(it->vertex(0)->point()),
                          point(it->vertex(1)->point()),
                          point(it->vertex(2)->point()),
                          point(it->vertex(3)->point()));
      }
    };

  }//namespace internal


  template < class GeomTraits,
             class Tr,
             class CacheDatum = Tag_false,
             class Handle = typename Tr::Cell_handle>
  class AABB_triangulation_3_cell_primitive
#ifndef DOXYGEN_RUNNING
    : public AABB_primitive<  Handle,
          internal::Tet_from_cell_iterator_proprety_map<GeomTraits, Handle>,
          internal::Point_from_cell_iterator_proprety_map<GeomTraits, Handle>,
          Tag_false,
          CacheDatum >
#endif
  {
    typedef AABB_primitive< Handle,
          internal::Tet_from_cell_iterator_proprety_map<GeomTraits, Handle>,
          internal::Point_from_cell_iterator_proprety_map<GeomTraits, Handle>,
          Tag_false,
          CacheDatum > Base;
  public:
    AABB_triangulation_3_cell_primitive(Handle h) : Base(h){}  };

}  // end namespace CGAL


#endif // CGAL_AABB_TRIANGULATION_3_CELL_PRIMITIVE_H_
