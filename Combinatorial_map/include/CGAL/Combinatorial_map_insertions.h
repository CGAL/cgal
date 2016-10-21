// Copyright (c) 2010-2011 CNRS and LIRIS' Establishments (France).
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
#ifndef CGAL_COMBINATORIAL_MAP_INSERTIONS_H
#define CGAL_COMBINATORIAL_MAP_INSERTIONS_H

#include <CGAL/config.h>

namespace CGAL
{
/** @file Combinatorial_map_insertions.h
 * Insertion operations on combinatorial map.
 */

#ifndef CGAL_NO_DEPRECATED_CODE

/** Insert a vertex in a given edge.
 * @param amap the used combinatorial map.
 * @param adart a dart of the edge (!=NULL && !=null_dart_handle).
 * @param update_attributes a boolean to update the enabled attributes
 *        (deprecated, now we use are_attributes_automatically_managed())
 * @return a dart of the new vertex.
 */
template<class CMap>
CGAL_DEPRECATED typename CMap::Dart_handle
insert_cell_0_in_cell_1( CMap& amap, typename CMap::Dart_handle adart,
                         typename CMap::template
                         Attribute_handle<0>::type ah=CMap::null_handle,
                         bool update_attributes=true )
{
  return amap.insert_cell_0_in_cell_1(adart, ah, update_attributes);
}

/** Insert a vertex in the given 2-cell which is splitted in triangles,
 * once for each inital edge of the facet.
 * @param amap the used combinatorial map.
 * @param adart a dart of the facet to triangulate.
 * @param update_attributes a boolean to update the enabled attributes
 *        (deprecated, now we use are_attributes_automatically_managed())
 * @return A dart incident to the new vertex.
 */
template < class CMap >
CGAL_DEPRECATED typename CMap::Dart_handle
insert_cell_0_in_cell_2( CMap& amap, typename CMap::Dart_handle adart,
                         typename CMap::template
                         Attribute_handle<0>::type ah=CMap::null_handle,
                         bool update_attributes=true )
{
  return amap.insert_cell_0_in_cell_2(adart, ah, update_attributes);
}
/** Insert a dangling edge in a 2-cell between given by a dart.
 * @param amap the used combinatorial map.
 * @param adart1 a first dart of the facet (!=NULL && !=null_dart_handle).
 * @param update_attributes a boolean to update the enabled attributes
 *        (deprecated, now we use are_attributes_automatically_managed())
 * @return a dart of the new edge, not incident to the vertex of adart1.
 */
template<class CMap>
CGAL_DEPRECATED typename CMap::Dart_handle
insert_dangling_cell_1_in_cell_2( CMap& amap,
                                  typename CMap::Dart_handle adart1,
                                  typename CMap::template
                                  Attribute_handle<0>::type ah=CMap::null_handle,
                                  bool update_attributes=true )
{
  return amap.insert_dangling_cell_1_in_cell_2(adart1,
                                               ah, update_attributes);
}

/** Test if an edge can be inserted onto a 2-cell between two given darts.
 * @param amap the used combinatorial map.
 * @param adart1 a first dart.
 * @param adart2 a second dart.
 * @return true iff an edge can be inserted between adart1 and adart2.
 */
template < class CMap >
CGAL_DEPRECATED bool is_insertable_cell_1_in_cell_2
(const CMap& amap, typename CMap::Dart_const_handle adart1,
 typename CMap::Dart_const_handle adart2)
{
  return amap.is_insertable_cell_1_in_cell_2(adart1, adart2);
}

/** Insert an edge in a 2-cell between two given darts.
 * @param amap the used combinatorial map.
 * @param adart1 a first dart of the facet (!=NULL && !=null_dart_handle).
 * @param adart2 a second dart of the facet. If NULL insert a dangling edge.
 * @param update_attributes a boolean to update the enabled attributes
 *        (deprecated, now we use are_attributes_automatically_managed())
 * @return a dart of the new edge, and not incident to the
 *         same vertex than adart1.
 */
template<class CMap>
CGAL_DEPRECATED typename CMap::Dart_handle
insert_cell_1_in_cell_2(CMap& amap,
                        typename CMap::Dart_handle adart1,
                        typename CMap::Dart_handle adart2,
                        bool update_attributes=true)
{
  return amap.insert_cell_1_in_cell_2(adart1, adart2,
                                      update_attributes);
}

/** Test if a 2-cell can be inserted onto a given 3-cell along
 * a path of edges.
 * @param amap the used combinatorial map.
 * @param afirst iterator on the begining of the path.
 * @param alast  iterator on the end of the path.
 * @return true iff a 2-cell can be inserted along the path.
 */
template <class CMap, class InputIterator>
CGAL_DEPRECATED bool is_insertable_cell_2_in_cell_3
(const CMap& amap, InputIterator afirst, InputIterator alast)
{
  return amap.is_insertable_cell_2_in_cell_3(afirst, alast);
}

/** Insert a 2-cell in a given 3-cell along a path of darts.
 * @param amap the used combinatorial map.
 * @param afirst iterator on the begining of the path.
 * @param alast  iterator on the end of the path.
 * @param update_attributes a boolean to update the enabled attributes
 *        (deprecated, now we use are_attributes_automatically_managed())
 * @return a dart of the new 2-cell.
 */
template<class CMap, class InputIterator>
CGAL_DEPRECATED typename CMap::Dart_handle
insert_cell_2_in_cell_3(CMap& amap, InputIterator afirst, InputIterator alast,
                        bool update_attributes=true)
{
  return amap.insert_cell_2_in_cell_3(afirst, alast, update_attributes);
}

#endif // CGAL_NO_DEPRECATED_CODE

} // namespace CGAL

#endif // CGAL_COMBINATORIAL_MAP_INSERTIONS_H
