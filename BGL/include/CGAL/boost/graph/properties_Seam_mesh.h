// Copyright (c) 2015  GeometryFactory (France).  All rights reserved.
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
// Author(s)     : Andreas Fabri


#ifndef CGAL_PROPERTIES_SEAM_MESH_H
#define CGAL_PROPERTIES_SEAM_MESH_H

#include <CGAL/boost/graph/Seam_mesh.h>

#include <CGAL/boost/graph/properties.h>

namespace CGAL {

template <class TM>
class Seam_mesh;

} // namespace CGAL


namespace boost {

template<typename P>
struct property_map<CGAL::Seam_mesh<P>, CGAL::vertex_point_t >
{
  typedef CGAL::Seam_mesh<P> SM; 

  typedef typename property_map<P,CGAL::vertex_point_t >::type type;
  
  typedef type const_type;
  
};

template<typename K>
typename property_map<CGAL::Seam_mesh<K>, CGAL::vertex_point_t >::const_type
get(CGAL::vertex_point_t, const CGAL::Seam_mesh<K>& g) {
  return get(CGAL::vertex_point_t(), const_cast<K&>(g.mesh()));
}


} // namespace boost 

#endif // CGAL_PROPERTIES_SEAM_MESH_H
