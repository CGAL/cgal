// Copyright (c) 2005, 2006 Fernando Luis Cacciola Carballal. All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: $
// $Id: $
//
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>
//
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_LINDSTROMTURK_SELECTION_MAP_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_LINDSTROMTURK_SELECTION_MAP_H 1

#include <CGAL/Surface_mesh_simplification/TSMS_common.h>


CGAL_BEGIN_NAMESPACE

namespace Triangulated_surface_mesh { namespace Simplification 
{

//*******************************************************************************************************************
//                               -= vertex-pair selection property maps =-
//
// Determines up-front whether a pair must be considered for collapsation
//
//*******************************************************************************************************************


//
// Lindstrom-Turk selection: only vertex-pairs which are edges of the triangulated surface mesh
//
template<class TSM_>
class Lindstrom_Turk_selection : public boost::put_get_helper< bool, Lindstrom_Turk_selection<TSM_> >
{
public:
    
    typedef TSM_ TSM ;
    
    typedef typename boost::graph_traits<TSM>::vertex_descriptor vertex_descriptor ;
     
public:

    typedef boost::readable_property_map_tag category;
    
    typedef bool value_type;
    
    typedef value_type reference ;
    
    typedef boost::tuple<vertex_descriptor,vertex_descriptor,bool,TSM*> key_type;

    value_type operator[](key_type const& e) const
    {
      bool lIsEdge ;
      boost::tie(boost::tuples::ignore,boost::tuples::ignore,lIsEdge,boost::tuples::ignore) = e ;
      return lIsEdge ;
    }
};

} } // namespace Triangulated_surface_mesh::Simplification

CGAL_END_NAMESPACE

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_LINDSTROMTURK_SELECTION_MAP_H //
// EOF //
 
