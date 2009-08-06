// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// All rights reserved.
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
// $URL$
// $Id$
//
//
// Author(s)     : Stephane Tayeb
//
//******************************************************************************
// File Description : refine_mesh_3 function declaration and implementation.
//******************************************************************************

#ifndef CGAL_REFINE_MESH_3_H
#define CGAL_REFINE_MESH_3_H

#include <CGAL/Mesh_3/Slivers_exuder.h>
#include <CGAL/Mesh_3/Mesher_3.h>

namespace CGAL {

  namespace details {
    
    /**
     * @class Insert_vertex_in_c3t3
     *
     * A functor designed to insert unweighted points into the triangulation
     * of a C3T3 from C3T3::Tr::Vertex , keeping the dimension and indices.
     */
    template <typename C3T3>
    class Insert_vertex_in_c3t3
    {
    private:
      typedef typename C3T3::Vertex_handle          Vertex_handle;
      typedef typename C3T3::Index                  Index;
      
      typedef typename C3T3::Triangulation          Tr;
      typedef typename Tr::Vertex                   Vertex;
      typedef typename Tr::Bare_point               Bare_point;
      
    public:
      Insert_vertex_in_c3t3(C3T3& c3t3)
        : r_c3t3_(c3t3) { };
      
      void operator()(const Vertex& vertex) const
      {
        // Get vh properties
        Bare_point point = vertex.point();
        int dimension = vertex.in_dimension();
        Index index = vertex.index();
        
        // Insert point and restore handle properties
        Vertex_handle new_vertex = r_c3t3_.triangulation().insert(point);
        r_c3t3_.set_index(new_vertex, index);
        r_c3t3_.set_dimension(new_vertex, dimension);
      }

    private:
      C3T3& r_c3t3_;
    };
  }
  
  
/**
 * @brief This function refines the mesh c3t3 wrt domain & criteria
 *
 * @param c3t3 the mesh to be refined.
 * @param domain the domain to be discretized
 * @param criteria the criteria
 * @param exude if \c true, an exudation step will be done at
 *   the end of the Delaunay refinment process
 * @param reset_c3t3 if \c true, a new C3T3 will be construct from param c3t3.
 *   The new c3t3 keeps only the vertices (as NON-weighted points with their
 *   dimension and Index) of the triangulation. That allows to refine a mesh
 *   which has been exuded.
 */
template<class C3T3, class MeshDomain, class MeshCriteria>
void refine_mesh_3(C3T3& c3t3,
                   const MeshDomain&   domain,
                   const MeshCriteria& criteria,
                   bool exude = true,
                   bool reset_c3t3 = true)
{
  typedef Mesh_3::Mesher_3<C3T3, MeshCriteria, MeshDomain> Mesher;
  typedef typename C3T3::Triangulation::Geom_traits Gt;
  
  // Reset c3t3 (i.e. remove weights) if needed
  if ( reset_c3t3 )
  {
    C3T3 tmp_c3t3;
    std::for_each(c3t3.triangulation().finite_vertices_begin(),
                  c3t3.triangulation().finite_vertices_end(),
                  details::Insert_vertex_in_c3t3<C3T3>(tmp_c3t3));
    c3t3.swap(tmp_c3t3);
  }
  
  // Build mesher and launch refinement process
  Mesher mesher (c3t3, domain, criteria);
  mesher.refine_mesh();

  // Exudation
  if ( exude )
  {
    typedef Mesh_3::Min_dihedral_angle_criterion<Gt> Sliver_criterion;
    //typedef Mesh_3::Radius_radio_criterion<Gt> Sliver_criterion;
    typedef typename Mesh_3::Slivers_exuder<C3T3, Sliver_criterion> Exuder;
    
    Exuder exuder(c3t3);
  
#ifdef CGAL_MESH_3_VERBOSE
    exuder.print_stats();
#endif // CGAL_MESH_3_VERBOSE
  
    exuder.pump_vertices();
  
#ifdef CGAL_MESH_3_VERBOSE
    exuder.print_stats();
#endif // CGAL_MESH_3_VERBOSE
  }
  
};

}  // end namespace CGAL


#endif // CGAL_REFINE_MESH_3_H
