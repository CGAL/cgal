// Copyright (c) 2010 INRIA Sophia-Antipolis (France).
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
// $URL: svn+ssh://mbogdanov@scm.gforge.inria.fr/svn/cgal/trunk/Mesh_3/include/CGAL/refine_mesh_3.h $
// $Id: refine_mesh_3.h 54197 2010-02-16 15:55:09Z stayeb $
//
//
// Author(s)     : Mikhail Bogdanov
//
//*****************************************************************************************
// File Description : refine_periodic_mesh_3 function declaration and implementation.
//*****************************************************************************************

#ifndef CGAL_REFINE_PERIODIC_MESH_3_H
#define CGAL_REFINE_PERIODIC_MESH_3_H

#include <CGAL/refine_mesh_3.h>

namespace CGAL {
  

// to be deleted!
template<class C3T3>
void dump_dimension(C3T3& c3t3)
{
  typedef typename C3T3::Triangulation::Finite_vertices_iterator Finite_vertices_iterator;
  
  int total_on_surface = 0;
  int total_default = 0;
  int total_volume = 0;
  int other = 0;
  
  for(Finite_vertices_iterator vit = c3t3.triangulation().finite_vertices_begin(); 
      vit != c3t3.triangulation().finite_vertices_end(); vit++) {
    //std::cout << vit->in_dimension() << std::endl;
    
    if(vit->in_dimension() == -1) {
      total_default += 1;
    }
    
    if(vit->in_dimension() == 2) {
      total_on_surface += 1;
    }
    
    if(vit->in_dimension() == 3) {
      total_volume += 1;
    }
    
    if(vit->in_dimension() != -1 && vit->in_dimension() != 2 && vit->in_dimension() != 3) {
      other++;
    }
  }
  
  std::cout << "Total number on surface: " << total_on_surface << std::endl; 
  std::cout << "Total number with -1: " << total_default << std::endl; 
  std::cout << "Total number of volume vertices: " << total_volume << std::endl;
  std::cout << "other " << other << std::endl;
}  

// to be deleted!
template<class C3T3>
void mark_vertices_in_complex(C3T3& c3t3)
{
  typedef typename C3T3::Cell_iterator Iterator;
  
  int number_of_changes = 0;
  for(Iterator cit = c3t3.cells_begin(); cit != c3t3.cells_end(); cit++) {
    
    for(int i = 0; i < 3; i++) {
      
      if(cit->vertex(i)->in_dimension() == -1) {
        cit->vertex(i)->set_dimension(3);
        number_of_changes += 1;
      }
      
    }
  }
  
  std::cout << "number of marked vertices: " << number_of_changes << std::endl;
}
    
/**
 * @brief This function refines the mesh c3t3 wrt domain & criteria
 *
 * @param c3t3 the mesh to be refined.
 * @param domain the domain to be discretized
 * @param criteria the criteria
 * @param reset_c3t3 if \c true, a new C3T3 will be construct from param c3t3.
 *   The new c3t3 keeps only the vertices (as NON-weighted points with their
 *   dimension and Index) of the triangulation. That allows to refine a mesh
 *   which has been exuded.
 */
template<class C3T3, class MeshDomain, class MeshCriteria>
void refine_periodic_mesh_3(C3T3& c3t3,
								 const MeshDomain&   domain,
								 const MeshCriteria& criteria,
                                 bool reset_c3t3)
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
  
  dump_dimension(c3t3);
  
  // Build mesher and launch refinement process
  Mesher mesher (c3t3, domain, criteria);
  
  dump_dimension(c3t3);
  
  double refine_time = mesher.refine_mesh();
  
  dump_dimension(c3t3);
  
  mark_vertices_in_complex(c3t3);
  
  dump_dimension(c3t3);
  
  std::cout << "the end" << std::endl;
  
  // Odt
  if ( /*odt*/false )
  {
    odt_optimize_mesh_3(c3t3,
                        domain,
                        0,//parameters::time_limit = odt.time_limit(),
                        0,//parameters::max_iteration_number = odt.max_iteration_number(),
                        0.02,//parameters::convergence = odt.convergence(),
                        0.01//parameters::freeze_bound = odt.bound());
                        );
  }
  
  // Lloyd
  if ( /*lloyd*/false )
  {
    lloyd_optimize_mesh_3(c3t3,
                          domain,
                          0,//parameters::time_limit = lloyd.time_limit(),
                          0,//parameters::max_iteration_number = lloyd.max_iteration_number(),
                          0.02,//parameters::convergence = lloyd.convergence(),
                          0.01//parameters::freeze_bound = lloyd.bound());
                          );
  }
  
  // Perturbation
  if ( /*perturb*/false )
  {
    double perturb_time_limit = refine_time;
    
    //if ( perturb.is_time_limit_set() )
      //perturb_time_limit = perturb.time_limit();
    
    perturb_mesh_3(c3t3,
                   domain,
                   20,//parameters::time_limit = perturb_time_limit,
                   0//parameters::sliver_bound = perturb.bound());
                   );
  }
  
  dump_dimension(c3t3);
  
  // Exudation
  //if ( /*exude*/false )
  /*{
    double exude_time_limit = refine_time;
    
    //if ( exude.is_time_limit_set() )
      //exude_time_limit = exude.time_limit();
    
    exude_mesh_3(c3t3,
                 0,//parameters::time_limit = exude_time_limit,
                 0//parameters::sliver_bound = exude.bound());
                 );
  }*/
};

}  // end namespace CGAL


#endif // CGAL_REFINE_PERIODIC_MESH_3_H