// Copyright (c) 2009, 2014 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
// File Description : Test meshing utilities.
//******************************************************************************

#ifndef CGAL_MESH_3_TEST_TEST_MESHING_UTILITIES
#define CGAL_MESH_3_TEST_TEST_MESHING_UTILITIES

#include "test_utilities.h"

#include <CGAL/Mesh_complex_3_in_triangulation_3.h>

#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/refine_mesh_3.h>
#include <CGAL/optimize_mesh_3.h>
#include <CGAL/remove_far_points_in_mesh_3.h>

#include <CGAL/Mesh_3/Triangle_accessor_primitive.h>
#include <CGAL/Triangle_accessor_3.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>

#include <vector>
#include <boost/optional/optional_io.hpp>

// IO
#include <fstream>
#include <iostream>

#include <climits>
#define STD_SIZE_T_MAX UINT_MAX

struct Bissection_tag {};
struct Polyhedral_tag {};


template <typename K>
struct Tester
{
  template<typename C3t3,
           typename Domain,
           typename Criteria,
           typename Domain_type_tag
           >
  void verify(C3t3& c3t3,
              const Domain& domain,
              const Criteria& criteria,
              const Domain_type_tag domain_type,
              const std::size_t min_vertices_expected = 0,
              const std::size_t max_vertices_expected = STD_SIZE_T_MAX,
              const std::size_t min_facets_expected = 0,
              const std::size_t max_facets_expected = STD_SIZE_T_MAX,
              const std::size_t min_cells_expected = 0,
              const std::size_t max_cells_expected = STD_SIZE_T_MAX) const
  {
    typedef typename C3t3::size_type size_type;

    // Store mesh properties
    size_type v = c3t3.triangulation().number_of_vertices();
    size_type f = c3t3.number_of_facets_in_complex();
    size_type c = c3t3.number_of_cells_in_complex();

    // Verify
    verify_c3t3(c3t3,domain,domain_type,
                min_vertices_expected,
                max_vertices_expected,
                min_facets_expected,
                max_facets_expected,
                min_cells_expected,
                max_cells_expected);

    double volume = compute_volume(c3t3);
    double hdist = compute_hausdorff_distance(c3t3, domain, domain_type);

    // Refine again and verify nothing changed
    std::cerr << "Refining again...\n";
    refine_mesh_3(c3t3,domain,criteria,
                  CGAL::parameters::no_exude(),
                  CGAL::parameters::no_perturb(),
                  CGAL::parameters::no_reset_c3t3());

#ifndef CGAL_MESH_3_USE_OLD_SURFACE_RESTRICTED_DELAUNAY_UPDATE
    // Using adjacencies instead of calling oracle to update restricted
    // Delaunay of the surface during the refinement of the volume
    // does not ensure that refinement is idempotent
    // BUT it should be after a few runs (all connected components should
    // have been discovered...)
    int n = 0;
    while ( c3t3.triangulation().number_of_vertices() != v && ++n < 11 )
    {
      v = c3t3.triangulation().number_of_vertices();

      refine_mesh_3(c3t3,domain,criteria,
                    CGAL::parameters::no_exude(),
                    CGAL::parameters::no_perturb(),
                    CGAL::parameters::no_reset_c3t3());
    }

    f = c3t3.number_of_facets_in_complex();
    c = c3t3.number_of_cells_in_complex(); 
    assert ( n < 11 );
#endif
    
    verify_c3t3(c3t3,domain,domain_type,v,v,f,f,c,c);
    verify_c3t3_hausdorff_distance(c3t3, domain, domain_type, hdist);

    // Exude.
    // Vertex number should not change (obvious)
    // Facet number should not change as exuder preserves boundary facets
    // Quality should increase
    C3t3 exude_c3t3(c3t3);
    std::cerr << "Exude...\n";
    CGAL::exude_mesh_3(exude_c3t3);
    verify_c3t3(exude_c3t3,domain,domain_type,v,v,f,f);
    verify_c3t3_quality(c3t3,exude_c3t3);
    verify_c3t3_volume(exude_c3t3, volume*0.95, volume*1.05);
    verify_c3t3_hausdorff_distance(exude_c3t3, domain, domain_type, hdist);
    
    // Perturb.
    // Vertex number should not change (obvious)
    // Quality should increase
    C3t3 perturb_c3t3(c3t3);
    std::cerr << "Perturb...\n";
    CGAL::perturb_mesh_3(perturb_c3t3, domain, CGAL::parameters::time_limit=5);
    verify_c3t3(perturb_c3t3,domain,domain_type,v,v);
    verify_c3t3_quality(c3t3,perturb_c3t3);
    verify_c3t3_volume(perturb_c3t3, volume*0.95, volume*1.05);
    verify_c3t3_hausdorff_distance(perturb_c3t3, domain, domain_type, hdist);

    // Odt-smoothing
    // Vertex number should not change (obvious)
    C3t3 odt_c3t3(c3t3);
    std::cerr << "Odt...\n";
    CGAL::odt_optimize_mesh_3(odt_c3t3, domain, CGAL::parameters::time_limit=5,
                              CGAL::parameters::convergence=0.001, CGAL::parameters::freeze_bound=0.0005);
    verify_c3t3(odt_c3t3,domain,domain_type,v,v);
    verify_c3t3_volume(odt_c3t3, volume*0.95, volume*1.05);
    verify_c3t3_hausdorff_distance(odt_c3t3, domain, domain_type, hdist);
    
    // Lloyd-smoothing
    // Vertex number should not change (obvious)
    C3t3 lloyd_c3t3(c3t3);
    std::cerr << "Lloyd...\n";
    CGAL::lloyd_optimize_mesh_3(lloyd_c3t3, domain, CGAL::parameters::time_limit=5,
                                CGAL::parameters::convergence=0.001, CGAL::parameters::freeze_bound=0.0005);
    verify_c3t3(lloyd_c3t3,domain,domain_type,v,v);
    verify_c3t3_volume(lloyd_c3t3, volume*0.95, volume*1.05);
    verify_c3t3_hausdorff_distance(lloyd_c3t3, domain, domain_type, hdist);
  }

  template<typename C3t3, typename Domain, typename Domain_type_tag>
  void verify_c3t3(const C3t3& c3t3,
                   const Domain& domain,
                   const Domain_type_tag domain_type,
                   const std::size_t min_vertices_expected = 0,
                   const std::size_t max_vertices_expected = STD_SIZE_T_MAX,
                   const std::size_t min_facets_expected = 0,
                   const std::size_t max_facets_expected = STD_SIZE_T_MAX,
                   const std::size_t min_cells_expected = 0,
                   const std::size_t max_cells_expected = STD_SIZE_T_MAX ) const
  {
    //-------------------------------------------------------
    // Verifications
    //-------------------------------------------------------
    std::cerr << "\tNumber of cells: " << c3t3.number_of_cells_in_complex() << "\n";
    std::cerr << "\tNumber of facets: " << c3t3.number_of_facets_in_complex() << "\n";
    std::cerr << "\tNumber of vertices: " << c3t3.triangulation().number_of_vertices() << "\n";
        
    std::size_t dist_facets ( std::distance(c3t3.facets_in_complex_begin(), 
                                            c3t3.facets_in_complex_end()) );
    std::size_t dist_cells ( std::distance(c3t3.cells_in_complex_begin(), 
                                            c3t3.cells_in_complex_end()) );

    assert(min_vertices_expected <= c3t3.triangulation().number_of_vertices());
    assert(max_vertices_expected >= c3t3.triangulation().number_of_vertices());

    assert(min_facets_expected <= c3t3.number_of_facets_in_complex());
    assert(max_facets_expected >= c3t3.number_of_facets_in_complex());
    assert(dist_facets == c3t3.number_of_facets_in_complex());
    
    assert(min_cells_expected <= c3t3.number_of_cells_in_complex());
    assert(max_cells_expected >= c3t3.number_of_cells_in_complex());
    assert(dist_cells == c3t3.number_of_cells_in_complex());
    verify_c3t3_combinatorics(c3t3, domain, domain_type);
  }
  
  template<typename C3t3>
  void verify_c3t3_quality(const C3t3& original_c3t3,
                           const C3t3& modified_c3t3) const
  {
    double original = min_value(original_c3t3);
    double modified = min_value(modified_c3t3);
    
    std::cout << "\tQuality before optimization: " << original
              << " - Quality after optimization: " << modified << std::endl;
    
    assert(original <= modified);
  }
  
  template<typename C3t3>
  double min_value(const C3t3& c3t3) const
  {
    typedef typename C3t3::Triangulation                              Tr;
    typedef typename CGAL::Mesh_3::Min_dihedral_angle_criterion<Tr>   Criterion;
    typedef typename C3t3::Cells_in_complex_iterator                  Cell_iterator;
    
    double min_value = Criterion::max_value;
    Criterion criterion(min_value, c3t3.triangulation());
    
    for ( Cell_iterator cit = c3t3.cells_in_complex_begin(),
         end = c3t3.cells_in_complex_end() ; cit != end ; ++cit )
    {
      min_value = (std::min)(min_value, criterion(c3t3.triangulation().tetrahedron(cit)));
    }
    
    return min_value;
  }

  template<typename C3t3>
  double compute_volume(const C3t3& c3t3) const
  {
    typedef typename C3t3::Cells_in_complex_iterator  Cell_iterator;

    double volume = 0.;
    for ( Cell_iterator cit = c3t3.cells_in_complex_begin(),
         end = c3t3.cells_in_complex_end() ; cit != end ; ++cit )
    {
      volume += c3t3.triangulation().tetrahedron(cit).volume();
    }
    return volume;
  }

  template<typename C3t3>
  void verify_c3t3_volume(const C3t3& c3t3,
                          const double volume_min,
                          const double volume_max) const
  {
    double volume = compute_volume(c3t3);
    std::cerr.precision(15);
    std::cerr << "\tVolume is " << volume << std::endl;

    assert(volume >= volume_min);
    assert(volume <= volume_max);
  }

  // For polyhedral domains, do nothing.
  template<typename C3t3, typename MeshDomain>
  void verify_c3t3_combinatorics(const C3t3&,
                                 const MeshDomain&,
                                 const Polyhedral_tag) const
  {}

  // For bissection domains, check the consistency between the subdomain
  // indices and the surface patch indices.
  template<typename C3t3, typename MeshDomain>
  void verify_c3t3_combinatorics(const C3t3& c3t3,
                                 const MeshDomain& domain,
                                 const Bissection_tag) const
  {
    typedef typename C3t3::Triangulation        Tr;
    typedef typename Tr::Facet                  Facet;
    typedef typename Tr::Cell_handle            Cell_handle;
    typedef typename Tr::Finite_facets_iterator Finite_facets_iterator;
    typedef typename C3t3::Surface_patch_index  Surface_patch_index;
    const Tr& tr = c3t3.triangulation();
    std::cerr.precision(17);

    for(Finite_facets_iterator fit = tr.finite_facets_begin();
      fit != tr.finite_facets_end();
      ++fit)
    {
      Facet f = *fit;
      Surface_patch_index index = c3t3.surface_patch_index(f);
      if(!c3t3.is_in_complex(f))
        assert(index == Surface_patch_index());
      else
      {
        assert(index.first != index.second);
        Cell_handle c1 = f.first;
        Cell_handle c2 = f.first->neighbor(f.second);
        if( c1->subdomain_index() == c2->subdomain_index() )
        {
          std::cerr << "ERROR:"
                    << "\nc1->subdomain_index(): " << c1->subdomain_index()
                    << "\nc2->subdomain_index(): " << c2->subdomain_index()
                    << "\nc3t3.surface_patch_index(f).first:  " << index.first
                    << "\nc3t3.surface_patch_index(f).second: " << index.second;
          if(tr.is_infinite(c1)) std::cerr << "\nc1 is infinite";
          else {
            std::cerr << "\nc1 circumcenter: " << tr.dual(c1);
            std::cerr << "\nc1 is in domain: "
                      << domain.is_in_domain_object()(tr.dual(c1));
          }
          if(tr.is_infinite(c2)) std::cerr << "\nc2 is infinite";
          else {
            std::cerr << "\nc2 circumcenter: "<< tr.dual(c2);
            std::cerr << "\nc2 is in domain: "
                      << domain.is_in_domain_object()(tr.dual(c2));
          }
          std::cerr << std::endl;
          assert(false);
        }
      }
    }
  }

  // For bissection domains, do nothing.
  template<typename C3t3, typename MeshDomain>
  double compute_hausdorff_distance(const C3t3&,
                                    const MeshDomain&,
                                    const Bissection_tag) const
  {
    return 0.;
  }

  // For polyhedral domains, compute the distance to polyhedron 
  // using the domain's AABBtree
  template<typename C3t3, typename MeshDomain>
  double compute_hausdorff_distance(const C3t3& c3t3,
                                    const MeshDomain& domain,
                                    const Polyhedral_tag) const
  {
    std::cerr.precision(17);
    typedef typename C3t3::Triangulation        Tr;
    typedef typename Tr::Facet                  Facet;
    typedef typename Tr::Finite_facets_iterator Finite_facets_iterator;
    const Tr& tr = c3t3.triangulation();
    const typename MeshDomain::AABB_tree& aabb_tree = domain.aabb_tree();

    double max_sqd = 0.;
    for(Finite_facets_iterator fit = tr.finite_facets_begin();
        fit != tr.finite_facets_end();
        ++fit)
    {
      Facet f = *fit;
      if(!c3t3.is_in_complex(f))
        continue;

      max_sqd = (std::max)(max_sqd, 
        aabb_tree.squared_distance(CGAL::centroid(tr.triangle(f))));
    }
    double hdist = std::sqrt(max_sqd);
    std::cout << "\tHausdorff distance to polyhedron is " << hdist << std::endl;
    return hdist;
  }

  template<typename C3t3, typename MeshDomain>
  void verify_c3t3_hausdorff_distance(const C3t3& c3t3,
                                      const MeshDomain& domain,
                                      const Polyhedral_tag,
                                      const double reference_value) const
  {
    double hdist = compute_hausdorff_distance(c3t3, domain, Polyhedral_tag());
#ifdef CGAL_LINKED_WITH_TBB
    // Parallel
    typedef typename C3t3::Concurrency_tag Concurrency_tag;

    if (boost::is_convertible<Concurrency_tag, CGAL::Parallel_tag>::value)
      assert(hdist <= reference_value*4.);
    else
#endif //CGAL_LINKED_WITH_TBB
      assert(hdist <= reference_value*1.01);
  }

  template<typename C3t3, typename MeshDomain>
  void verify_c3t3_hausdorff_distance(const C3t3&,
                                      const MeshDomain&,
                                      const Bissection_tag,
                                      const double) const
  { //nothing to do
  }
};

#endif // CGAL_MESH_3_TEST_TEST_MESHING_UTILITIES
