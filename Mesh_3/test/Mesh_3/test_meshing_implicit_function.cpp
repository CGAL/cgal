// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
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
// File Description : 
//******************************************************************************

#include "test_meshing_utilities.h"
#include <CGAL/Implicit_mesh_domain_3.h>

template <typename K, typename Concurrency_tag = CGAL::Sequential_tag>
struct Implicit_tester : public Tester<K>
{
  typedef typename K::Point_3 Point;
  typedef typename K::FT FT;
  static FT sphere_function (const Point& p)
  {
    const FT x2=p.x()*p.x(), y2=p.y()*p.y(), z2=p.z()*p.z();
    return x2+y2+z2-1;
  }
  
  void implicit() const
  {
    
    typedef FT (Function)(const Point&);
    
    typedef CGAL::Implicit_mesh_domain_3<Function, K> Mesh_domain;
    
    typedef typename CGAL::Mesh_triangulation_3<
      Mesh_domain,
      typename CGAL::Kernel_traits<Mesh_domain>::Kernel,
      Concurrency_tag>::type Tr;
    typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;
    
    typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
    typedef typename Mesh_criteria::Facet_criteria Facet_criteria;
    typedef typename Mesh_criteria::Cell_criteria Cell_criteria;
    
    typedef typename K::Sphere_3 Sphere_3;
    
    typedef typename Mesh_domain::Surface_patch_index Surface_patch_index;
    
    //-------------------------------------------------------
    // Data generation
    //-------------------------------------------------------
    std::cout << "\tSeed is\t" 
      << CGAL::get_default_random().get_seed() << std::endl;
    Mesh_domain domain(Implicit_tester<K>::sphere_function,
                       Sphere_3(CGAL::ORIGIN,2.),
                       1e-3,
                       &CGAL::get_default_random());
    
    // Set mesh criteria
    Facet_criteria facet_criteria(0, 0, 0.3);
    Cell_criteria cell_criteria(0, 0.5);
    Mesh_criteria criteria(facet_criteria, cell_criteria);
    
    std::vector<Point> initial_points;
    initial_points.push_back(Point(1,0,0));
    initial_points.push_back(Point(0,1,0));
    initial_points.push_back(Point(0,0,1));
    initial_points.push_back(Point(-1,0,0));
    initial_points.push_back(Point(0,-1,0));
    initial_points.push_back(Point(0,0,-1));
    
    // Mesh generation
    C3t3 c3t3;
    c3t3.insert_surface_points(initial_points.begin(),
                               initial_points.end(),
                               domain.index_from_surface_patch_index(Surface_patch_index(0,1)));
    
    CGAL::refine_mesh_3(c3t3, domain, criteria,
                        CGAL::parameters::no_exude(),
                        CGAL::parameters::no_perturb());

    CGAL::remove_far_points_in_mesh_3(c3t3);

#ifdef CGAL_LINKED_WITH_TBB
    // Parallel
    if (boost::is_convertible<Concurrency_tag, CGAL::Parallel_tag>::value)
    {
      this->verify(c3t3, domain, criteria, Bissection_tag(), 40, 65, 60, 110);
    }
    else
#endif //CGAL_LINKED_WITH_TBB
    {
      // Verify
      this->verify(c3t3, domain, criteria, Bissection_tag(), 50, 58, 80, 90);
    }
  }
};


int main()
{
  Implicit_tester<K_e_i> test_epic;
  std::cerr << "Mesh generation from an implicit function:\n";
  test_epic.implicit();
  
#ifdef CGAL_LINKED_WITH_TBB
  Implicit_tester<K_e_i, CGAL::Parallel_tag> test_epic_p;
  std::cerr << "Parallel mesh generation from an implicit function:\n";
  test_epic_p.implicit();
#endif
  return EXIT_SUCCESS;
}
