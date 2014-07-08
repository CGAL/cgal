// Copyright (c) 2010  GeometryFactory Sarl (France).
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
// Author(s)     : Laurent Rineau

#ifndef CGAL_NO_DEPRECATED_CODE

#include <boost/optional.hpp>

// The following Mesh_cell_criteria<Tr> and Mesh_facet_criteria<Tr>
// implements models of Mesh(Cell|Facet)Criteria_3 (from specification of
// CGAL-3.7), that always return "not bad".

template <typename Tr>
struct Mesh_cell_criteria {
  typedef int Cell_quality;
  typedef boost::optional<Cell_quality> Cell_badness;

  Cell_badness operator()(const typename Tr::Cell_handle) const {
    return Cell_badness();
  }
};

template <typename Tr>
struct Mesh_facet_criteria {
  typedef int Facet_quality;
  typedef boost::optional<Facet_quality> Facet_badness;

  Facet_badness operator()(const typename Tr::Facet&) const {
    return Facet_badness();
  }

  Facet_badness operator()(const typename Tr::Cell_handle, int) const {
    return Facet_badness();
  }
};

// The following Mesh_new_cell_criteria<Tr> and Mesh_new_facet_criteria<Tr>
// implements models of Mesh(Cell|Facet)Criteria_3 (from specification of
// CGAL-3.8), that always return "not bad".

template <typename Tr>
struct Mesh_new_cell_criteria {
  typedef int Cell_quality;
  typedef boost::optional<Cell_quality> Is_cell_bad;

  Is_cell_bad operator()(const typename Tr::Cell_handle) const {
    return Is_cell_bad();
  }
};

template <typename Tr>
struct Mesh_new_facet_criteria {
  typedef int Facet_quality;
  typedef boost::optional<Facet_quality> Is_facet_bad;

  Is_facet_bad operator()(const typename Tr::Facet&) const {
    return Is_facet_bad();
  }

  Is_facet_bad operator()(const typename Tr::Cell_handle, int) const {
    return Is_facet_bad();
  }
};

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Implicit_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/remove_far_points_in_mesh_3.h>

// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::FT FT;
typedef K::Point_3 Point;
typedef FT (Function)(const Point&);
typedef CGAL::Implicit_mesh_domain_3<Function,K> Mesh_domain;

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;

// Function
FT sphere_function (const Point& p)
{ return CGAL::squared_distance(p, Point(CGAL::ORIGIN))-1; }

template <typename Tr>
void test()
{
  typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

  // Criteria
  typedef CGAL::Mesh_criteria_3<Tr, 
                                CGAL::Mesh_edge_criteria_3<Tr>,
                                Mesh_facet_criteria<Tr>,
                                Mesh_cell_criteria<Tr> > Mesh_criteria;
  typedef CGAL::Mesh_criteria_3<Tr, 
                                CGAL::Mesh_edge_criteria_3<Tr>,
                                Mesh_new_facet_criteria<Tr>,
                                Mesh_new_cell_criteria<Tr> > Mesh_new_criteria;

  // Domain (Warning: Sphere_3 constructor uses squared radius !)
  Mesh_domain domain(sphere_function, K::Sphere_3(CGAL::ORIGIN, 2.));

  // Mesh criteria
  Mesh_criteria criteria = Mesh_criteria(Mesh_facet_criteria<Tr>(),
                                         Mesh_cell_criteria<Tr>());
  Mesh_new_criteria new_criteria = 
    Mesh_new_criteria(Mesh_new_facet_criteria<Tr>(),
                      Mesh_new_cell_criteria<Tr>());
  
  // Mesh generation
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, no_exude(), no_perturb());
  CGAL::refine_mesh_3<C3t3>(c3t3, domain, new_criteria, no_exude(), no_perturb());
  CGAL::remove_far_points_in_mesh_3(c3t3);

  std::cout << "Number of vertices: "
            << c3t3.triangulation().number_of_vertices() << "\n"
            << " (should be a small number, because there are no criteria"
            << " on facets and cells)\n";
}

int main()
{
  std::cout << "==== Test sequential meshing ====" << std::endl;
  typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
  test<Tr>();

#ifdef CGAL_LINKED_WITH_TBB
  std::cout << "==== Test parallel meshing ====" << std::endl;
  typedef CGAL::Mesh_triangulation_3<
    Mesh_domain,
    CGAL::Kernel_traits<Mesh_domain>::Kernel,
    CGAL::Parallel_tag>::type TrP;
  test<TrP>();
#endif

  return 0;
}

#else // CGAL_NO_DEPRECATED_CODE
#include <iostream>
int main() { 
  std::cerr << "CGAL_NO_DEPRECATED_CODE is defined. Nothing to test.\n";
  return 0; 
}
#endif // CGAL_NO_DEPRECATED_CODE
