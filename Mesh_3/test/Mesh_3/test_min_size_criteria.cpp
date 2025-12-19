// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Jane Tournois
//
//******************************************************************************
// File Description :
//
//******************************************************************************

#include "test_meshing_utilities.h"

#include <CGAL/Image_3.h>
#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/use.h>
#include <CGAL/Mesh_criteria_3.h>

template <typename Concurrency_tag = CGAL::Sequential_tag>
struct Image_tester : public Tester<K_e_i>
{
public:
  void image() const
  {
    using Image = CGAL::Image_3;
    using Mesh_domain = CGAL::Labeled_mesh_domain_3<K_e_i>;

    using Tr = typename CGAL::Mesh_triangulation_3<
      Mesh_domain,
      CGAL::Kernel_traits<Mesh_domain>::Kernel,
      Concurrency_tag>::type;
    using C3t3 = CGAL::Mesh_complex_3_in_triangulation_3<Tr>;

    using GT             = typename Tr::Geom_traits;
    using FT             = typename Tr::Geom_traits::FT;
    using Bare_point     = typename Tr::Bare_point;

    using Mesh_criteria  = CGAL::Mesh_criteria_3<Tr>;
    using Facet_criteria = typename Mesh_criteria::Facet_criteria;
    using Cell_criteria  = typename Mesh_criteria::Cell_criteria;

    //-------------------------------------------------------
    // Data generation
    //-------------------------------------------------------
    Image image;
    image.read(CGAL::data_file_path("images/liver.inr.gz"));

    std::cout << "\tSeed is\t"
      << CGAL::get_default_random().get_seed() << std::endl;
    Mesh_domain domain = Mesh_domain::create_labeled_image_mesh_domain
       (image,
        CGAL::parameters::relative_error_bound = 1e-9,
        CGAL::parameters::p_rng = &CGAL::get_default_random());

    // Set mesh criteria
    const double fangle = 25;
    const double fsize = 20 * image.vx();
    const double fapprox = image.vx();
    const CGAL::Mesh_facet_topology ftopo = CGAL::FACET_VERTICES_ON_SURFACE;
    const double fminsize = 0.5 * image.vx();

    const double cshape = 4.;
    const double csize = 25 * image.vx();
    const double cminsize = 0.5 * image.vx();

    Facet_criteria facet_criteria(fangle, fsize, fapprox, ftopo, fminsize);
    Cell_criteria cell_criteria(cshape, csize, cminsize);
    Mesh_criteria criteria(facet_criteria, cell_criteria);

    // Mesh generation
    C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria,
                                        CGAL::parameters::no_exude(),
                                        CGAL::parameters::no_perturb());

    const Tr& tr = c3t3.triangulation();
    typename GT::Construct_point_3 cp = tr.geom_traits().construct_point_3_object();
    typename GT::Compute_squared_radius_3 sq_radius = tr.geom_traits().compute_squared_radius_3_object();

    double max_sq_facet_radius = 0.;
    double max_sq_cell_radius = 0.;

    for (auto f : c3t3.facets_in_complex())
    {
      const Bare_point p1 = cp(tr.point(f.first, (f.second + 1) & 3));
      const Bare_point& ball_center = f.first->get_facet_surface_center(f.second);

      const FT sqr = tr.min_squared_distance(p1, ball_center);
      max_sq_facet_radius = (std::max)(max_sq_facet_radius, sqr);
    }

    for (auto c : c3t3.cells_in_complex())
    {
      const Bare_point& p = cp(tr.point(c, 0));
      const Bare_point& q = cp(tr.point(c, 1));
      const Bare_point& r = cp(tr.point(c, 2));
      const Bare_point& s = cp(tr.point(c, 3));

      const FT sqr = sq_radius(p, q, r, s);
      max_sq_cell_radius = (std::max)(max_sq_cell_radius, sqr);
    }

    const std::size_t nbv = c3t3.triangulation().number_of_vertices();
    std::cout << "C3t3 initial = " << nbv << std::endl;

    //new criteria with really small facet_size and cell_size
    //compared to actual elements size
    const double c3t3_fsize = CGAL::approximate_sqrt(max_sq_facet_radius);
    const double c3t3_csize = CGAL::approximate_sqrt(max_sq_cell_radius);

    //use the other construction API
    Mesh_criteria new_criteria(CGAL::parameters::facet_size     = 0.01 * c3t3_fsize,
                               CGAL::parameters::facet_min_size = 1.00001 * c3t3_fsize,
                               CGAL::parameters::cell_size      = 0.01 * c3t3_csize,
                               CGAL::parameters::cell_min_size  = 1.00001 * c3t3_csize);

    CGAL::refine_mesh_3(c3t3, domain, new_criteria,
                        CGAL::parameters::no_perturb(),
                        CGAL::parameters::no_exude());

    const std::size_t nbv2 = c3t3.triangulation().number_of_vertices();
    std::cout << "C3t3 after refinement = " << nbv2 << std::endl;

    assert(nbv == nbv2);
  }
};



int main()
{
  Image_tester<> test_epic;
  std::cerr << "Mesh generation from a 3D image:\n";
  test_epic.image();

#ifdef CGAL_LINKED_WITH_TBB
  Image_tester<CGAL::Parallel_tag> test_epic_p;
  std::cerr << "Parallel mesh generation from a 3D image:\n";
  test_epic_p.image();
#endif

  return EXIT_SUCCESS;
}
