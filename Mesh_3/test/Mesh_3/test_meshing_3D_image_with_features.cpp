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
// Author(s)     : Stephane Tayeb, Jane Tournois
//
//******************************************************************************
// File Description :
//******************************************************************************

#include "test_meshing_utilities.h"
#include <CGAL/Image_3.h>
#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/Mesh_domain_with_polyline_features_3.h>
#include <CGAL/Mesh_3/Detect_features_in_image.h>
#include <CGAL/Mesh_3/Detect_features_on_image_bbox.h>
#include <CGAL/use.h>

template <typename Concurrency_tag = CGAL::Sequential_tag>
struct Image_tester : public Tester<K_e_i>
{
    typedef CGAL::Image_3 Image;
    typedef CGAL::Labeled_mesh_domain_3<K_e_i> Domain;
    typedef CGAL::Mesh_domain_with_polyline_features_3<Domain> Mesh_domain;

    typedef typename CGAL::Mesh_triangulation_3<
      Mesh_domain,
      CGAL::Kernel_traits<Mesh_domain>::Kernel,
      Concurrency_tag>::type Tr;
    typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

    typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

  void mesh_and_verify(Mesh_domain& domain, const Image& image, const double volume) const
  {
    namespace p = CGAL::parameters;
    // Set mesh criteria
    Mesh_criteria criteria(p::edge_size = 2 * image.vx(),
                           p::facet_angle = 30,
                           p::facet_size = 20 * image.vx(),
                           p::facet_distance = 5 * image.vx(),
                           p::cell_radius_edge_ratio = 3.,
                           p::cell_size = 25 * image.vx());

    // Mesh generation
    C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria,
                                        CGAL::parameters::no_exude(),
                                        CGAL::parameters::no_perturb());

    c3t3.remove_isolated_vertices();

    // Verify
    this->verify_c3t3_volume(c3t3, volume * 0.95, volume * 1.05);
    this->verify(c3t3, domain, criteria, Bissection_tag());

    typedef typename Mesh_domain::Surface_patch_index Patch_id;
    CGAL_static_assertion(CGAL::Output_rep<Patch_id>::is_specialized);
    CGAL_USE_TYPE(Patch_id);
  }

public:
  void image() const
  {
    namespace p = CGAL::parameters;
    std::cout << "\tSeed is\t"
      << CGAL::get_default_random().get_seed() << std::endl;

    //-------------------------------------------------------
    // Data generation
    //-------------------------------------------------------
    Image image;
    image.read(CGAL::data_file_path("images/liver.inr.gz"));

    Mesh_domain domain = Mesh_domain::create_labeled_image_mesh_domain
     (p::image = image,
      p::relative_error_bound = 1e-9,
      CGAL::parameters::p_rng = &CGAL::get_default_random(),
      CGAL::parameters::detect_features = CGAL::Mesh_3::Detect_features_in_image());

    mesh_and_verify(domain, image, 1772330.);
  }

  void image_in_bbox() const
  {
    namespace p = CGAL::parameters;
    std::cout << "\tSeed is\t"
      << CGAL::get_default_random().get_seed() << std::endl;

    //-------------------------------------------------------
    // Data generation
    //-------------------------------------------------------
    Image image;
    image.read(CGAL::data_file_path("images/40420.inr"));

    Mesh_domain domain = Mesh_domain::create_labeled_image_mesh_domain
     (p::image = image,
      p::relative_error_bound = 1e-9,
      CGAL::parameters::p_rng = &CGAL::get_default_random(),
      CGAL::parameters::detect_features = CGAL::Mesh_3::Detect_features_on_image_bbox());

    mesh_and_verify(domain, image, 625044.);
  }
};



int main()
{
  Image_tester<> test_epic;
  std::cerr << "Mesh generation from a 3D image"
            << " with detection of triple lines:\n";
  test_epic.image();

  std::cerr << "Mesh generation from a 3D image"
            << " with detection of triple lines on bbox only:\n";
  test_epic.image_in_bbox();

#ifdef CGAL_LINKED_WITH_TBB
  Image_tester<CGAL::Parallel_tag> test_epic_p;
  std::cerr << "Parallel mesh generation from a 3D image"
            << " with detection of triple lines:\n";
  test_epic_p.image();

  std::cerr << "Parallel mesh generation from a 3D image"
            << " with detection of triple lines on bbox only:\n";
  test_epic_p.image_in_bbox();
#endif

  return EXIT_SUCCESS;
}
