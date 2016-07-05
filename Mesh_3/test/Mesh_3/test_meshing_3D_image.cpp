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
#include <CGAL/Image_3.h>
#include <CGAL/Labeled_image_mesh_domain_3.h>
#include <CGAL/use.h>

template <typename Concurrency_tag = CGAL::Sequential_tag>
struct Image_tester : public Tester<K_e_i>
{
public:
  void image() const
  {
    typedef CGAL::Image_3 Image;
    typedef CGAL::Labeled_image_mesh_domain_3<Image, K_e_i> Mesh_domain;

    typedef typename CGAL::Mesh_triangulation_3<
      Mesh_domain,
      CGAL::Kernel_traits<Mesh_domain>::Kernel,
      Concurrency_tag>::type Tr;
    typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

    typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
    typedef typename Mesh_criteria::Facet_criteria Facet_criteria;
    typedef typename Mesh_criteria::Cell_criteria Cell_criteria;

    //-------------------------------------------------------
    // Data generation
    //-------------------------------------------------------
    Image image;
    image.read("data/liver.inr.gz");

    std::cout << "\tSeed is\t"
      << CGAL::get_default_random().get_seed() << std::endl;
    Mesh_domain domain(image, 1e-9, &CGAL::get_default_random());

    // Set mesh criteria
    Facet_criteria facet_criteria(25, 20*image.vx(), 5*image.vx());
    Cell_criteria cell_criteria(4, 25*image.vx());
    Mesh_criteria criteria(facet_criteria, cell_criteria);

    // Mesh generation
    C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria,
                                        CGAL::parameters::no_exude(),
                                        CGAL::parameters::no_perturb());

    // Verify
    this->verify_c3t3_volume(c3t3, 1772330*0.95, 1772330*1.05);
    this->verify(c3t3,domain,criteria, Bissection_tag());

    typedef typename Mesh_domain::Surface_patch_index Patch_id;
    CGAL_static_assertion(CGAL::Output_rep<Patch_id>::is_specialized);
    CGAL_USE_TYPE(Patch_id);
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
