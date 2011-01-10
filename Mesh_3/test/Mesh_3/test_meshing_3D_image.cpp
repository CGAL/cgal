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
// File Description : 
//******************************************************************************

#include "test_meshing_utilities.h"
#include <CGAL/Image_3.h>
#include <CGAL/Labeled_image_mesh_domain_3.h>

template <typename K>
struct Image_tester : public Tester<K>
{
public:
  void image() const
  {
    typedef CGAL::Image_3 Image;
    typedef CGAL::Labeled_image_mesh_domain_3<Image, K> Mesh_domain;
    
    typedef typename CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
    typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;
    
    typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
    typedef typename Mesh_criteria::Facet_criteria Facet_criteria;
    typedef typename Mesh_criteria::Cell_criteria Cell_criteria;
    
    typedef typename Mesh_domain::Surface_patch_index Surface_patch_index;
    
    //-------------------------------------------------------
    // Data generation
    //-------------------------------------------------------
    Image image;
    image.read("data/liver.inr.gz");
    Mesh_domain domain(image,1e-9);
    
    // Set mesh criteria
    Facet_criteria facet_criteria(25, 20*image.vx(), 5*image.vx());
    Cell_criteria cell_criteria(4, 25*image.vx());
    Mesh_criteria criteria(facet_criteria, cell_criteria);
    
    // Mesh generation
    C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria,
                                        CGAL::parameters::no_exude(),
                                        CGAL::parameters::no_perturb());
    
    // Verify
    this->verify(c3t3,domain,criteria);
  }
};



int main()
{
  Image_tester<K_e_i> test_epic;
  std::cerr << "Mesh generation from a 3D image:\n";
  test_epic.image();
  
  return EXIT_SUCCESS;  
}
