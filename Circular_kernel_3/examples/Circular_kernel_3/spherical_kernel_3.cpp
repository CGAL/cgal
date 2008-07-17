// Copyright (c) 2005-2006  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473
// (ECG - Effective Computational Geometry for Curves and Surfaces)
// and a STREP (FET Open) Project under Contract No  IST-006413
// (ACS -- Algorithms for Complex Shapes)
//
// $URL$
// $Id$
//
// Author(s) : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//             Sylvain Pion     <Sylvain.Pion@sophia.inria.fr>
//             Pedro Machado    <tashimir@gmail.com>

#include <CGAL/Cartesian.h>
#include <CGAL/Algebraic_kernel_for_spheres_2_3.h>
#include <CGAL/Spherical_kernel_3.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Random.h>

typedef CGAL::Quotient< CGAL::MP_Float >                    NT;
typedef CGAL::Cartesian<NT>                                 Linear_k1;
typedef CGAL::Algebraic_kernel_for_spheres_2_3<NT>          Algebraic_k1;
typedef CGAL::Spherical_kernel_3<Linear_k1,Algebraic_k1>    SK;
typedef SK::FT                                              FT;
typedef SK::Sphere_3                                        Sphere_3;
typedef SK::Circle_3                                        Circle_3;
typedef SK::Intersect_3                                     Intersect_3;
typedef SK::Construct_sphere_3                              Construct_sphere_3;
typedef SK::Algebraic_kernel                                AK;
typedef AK::Polynomial_for_spheres_2_3                      Polynomial_for_spheres_2_3;

int main() {

  Intersect_3 theIntersect_3 = SK().intersect_3_object();
  Construct_sphere_3 theConstruct_sphere_3 = SK().construct_sphere_3_object();

  CGAL::Random generatorOfgenerator;
  int random_seed = generatorOfgenerator.get_int(0, 123456);
  CGAL::Random theRandom(random_seed);
  int count = 0;

  std::cout << "We will calcule the approximate probability that 3 spheres with radius 1 intersect on a 5x5x5 box, it may take some time." << std::endl;

  for(int i=0; i<10000; i++) {
    double x1 = theRandom.get_double(0.0,5.0);
    double y1 = theRandom.get_double(0.0,5.0);
    double z1 = theRandom.get_double(0.0,5.0);
    double r = 1.0;
    double x2 = theRandom.get_double(0.0,5.0);
    double y2 = theRandom.get_double(0.0,5.0);
    double z2 = theRandom.get_double(0.0,5.0);
    double x3 = theRandom.get_double(0.0,5.0);
    double y3 = theRandom.get_double(0.0,5.0);
    double z3 = theRandom.get_double(0.0,5.0);
    Sphere_3 s1 = theConstruct_sphere_3(Polynomial_for_spheres_2_3(FT(x1),FT(y1),FT(z1),FT(r)));
    Sphere_3 s2 = theConstruct_sphere_3(Polynomial_for_spheres_2_3(FT(x2),FT(y2),FT(z2),FT(r)));
    Sphere_3 s3 = theConstruct_sphere_3(Polynomial_for_spheres_2_3(FT(x3),FT(y3),FT(z3),FT(r)));
    std::vector< CGAL::Object > intersection_1;
    theIntersect_3(s1, s2, s3, std::back_inserter(intersection_1));
    if(intersection_1.size() > 0) count++;
  }

  std::cout << "The approximate probability that 3 spheres with radius 1"
            << std::endl;
  std::cout << "choosen (uniformly) randomly on a 5x5x5 box intersect is: "
            << ((double)count)/((double)(10000)) << std::endl;

  return 0;
}
