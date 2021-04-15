// Copyright (c) 2008  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Monique Teillaud, Sylvain Pion, Pedro Machado

// Partially supported by the IST Programme of the EU as a
// STREP (FET Open) Project under Contract No  IST-006413
// (ACS -- Algorithms for Complex Shapes)

template <class SK>
void _test_spherical_kernel_compute(SK sk)
{

  typedef typename SK::FT                                          FT;

  typedef typename SK::Circular_arc_point_3                        Circular_arc_point_3;
  typedef typename SK::Circular_arc_3                              Circular_arc_3;


  typedef typename SK::Circle_3                                    Circle_3;

  typedef typename SK::Algebraic_kernel                            AK;


  typedef typename SK::Construct_circle_3                          Construct_circle_3;

  typedef typename SK::Construct_circular_arc_3                    Construct_circular_arc_3;
  typedef typename SK::Construct_circular_arc_point_3              Construct_circular_arc_point_3;



  typedef typename SK::Compute_approximate_squared_length_3         Compute_approximate_squared_length_3;
  typedef typename SK::Compute_approximate_angle_3                   Compute_approximate_angle_3;
  typedef typename SK::Polynomials_for_circle_3                    Polynomials_for_circle_3;
  typedef typename AK::Polynomial_for_spheres_2_3                  Polynomial_for_spheres_2_3;
  typedef typename AK::Polynomial_1_3                              Polynomial_1_3;

  typedef typename AK::Root_for_spheres_2_3                        Root_for_spheres_2_3;

  std::cout << "TESTING COMPUTATIONS" << std::endl;

  (void)/* Equal_3 theEqual_3 = */ sk.equal_3_object();
  (void)/* Get_equation theGet_equation = */ sk.get_equation_object();
  Construct_circle_3 theConstruct_circle_3 = sk.construct_circle_3_object();
  (void)/* Construct_sphere_3 theConstruct_sphere_3 = */ sk.construct_sphere_3_object();
  Construct_circular_arc_3 theConstruct_circular_arc_3 = sk.construct_circular_arc_3_object();
  Construct_circular_arc_point_3 theConstruct_circular_arc_point_3 = sk.construct_circular_arc_point_3_object();
  (void)/* Compute_area_divided_by_pi_3 theCompute_area_divided_by_pi_3 = */ sk.compute_area_divided_by_pi_3_object();
  (void)/* Compute_squared_length_divided_by_pi_square_3 theCompute_squared_length_divided_by_pi_square_3 = */
    sk.compute_squared_length_divided_by_pi_square_3_object();
  (void)/* Compute_approximate_area_3 theCompute_approximate_area_3 = */ sk.compute_approximate_area_3_object();
  Compute_approximate_squared_length_3 theCompute_approximate_squared_length_3 =
    sk.compute_approximate_squared_length_3_object();
  Compute_approximate_angle_3 theCompute_approximate_angle_3 = sk.compute_approximate_angle_3_object();

  std::cout << "Testing Approximate_angle of a Circular_arc_3" << std::endl;
  std::cout << "Testing Approximate_squared_length of a Circular_arc_3" << std::endl;

  Root_for_spheres_2_3 rt[8];

  rt[0] = Root_for_spheres_2_3(0,1,0);
  rt[1] = Root_for_spheres_2_3(CGAL::make_root_of_2(FT(0),FT(-FT(1)/FT(2)),FT(2)), CGAL::make_root_of_2(FT(0),FT(FT(1)/FT(2)),FT(2)),0);
  rt[2] = Root_for_spheres_2_3(-1,0,0);
  rt[3] = Root_for_spheres_2_3(CGAL::make_root_of_2(FT(0),FT(-FT(1)/FT(2)),FT(2)), CGAL::make_root_of_2(FT(0),FT(-FT(1)/FT(2)),FT(2)),0);
  rt[4] = Root_for_spheres_2_3(0,-1,0);
  rt[5] = Root_for_spheres_2_3(CGAL::make_root_of_2(FT(0),FT(FT(1)/FT(2)),FT(2)), CGAL::make_root_of_2(FT(0),FT(-FT(1)/FT(2)),FT(2)),0);
  rt[6] = Root_for_spheres_2_3(1,0,0);
  rt[7] = Root_for_spheres_2_3(CGAL::make_root_of_2(FT(0),FT(FT(1)/FT(2)),FT(2)), CGAL::make_root_of_2(FT(0),FT(FT(1)/FT(2)),FT(2)),0);

  Circular_arc_point_3 cp[8];
  for(int i=0; i<8; i++) {
    cp[i] = theConstruct_circular_arc_point_3(rt[i]);
  }

  const double pi = CGAL_PI;

  const Polynomials_for_circle_3 pcc_test =
      std::make_pair(Polynomial_for_spheres_2_3(0,0,0,1),
                     Polynomial_1_3(0,0,1,0));
  Circle_3 cc = theConstruct_circle_3(pcc_test);
  for(int i=0; i<8; i++) {
    for(int j=i+1; j<8; j++) {
      Circular_arc_3 ca = theConstruct_circular_arc_3(cc,cp[i],cp[j]);
      Circular_arc_3 cb = theConstruct_circular_arc_3(cc,cp[j],cp[i]);
      const double num = (double) (j-i);
      const double ang1 = (pi/4.0) * num;
      const double ang2 = 2*pi - ang1;
      const double app_ang1 = theCompute_approximate_angle_3(ca);
      const double app_ang2 = theCompute_approximate_angle_3(cb);
      const double v1 = app_ang1 - ang1;
      const double v2 = app_ang2 - ang2;
      const double diff1 = ((v1 > 0) ? (v1) : (-v1));
      const double diff2 = ((v2 > 0) ? (v2) : (-v2));
      // we suppose at least a precision of 10e-4, but it is not necessarily true
      assert(diff1 < 10e-4);
      assert(diff2 < 10e-4);

      const double sql1 = ang1 * ang1;
      const double sql2 = ang2 * ang2;
      const double app_sql1 = theCompute_approximate_squared_length_3(ca);
      const double app_sql2 = theCompute_approximate_squared_length_3(cb);
      const double vv1 = app_sql1 - sql1;
      const double vv2 = app_sql2 - sql2;
      const double diffv1 = ((vv1 > 0) ? (vv1) : (-vv1));
      const double diffv2 = ((vv2 > 0) ? (vv2) : (-vv2));
      std::cout << sql1 << " " << app_sql1 << std::endl;
      std::cout << sql2 << " " << app_sql2 << std::endl;
      // we suppose at least a precision of 10e-4, but it is not necessarily true
      assert(diffv1 < 10e-4);
      assert(diffv2 < 10e-4);
    }
  }

  std::cout << "All tests on computations are OK." << std::endl;
}
