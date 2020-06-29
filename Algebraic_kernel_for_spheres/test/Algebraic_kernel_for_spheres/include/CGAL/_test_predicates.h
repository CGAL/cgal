// Copyright (c) 2005-2006  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473
// (ECG - Effective Computational Geometry for Curves and Surfaces)
// and a STREP (FET Open) Project under Contract No  IST-006413
// (ACS -- Algorithms for Complex Shapes)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//             Sylvain Pion
//             Pedro Machado

#include <CGAL/Random.h>
#include <cassert>
#include <CGAL/Polynomials_1_3.h>

template <class AK>
void _test_solve(AK ak)
{
  typedef typename AK::FT FT;
  typedef typename AK::Root_of_2 Root_of_2;
  typedef typename AK::Root_for_spheres_2_3 Root_for_spheres_2_3;
  typename AK::Solve theSolve = ak.solve_object();
  typename AK::Construct_polynomial_for_spheres_2_3 theConstruct_2_3 =
    ak.construct_polynomial_for_spheres_2_3_object();
  typename AK::Construct_polynomial_1_3 theConstruct_1_3 =
    ak.construct_polynomial_1_3_object();
  typename AK::Sign_at theSigh_at =
    ak.sign_at_object();

  //Polynomial_for_spheres_2_3

  std::vector< std::pair<Root_for_spheres_2_3, size_t> > res_sss1;
  theSolve(theConstruct_2_3(5, 5, 0, 25),
           theConstruct_2_3(0, 5, 0, 100),
           theConstruct_2_3(7, 5, 0, 9),
           std::back_inserter(res_sss1));
  assert(res_sss1.size() == 1);
  assert(res_sss1[0].second == 2u);
  assert(res_sss1[0].first == Root_for_spheres_2_3(10, 5, 0));

  std::vector< std::pair<Root_for_spheres_2_3, size_t> > res_sss2;
  theSolve(theConstruct_2_3(0, 5, 0, 100),
           theConstruct_2_3(20, 5, 0, 100),
           theConstruct_2_3(10, 5, 0, 25),
           std::back_inserter(res_sss2));
  assert(res_sss2.size() == 0);

  std::vector< std::pair<Root_for_spheres_2_3, size_t> > res_sss3;
  theSolve(theConstruct_2_3(0, 5, 0, 100),
           theConstruct_2_3(20, 5, 0, 100),
           theConstruct_2_3(10, 10, 0, 25),
           std::back_inserter(res_sss3));
  assert(res_sss3.size() == 1);
  assert(res_sss3[0].second == 2u);
  assert(res_sss3[0].first == Root_for_spheres_2_3(10, 5, 0));

  std::vector< std::pair<Root_for_spheres_2_3, size_t> > res_sss4;
  theSolve(theConstruct_2_3(0, 5, 0, 100),
           theConstruct_2_3(20, 5, 0, 100),
           theConstruct_2_3(10, 5, 5, 25),
           std::back_inserter(res_sss4));
  assert(res_sss4.size() == 1);
  assert(res_sss4[0].second == 2u);
  assert(res_sss4[0].first == Root_for_spheres_2_3(10, 5, 0));

  std::vector< std::pair<Root_for_spheres_2_3, size_t> > res_sss5;
  theSolve(theConstruct_2_3(0, 5, 0, 100),
           theConstruct_2_3(20, 5, 0, 100),
           theConstruct_2_3(10, 0, 5, 25),
           std::back_inserter(res_sss5));
  assert(res_sss5.size() == 0);

  std::vector< std::pair<Root_for_spheres_2_3, size_t> > res_sss6;
  theSolve(theConstruct_2_3(5, 0, 0, 100),
           theConstruct_2_3(0, 5, 0, 100),
           theConstruct_2_3(0, 0, 5, 100),
           std::back_inserter(res_sss6));
  assert(res_sss6.size() == 2);
  Root_of_2 res_sss6_r1 = make_root_of_2(FT(5,3),-FT(5,3),FT(10));
  Root_of_2 res_sss6_r2 = make_root_of_2(FT(5,3),FT(5,3),FT(10));
  assert(res_sss6[0].second == 1u);
  assert(res_sss6[0].first == Root_for_spheres_2_3(res_sss6_r1, res_sss6_r1, res_sss6_r1));
  assert(res_sss6[1].second == 1u);
  assert(res_sss6[1].first == Root_for_spheres_2_3(res_sss6_r2, res_sss6_r2, res_sss6_r2));

  std::vector< std::pair<Root_for_spheres_2_3, size_t> > res_sss_7;
  theSolve(theConstruct_2_3(5, 0, 0, 25),
           theConstruct_2_3(0, 5, 0, 25),
           theConstruct_2_3(0, 0, 5, 25),
           std::back_inserter(res_sss_7));
  assert(res_sss_7.size() == 2);
  assert(res_sss_7[0].second == 1u);
  assert(res_sss_7[0].first == Root_for_spheres_2_3(0, 0, 0));
  assert(res_sss_7[1].second == 1u);
  assert(res_sss_7[1].first == Root_for_spheres_2_3(Root_of_2(FT(10,3)),
                                               Root_of_2(FT(10,3)),
                                               Root_of_2(FT(10,3))));

  std::vector< std::pair<Root_for_spheres_2_3, size_t> > res_sss_8;
  theSolve(theConstruct_2_3(0, 5, 0, 100),
           theConstruct_2_3(20, 5, 0, 100),
           theConstruct_2_3(20, 5, 0, 100),
           std::back_inserter(res_sss_8));
  assert(res_sss_8.size() == 1);
  assert(res_sss_8[0].second == 2u);
  assert(res_sss_8[0].first == Root_for_spheres_2_3(10, 5, 0));

  std::vector< std::pair<Root_for_spheres_2_3, size_t> > res_sss9;
  theSolve(theConstruct_2_3(0, 5, 0, 100),
           theConstruct_2_3(20, 5, 0, 100),
           theConstruct_2_3(0, 5, 0, 125),
           std::back_inserter(res_sss9));
  assert(res_sss9.size() == 0);

  //Polynomial_1_3 Polynomial_for_spheres_2_3

  // 2 spheres and a plane
  std::vector< std::pair<Root_for_spheres_2_3, size_t> > res_ssp_1;
  theSolve(theConstruct_2_3(0, 5, 0, 100),
           theConstruct_2_3(20, 5, 0, 100),
           theConstruct_1_3(0, 0, 1, 0),
           std::back_inserter(res_ssp_1));
  assert(res_ssp_1.size() == 1);
  assert(res_ssp_1[0].second == 2u);
  assert(res_ssp_1[0].first == Root_for_spheres_2_3(10, 5, 0));

  std::vector< std::pair<Root_for_spheres_2_3, size_t> > res_ssp_2;
  theSolve(theConstruct_2_3(0, 5, 0, 100),
           theConstruct_2_3(20, 5, 0, 100),
           theConstruct_1_3(1, -3, 3, 5),
           std::back_inserter(res_ssp_2));
  assert(res_ssp_2.size() == 1);
  assert(res_ssp_2[0].second == 2u);
  assert(res_ssp_2[0].first == Root_for_spheres_2_3(10, 5, 0));

  std::vector< std::pair<Root_for_spheres_2_3, size_t> > res_ssp_3;
  theSolve(theConstruct_2_3(0, 5, 0, 29),
           theConstruct_2_3(5, 0, 0, 29),
           theConstruct_1_3(0, 0, 1, 5),
           std::back_inserter(res_ssp_3));
  assert(res_ssp_3.size() == 0);

  std::vector< std::pair<Root_for_spheres_2_3, size_t> > res_ssp_4;
  theSolve(theConstruct_2_3(0, 5, 0, FT(75,2)),
           theConstruct_2_3(5, 0, 0, FT(75,2)),
           theConstruct_1_3(0, 0, 1, -5),
           std::back_inserter(res_ssp_4));
  assert(res_ssp_4.size() == 1);
  assert(res_ssp_4[0].second == 2u);
  assert(res_ssp_4[0].first == Root_for_spheres_2_3(FT(5,2), FT(5,2), 5));

  std::vector< std::pair<Root_for_spheres_2_3, size_t> > res_ssp_5;
  theSolve(theConstruct_2_3(0, 5, 0, 50),
           theConstruct_2_3(5, 0, 0, 50),
           theConstruct_1_3(0, 0, 1, -5),
           std::back_inserter(res_ssp_5));
  assert(res_ssp_5.size() == 2);
  assert(res_ssp_5[0].second == 1u);
  assert(res_ssp_5[0].first == Root_for_spheres_2_3(0,0,5));
  assert(res_ssp_5[1].second == 1u);
  assert(res_ssp_5[1].first == Root_for_spheres_2_3(5, 5, 5));

  std::vector< std::pair<Root_for_spheres_2_3, size_t> > res_ssp_6;
  theSolve(theConstruct_2_3(5, 5, 5, 25),
           theConstruct_2_3(5, 5, 5, 25),
           theConstruct_1_3(0, 0, 1, 0),
           std::back_inserter(res_ssp_6));
  assert(res_ssp_6.size() == 1);
  assert(res_ssp_6[0].second == 2u);
  assert(res_ssp_6[0].first == Root_for_spheres_2_3(5, 5, 0));

  std::vector< std::pair<Root_for_spheres_2_3, size_t> > res_ssp_7;
  theSolve(theConstruct_2_3(5, 5, 5, 25),
           theConstruct_2_3(3, 5, 7, 100),
           theConstruct_1_3(1, 1, 1, 10),
           std::back_inserter(res_ssp_7));
  assert(res_ssp_7.size() == 0);

  std::vector< std::pair<Root_for_spheres_2_3, size_t> > res_ssp_8;
  theSolve(theConstruct_2_3(5, 5, 5, 25),
           theConstruct_2_3(3, 5, 7, 40),
           theConstruct_1_3(1, 1, 1, -10),
           std::back_inserter(res_ssp_8));
  assert(res_ssp_8.size() == 2);
  assert(res_ssp_8[0].second == 1u);
  assert(res_ssp_8[1].second == 1u);
  Root_of_2 res_ssp_8_x1 = make_root_of_2(FT(101,24),FT(1,24),FT(1453));
  Root_of_2 res_ssp_8_x2 = make_root_of_2(FT(101,24),-FT(1,24),FT(1453));
  Root_of_2 res_ssp_8_y1 = make_root_of_2(FT(10,3),-FT(1,12),FT(1453));
  Root_of_2 res_ssp_8_y2 = make_root_of_2(FT(10,3),FT(1,12),FT(1453));
  Root_of_2 res_ssp_8_z1 = make_root_of_2(FT(59,24),FT(1,24),FT(1453));
  Root_of_2 res_ssp_8_z2 = make_root_of_2(FT(59,24),-FT(1,24),FT(1453));
  Root_for_spheres_2_3 res_ssp_8_sol2 = Root_for_spheres_2_3(res_ssp_8_x1,
                                                             res_ssp_8_y1,
                                                             res_ssp_8_z1);
  Root_for_spheres_2_3 res_ssp_8_sol1 = Root_for_spheres_2_3(res_ssp_8_x2,
                                                             res_ssp_8_y2,
                                                             res_ssp_8_z2);
  assert(res_ssp_8[0].first == res_ssp_8_sol1);
  assert(res_ssp_8[1].first == res_ssp_8_sol2);

  std::vector< std::pair<Root_for_spheres_2_3, size_t> > res_ssp_9;
  theSolve(theConstruct_2_3(5, 5, 5, 25),
           theConstruct_2_3(3, 5, 7, 100),
           theConstruct_1_3(1, 1, 1, -10),
           std::back_inserter(res_ssp_9));
  assert(res_ssp_9.size() == 0);

  std::vector< std::pair<Root_for_spheres_2_3, size_t> > res_ssp_10;
  theSolve(theConstruct_2_3(5, 0, 0, 100),
           theConstruct_2_3(0, 0, 0, 100),
           theConstruct_1_3(1, 2, 3, 100),
           std::back_inserter(res_ssp_10));
  assert(res_ssp_10.size() == 0);

  // 1 sphere and 2 plane
  std::vector< std::pair<Root_for_spheres_2_3, size_t> > res_spp_1;
  theSolve(theConstruct_2_3(5, 10, 10, 5000),
           theConstruct_1_3(2, 4, 6, 50),
           theConstruct_1_3(1, 2, 3, 100),
           std::back_inserter(res_spp_1));
  assert(res_spp_1.size() == 0);

  std::vector< std::pair<Root_for_spheres_2_3, size_t> > res_spp_2;
  theSolve(theConstruct_2_3(0, 0, 0, 9),
           theConstruct_1_3(1, 0, 0, -3),
           theConstruct_1_3(2, 0, 0, -6),
           std::back_inserter(res_spp_2));
  assert(res_spp_2.size() == 1);
  assert(res_spp_2[0].second == 2u);
  assert(res_spp_2[0].first == Root_for_spheres_2_3(3, 0, 0));

  std::vector< std::pair<Root_for_spheres_2_3, size_t> > res_spp_3;
  theSolve(theConstruct_2_3(0, 0, 0, 1),
           theConstruct_1_3(1, 0, 0, 0),
           theConstruct_1_3(0, 0, 1, 0),
           std::back_inserter(res_spp_3));
  assert(res_spp_3.size() == 2);
  assert(res_spp_3[0].second == 1u);
  assert(res_spp_3[0].first == Root_for_spheres_2_3(0,-1,0));
  assert(res_spp_3[1].second == 1u);
  assert(res_spp_3[1].first == Root_for_spheres_2_3(0,1,0));

  std::vector< std::pair<Root_for_spheres_2_3, size_t> > res_spp_4;
  theSolve(theConstruct_2_3(0, 0, 0, 1),
           theConstruct_1_3(1, 1, 0, -1),
           theConstruct_1_3(1, -1, 0, 0),
           std::back_inserter(res_spp_4));
  assert(res_spp_4.size() == 2);
  assert(res_spp_4[0].second == 1u);
  assert(res_spp_4[0].first == Root_for_spheres_2_3(FT(1,2),FT(1,2),
                                 make_root_of_2(FT(0),-FT(1,2),FT(2))));
  assert(res_spp_4[1].second == 1u);
  assert(res_spp_4[1].first == Root_for_spheres_2_3(FT(1,2),FT(1,2),
                                 make_root_of_2(FT(0),FT(1,2),FT(2))));

  std::vector< std::pair<Root_for_spheres_2_3, size_t> > res_spp_5;
  theSolve(theConstruct_2_3(10, 10, 10, 25),
           theConstruct_1_3(1, 0, 0, 0),
           theConstruct_1_3(0, 1, 0, 0),
           std::back_inserter(res_spp_5));
  assert(res_spp_5.size() == 0);

  // 1 circle and 1 plane && 1 plane and 1 circle

  // circle - plane
  std::vector< std::pair<Root_for_spheres_2_3, size_t> > res_cp_1;
  theSolve(std::make_pair(theConstruct_2_3(0,0,0,1),
                          theConstruct_1_3(0,1,0,0)),
                          theConstruct_1_3(0,1,0,10),
           std::back_inserter(res_cp_1));
  assert(res_cp_1.size() == 0);

  std::vector< std::pair<Root_for_spheres_2_3, size_t> > res_cp_2;
  theSolve(std::make_pair(theConstruct_2_3(0,0,0,1),
                          theConstruct_1_3(0,1,0,0)),
                          theConstruct_1_3(1,0,0,1),
           std::back_inserter(res_cp_2));
  assert(res_cp_2.size() == 1);
  assert(res_cp_2[0].second == 2u);
  assert(res_cp_2[0].first == Root_for_spheres_2_3(-1, 0, 0));

  std::vector< std::pair<Root_for_spheres_2_3, size_t> > res_cp_3;
  theSolve(std::make_pair(theConstruct_2_3(0,0,0,1),
                          theConstruct_1_3(0,1,0,0)),
                          theConstruct_1_3(1,0,0,0),
           std::back_inserter(res_cp_3));
  assert(res_cp_3.size() == 2);
  assert(res_cp_3[0].second == 1u);
  assert(res_cp_3[0].first ==
          Root_for_spheres_2_3(0, 0, -1));
  assert(res_cp_3[1].second == 1u);
  assert(res_cp_3[1].first ==
          Root_for_spheres_2_3(0, 0, 1));

  // plane - circle
  std::vector< std::pair<Root_for_spheres_2_3, size_t> > res_pc_1;
  theSolve(theConstruct_1_3(0,1,0,10),
               std::make_pair(theConstruct_2_3(0,0,0,1),
                                      theConstruct_1_3(0,1,0,0)),
                   std::back_inserter(res_pc_1));
  assert(res_pc_1.size() == 0);

  std::vector< std::pair<Root_for_spheres_2_3, size_t> > res_pc_2;
  theSolve(theConstruct_1_3(1,0,0,1),
           std::make_pair(theConstruct_2_3(0,0,0,1),
                                          theConstruct_1_3(0,1,0,0)),
                   std::back_inserter(res_pc_2));
  assert(res_pc_2.size() == 1);
  assert(res_pc_2[0].second == 2u);
  assert(res_pc_2[0].first == Root_for_spheres_2_3(-1, 0, 0));

  std::vector< std::pair<Root_for_spheres_2_3, size_t> > res_pc_3;
  theSolve(theConstruct_1_3(1,0,0,0),
           std::make_pair(theConstruct_2_3(0,0,0,1),
                                  theConstruct_1_3(0,1,0,0)),
                   std::back_inserter(res_pc_3));
  assert(res_pc_3.size() == 2);
  assert(res_pc_3[0].second == 1u);
  assert(res_pc_3[0].first ==
                  Root_for_spheres_2_3(0, 0, -1));
  assert(res_pc_3[1].second == 1u);
  assert(res_pc_3[1].first ==
                  Root_for_spheres_2_3(0, 0, 1));

  // 1 circle and 1 sphere && 1 sphere and 1 circle

  // circle - sphere
  std::vector< std::pair<Root_for_spheres_2_3, size_t> > res_cs_1;
  theSolve(std::make_pair(theConstruct_2_3(0,0,0,1),
                          theConstruct_1_3(0,1,0,0)),
                          theConstruct_2_3(2,2,2,1),
           std::back_inserter(res_cs_1));
  assert(res_cs_1.size() == 0);

  std::vector< std::pair<Root_for_spheres_2_3, size_t> > res_cs_2;
  theSolve(std::make_pair(theConstruct_2_3(0,0,0,1),
                          theConstruct_1_3(0,1,0,0)),
                          theConstruct_2_3(2,0,0,1),
           std::back_inserter(res_cs_2));
  assert(res_cs_2.size() == 1);
  assert(res_cs_2[0].second == 2u);
  assert(res_cs_2[0].first == Root_for_spheres_2_3(1, 0, 0));

  std::vector< std::pair<Root_for_spheres_2_3, size_t> > res_cs_3;
  theSolve(std::make_pair(theConstruct_2_3(0,0,0,1),
                          theConstruct_1_3(0,1,0,0)),
                          theConstruct_2_3(1,0,0,1),
           std::back_inserter(res_cs_3));
  assert(res_cs_3.size() == 2);
  assert(res_cs_3[0].second == 1u);
  assert(res_cs_3[0].first ==
          Root_for_spheres_2_3(FT(1,2), 0, make_root_of_2(FT(0),-FT(1,2),FT(3))));
  assert(res_cs_3[1].second == 1u);
  assert(res_cs_3[1].first ==
          Root_for_spheres_2_3(FT(1,2), 0, make_root_of_2(FT(0),FT(1,2),FT(3))));

  // sphere - circle
  std::vector< std::pair<Root_for_spheres_2_3, size_t> > res_sc_1;
  theSolve(theConstruct_2_3(2,2,2,1),
               std::make_pair(theConstruct_2_3(0,0,0,1),
                                      theConstruct_1_3(0,1,0,0)),
                   std::back_inserter(res_sc_1));
  assert(res_sc_1.size() == 0);

  std::vector< std::pair<Root_for_spheres_2_3, size_t> > res_sc_2;
  theSolve(theConstruct_2_3(2,0,0,1),
           std::make_pair(theConstruct_2_3(0,0,0,1),
                                          theConstruct_1_3(0,1,0,0)),
                   std::back_inserter(res_sc_2));
  assert(res_sc_2.size() == 1);
  assert(res_sc_2[0].second == 2u);
  assert(res_sc_2[0].first == Root_for_spheres_2_3(1, 0, 0));

  std::vector< std::pair<Root_for_spheres_2_3, size_t> > res_sc_3;
  theSolve(theConstruct_2_3(1,0,0,1),
           std::make_pair(theConstruct_2_3(0,0,0,1),
                                          theConstruct_1_3(0,1,0,0)),
                   std::back_inserter(res_sc_3));
  assert(res_sc_3.size() == 2);
  assert(res_sc_3[0].second == 1u);
  assert(res_sc_3[0].first ==
                  Root_for_spheres_2_3(FT(1,2), 0, make_root_of_2(FT(0),-FT(1,2),FT(3))));
  assert(res_sc_3[1].second == 1u);
  assert(res_sc_3[1].first ==
                  Root_for_spheres_2_3(FT(1,2), 0, make_root_of_2(FT(0),FT(1,2),FT(3))));

  // 2 circles
  std::vector< std::pair<Root_for_spheres_2_3, size_t> > res_cc_1;
  theSolve(std::make_pair(theConstruct_2_3(0,0,0,1),
                          theConstruct_1_3(0,1,0,0)),
           std::make_pair(theConstruct_2_3(3,0,0,1),
                          theConstruct_1_3(0,1,0,0)),
           std::back_inserter(res_cc_1));
  assert(res_cc_1.size() == 0);

  std::vector< std::pair<Root_for_spheres_2_3, size_t> > res_cc_2;
  theSolve(std::make_pair(theConstruct_2_3(0,0,0,1),
                          theConstruct_1_3(0,1,0,0)),
           std::make_pair(theConstruct_2_3(2,0,0,1),
                          theConstruct_1_3(0,1,0,0)),
           std::back_inserter(res_cc_2));
  assert(res_cc_2.size() == 1);
  assert(res_cc_2[0].second == 2u);
  assert(res_cc_2[0].first == Root_for_spheres_2_3(1, 0, 0));

  std::vector< std::pair<Root_for_spheres_2_3, size_t> > res_cc_3;
  theSolve(std::make_pair(theConstruct_2_3(0,0,0,4),
                          theConstruct_1_3(0,1,0,0)),
           std::make_pair(theConstruct_2_3(2,0,0,4),
                          theConstruct_1_3(0,1,0,0)),
           std::back_inserter(res_cc_3));
  assert(res_cc_3.size() == 2);
  assert(res_cc_3[0].second == 1u);
  assert(res_cc_3[0].first ==
          Root_for_spheres_2_3(1, 0, make_root_of_2(FT(0),FT(-1),FT(3))));
  assert(res_cc_3[1].second == 1u);
  assert(res_cc_3[1].first ==
          Root_for_spheres_2_3(1, 0, make_root_of_2(FT(0),FT(1),FT(3))));

  // A QUESTION: WHAT IS THE MULTIPLICITY OF A CIRCLE vs CIRCLE SOLUTION?
  // FOR NOW, IF THE NUMBER OF SOLUTION IS 1, THE MULTIPLICITY IS 2
  //          IF THE NUMBER OF SOLUTION IS 2, THE MULTIPLICITY IS 1
  std::vector< std::pair<Root_for_spheres_2_3, size_t> > res_cc_4;
  theSolve(std::make_pair(theConstruct_2_3(0,0,0,4),
                          theConstruct_1_3(0,1,0,0)),
           std::make_pair(theConstruct_2_3(2,0,0,4),
                          theConstruct_1_3(3,5,0,-3)),
           std::back_inserter(res_cc_4));
  assert(res_cc_4.size() == 2);
  assert(res_cc_4[0].second == 1u);
  assert(res_cc_4[0].first ==
          Root_for_spheres_2_3(1, 0, make_root_of_2(FT(0),FT(-1),FT(3))));
  assert(res_cc_4[1].second == 1u);
  assert(res_cc_4[1].first ==
          Root_for_spheres_2_3(1, 0, make_root_of_2(FT(0),FT(1),FT(3))));

  std::vector< std::pair<Root_for_spheres_2_3, size_t> > res_cc_5;
  theSolve(std::make_pair(theConstruct_2_3(0,0,0,169),
                          theConstruct_1_3(0,1,0,0)),
           std::make_pair(theConstruct_2_3(24,0,0,169),
                          theConstruct_1_3(0,1,0,0)),
           std::back_inserter(res_cc_5));
  assert(res_cc_5.size() == 2);
  assert(res_cc_5[0].second == 1u);
  assert(res_cc_5[0].first == Root_for_spheres_2_3(12, 0, -5));
  assert(res_cc_5[1].second == 1u);
  assert(res_cc_5[1].first == Root_for_spheres_2_3(12, 0, 5));

  std::vector< std::pair<Root_for_spheres_2_3, size_t> > res_cc_6;
  theSolve(std::make_pair(theConstruct_2_3(0,0,0,169),
                          theConstruct_1_3(0,1,0,0)),
           std::make_pair(theConstruct_2_3(24,0,0,169),
                          theConstruct_1_3(1,3,1,-17)),
           std::back_inserter(res_cc_6));
  assert(res_cc_6.size() == 1);
  assert(res_cc_6[0].second == 2u);
  assert(res_cc_6[0].first == Root_for_spheres_2_3(12, 0, 5));

  std::vector< std::pair<Root_for_spheres_2_3, size_t> > res_cc_7;
  theSolve(std::make_pair(theConstruct_2_3(0,0,0,169),
                          theConstruct_1_3(0,0,1,-1)),
           std::make_pair(theConstruct_2_3(0,0,0,169),
                          theConstruct_1_3(0,0,1,1)),
           std::back_inserter(res_cc_7));
  assert(res_cc_7.size() == 0);

  std::vector< std::pair<Root_for_spheres_2_3, size_t> > res_cc_8;
  theSolve(std::make_pair(theConstruct_2_3(0,0,0,169),
                          theConstruct_1_3(1,1,0,-5)),
           std::make_pair(theConstruct_2_3(0,0,0,169),
                          theConstruct_1_3(1,-1,0,5)),
           std::back_inserter(res_cc_8));
  assert(res_cc_8.size() == 2);
  assert(res_cc_8[0].second == 1u);
  assert(res_cc_8[0].first == Root_for_spheres_2_3(0, 5, -12));
  assert(res_cc_8[1].second == 1u);
  assert(res_cc_8[1].first == Root_for_spheres_2_3(0, 5, 12));

  std::vector< std::pair<Root_for_spheres_2_3, size_t> > res_cc_9;
  theSolve(std::make_pair(theConstruct_2_3(0,0,0,169),
                          theConstruct_1_3(0,1,0,-13)),
           std::make_pair(theConstruct_2_3(0,0,0,169),
                          theConstruct_1_3(1,0,0,0)),
           std::back_inserter(res_cc_9));
  assert(res_cc_9.size() == 1);
  assert(res_cc_9[0].second == 2u);
  assert(res_cc_9[0].first == Root_for_spheres_2_3(0, 13, 0));

  std::vector< std::pair<Root_for_spheres_2_3, size_t> > res_cc_10;
  theSolve(std::make_pair(theConstruct_2_3(0,0,0,169),
                          theConstruct_1_3(2,1,-1,-13)),
           std::make_pair(theConstruct_2_3(0,0,0,169),
                          theConstruct_1_3(1,0,0,0)),
           std::back_inserter(res_cc_10));
  assert(res_cc_10.size() == 2);
  assert(res_cc_10[0].second == 1u);
  assert(res_cc_10[0].first == Root_for_spheres_2_3(0, 0, -13));
  assert(res_cc_10[1].second == 1u);
  assert(res_cc_10[1].first == Root_for_spheres_2_3(0, 13, 0));

  std::vector< std::pair<Root_for_spheres_2_3, size_t> > res_cc_11;
  theSolve(std::make_pair(theConstruct_2_3(0,0,0,169),
                          theConstruct_1_3(2,1,0,-13)),
           std::make_pair(theConstruct_2_3(0,0,0,169),
                          theConstruct_1_3(1,0,0,0)),
           std::back_inserter(res_cc_11));
  assert(res_cc_11.size() == 1);
  assert(res_cc_11[0].second == 2u);
  assert(res_cc_11[0].first == Root_for_spheres_2_3(0, 13, 0));

  // only Polynomial_1_3
  std::vector< std::pair<Root_for_spheres_2_3, size_t> > res_p12;
  theSolve(theConstruct_1_3(1, 1, 0, -5),
           theConstruct_1_3(1, -1, 0, 5),
           theConstruct_1_3(0, 0, 1, 0),
           std::back_inserter(res_p12));
  assert(res_p12.size() == 1);
  assert(res_p12[0].second == 1u);
  assert(res_p12[0].first == Root_for_spheres_2_3(0, 5, 0));

  std::vector< std::pair<Root_for_spheres_2_3, size_t> > res_p13;
  theSolve(theConstruct_1_3(0, 1, 0, -5),
           theConstruct_1_3(1, -1, 0, 5),
           theConstruct_1_3(0, 0, 1, 0),
           std::back_inserter(res_p13));
  assert(res_p13.size() == 1);
  assert(res_p13[0].second == 1u);
  assert(res_p13[0].first == Root_for_spheres_2_3(0, 5, 0));

  std::vector< std::pair<Root_for_spheres_2_3, size_t> > res_p14;
  theSolve(theConstruct_1_3(1, -1, 0, 5),
           theConstruct_1_3(0, 1, 0, -5),
           theConstruct_1_3(0, 0, 1, 0),
           std::back_inserter(res_p14));
  assert(res_p14.size() == 1);
  assert(res_p14[0].second == 1u);
  assert(res_p14[0].first == Root_for_spheres_2_3(0, 5, 0));

  std::vector< std::pair<Root_for_spheres_2_3, size_t> > res_p15;
  theSolve(theConstruct_1_3(1, 0, 0, 0),
           theConstruct_1_3(1, -1, 0, 5),
           theConstruct_1_3(0, 0, 1, 0),
           std::back_inserter(res_p15));
  assert(res_p15.size() == 1);
  assert(res_p15[0].second == 1u);
  assert(res_p15[0].first == Root_for_spheres_2_3(0, 5, 0));

  std::vector< std::pair<Root_for_spheres_2_3, size_t> > res_p16;
  theSolve(theConstruct_1_3(1, -1, 0, 5),
           theConstruct_1_3(1, 0, 0, 0),
           theConstruct_1_3(0, 0, 1, 0),
           std::back_inserter(res_p16));
  assert(res_p16.size() == 1);
  assert(res_p16[0].second == 1u);
  assert(res_p16[0].first == Root_for_spheres_2_3(0, 5, 0));

  std::vector< std::pair<Root_for_spheres_2_3, size_t> > res_p17;
  theSolve(theConstruct_1_3(0, 1, 0, -5),
           theConstruct_1_3(1, 0, 0, 0),
           theConstruct_1_3(0, 0, 1, 0),
           std::back_inserter(res_p17));
  assert(res_p17.size() == 1);
  assert(res_p17[0].second == 1u);
  assert(res_p17[0].first == Root_for_spheres_2_3(0, 5, 0));

  std::vector< std::pair<Root_for_spheres_2_3, size_t> > res_p18;
  theSolve(theConstruct_1_3(1, 0, 0, 0),
           theConstruct_1_3(0, 1, 0, -5),
           theConstruct_1_3(0, 0, 1, 0),
           std::back_inserter(res_p18));
  assert(res_p18.size() == 1);
  assert(res_p18[0].second == 1u);
  assert(res_p18[0].first == Root_for_spheres_2_3(0, 5, 0));

  std::vector< std::pair<Root_for_spheres_2_3, size_t> > res_p20;
  theSolve(theConstruct_1_3(1, 0, 0, 0),
           theConstruct_1_3(1, 0, 0, 5),
           theConstruct_1_3(0, 0, 1, 0),
           std::back_inserter(res_p20));
  assert(res_p20.size() == 0);

  std::vector< std::pair<Root_for_spheres_2_3, size_t> > res_p22;
  theSolve(theConstruct_1_3(0, 1, 0, -5),
           theConstruct_1_3(0, 1, 0, 0),
           theConstruct_1_3(0, 0, 1, 0),
           std::back_inserter(res_p22));
  assert(res_p22.size() == 0);

  std::vector< std::pair<Root_for_spheres_2_3, size_t> > res_p24;
  theSolve(theConstruct_1_3(1, -1, 0, 5),
           theConstruct_1_3(2, -2, 0, 15),
           theConstruct_1_3(0, 0, 1, 0),
           std::back_inserter(res_p24));
  assert(res_p24.size() == 0);

  std::vector< std::pair<Root_for_spheres_2_3, size_t> > res_p25;
  theSolve(theConstruct_1_3(1, 0, 0, -5),
           theConstruct_1_3(0, 1, 0, -3),
           theConstruct_1_3(0, 0, 1, -2),
           std::back_inserter(res_p25));
  assert(res_p25.size() == 1);
  assert(res_p25[0].second == 1u);
  assert(res_p25[0].first == Root_for_spheres_2_3(5, 3, 2));

  std::vector< std::pair<Root_for_spheres_2_3, size_t> > res_p26;
  theSolve(theConstruct_1_3(3, 1, 1, -18),
           theConstruct_1_3(1, 2, 1, -16),
           theConstruct_1_3(1, 1, 1, -12),
           std::back_inserter(res_p26));
  assert(res_p26.size() == 1);
  assert(res_p26[0].second == 1u);
  assert(res_p26[0].first == Root_for_spheres_2_3(3, 4, 5));

  CGAL::Random generatorOfgenerator;
  int random_seed = generatorOfgenerator.get_int(0, 123456);
  std::cout << "random_seed = " << random_seed << std::endl;
  CGAL::Random theRandom(random_seed);
  int random_max = 5;
  int random_min = -5;

   typedef typename AK::Polynomial_for_spheres_2_3
      Polynomial_for_spheres_2_3;
   typedef typename AK::Polynomial_1_3
      Polynomial_1_3;



   typedef typename AK::Polynomials_for_line_3 Polynomials_for_line_3;

  // only Polynomial_1_3
  for(std::size_t i = 0; i < 500; i++){
    int a1,a2,a3,a4,b1,b2,b3,b4,c1,c2,c3,c4;
    Polynomial_1_3 p1, p2, p3;
    do {
      a1 = theRandom.get_int(random_min,random_max);
      a2 = theRandom.get_int(random_min,random_max);
      a3 = theRandom.get_int(random_min,random_max);
      a4 = theRandom.get_int(random_min,random_max);
      b1 = theRandom.get_int(random_min,random_max);
      b2 = theRandom.get_int(random_min,random_max);
      b3 = theRandom.get_int(random_min,random_max);
      b4 = theRandom.get_int(random_min,random_max);
      c1 = theRandom.get_int(random_min,random_max);
      c2 = theRandom.get_int(random_min,random_max);
      c3 = theRandom.get_int(random_min,random_max);
      c4 = theRandom.get_int(random_min,random_max);
      if(a1 == 0 && a2 == 0 && a3 == 0) continue;
      if(b1 == 0 && b2 == 0 && b3 == 0) continue;
      if(c1 == 0 && c2 == 0 && c3 == 0) continue;
      p1 = Polynomial_1_3(a1,a2,a3,a4);
      p2 = Polynomial_1_3(b1,b2,b3,b4);
      p3 = Polynomial_1_3(c1,c2,c3,c4);
    } while(CGAL::same_solutions<FT>(p1,p2) ||
            CGAL::same_solutions<FT>(p1,p3) ||
            CGAL::same_solutions<FT>(p2,p3));
                std::cout << a1 << " " << a2 << " " << a3 << " " << a4 << std::endl;
                std::cout << b1 << " " << b2 << " " << b3 << " " << b4 << std::endl;
          std::cout << c1 << " " << c2 << " " << c3 << " " << c4 << std::endl;
    std::vector< std::pair<Root_for_spheres_2_3, size_t> > res_ppp;
    theSolve(p1,p2,p3,std::back_inserter(res_ppp));
    if(res_ppp.size() != 0) {
      assert(theSigh_at(p1,res_ppp[0].first) == CGAL::ZERO);
      assert(theSigh_at(p2,res_ppp[0].first) == CGAL::ZERO);
      assert(theSigh_at(p3,res_ppp[0].first) == CGAL::ZERO);
    }
  }

  // only spheres
  for(std::size_t i = 0; i < 500; i++){
    int a1,a2,a3,a4,b1,b2,b3,b4,c1,c2,c3,c4;
    Polynomial_for_spheres_2_3 s1, s2, s3;
    do {
      a1 = theRandom.get_int(random_min,random_max);
      a2 = theRandom.get_int(random_min,random_max);
      a3 = theRandom.get_int(random_min,random_max);
      a4 = theRandom.get_int(random_min,random_max);
      b1 = theRandom.get_int(random_min,random_max);
      b2 = theRandom.get_int(random_min,random_max);
      b3 = theRandom.get_int(random_min,random_max);
      b4 = theRandom.get_int(random_min,random_max);
      c1 = theRandom.get_int(random_min,random_max);
      c2 = theRandom.get_int(random_min,random_max);
      c3 = theRandom.get_int(random_min,random_max);
      c4 = theRandom.get_int(random_min,random_max);
      s1 = Polynomial_for_spheres_2_3(a1,a2,a3,a4);
      s2 = Polynomial_for_spheres_2_3(b1,b2,b3,b4);
      s3 = Polynomial_for_spheres_2_3(c1,c2,c3,c4);
    } while((a4 <= 0) || (b4 <= 0) || (c4 <= 0) ||
            (s1 == s2) || (s2 == s3) || (s1 == s3));

                if(CGAL::AlgebraicSphereFunctors::intersect<AK>(s1,s2)) {
      Polynomial_1_3 p1 = CGAL::AlgebraicSphereFunctors::plane_from_2_spheres<AK>(s1,s2);
                        if(CGAL::AlgebraicSphereFunctors::intersect<AK>(s2,s3)) {
                          Polynomial_1_3 p2 = CGAL::AlgebraicSphereFunctors::plane_from_2_spheres<AK>(s2,s3);
                          if(CGAL::same_solutions<FT>(p1,p2)) {
                            const FT sq_d1 = CGAL::square(p1.a()*s1.a() + p1.b()*s1.b() +
                                                   p1.c()*s1.c() + p1.d()) /
                                     (CGAL::square(p1.a()) + CGAL::square(p1.b()) + CGAL::square(p1.c()));
                                        const FT r1_sqr = s1.r_sq() - sq_d1;
                const FT sq_d2 = CGAL::square(p2.a()*s2.a() + p2.b()*s2.b() +
                                                                p2.c()*s2.c() + p2.d()) /
                                                                                (CGAL::square(p2.a()) + CGAL::square(p2.b()) + CGAL::square(p2.c()));
                                        const FT r2_sqr = s2.r_sq() - sq_d2;
                                        if(r1_sqr == r2_sqr) if(r1_sqr != 0) { --i; continue; }
        }
      }
    }

                std::cout << a1 << " " << a2 << " " << a3 << " " << a4 << std::endl;
                std::cout << b1 << " " << b2 << " " << b3 << " " << b4 << std::endl;
                std::cout << c1 << " " << c2 << " " << c3 << " " << c4 << std::endl;
    std::vector< std::pair<Root_for_spheres_2_3, size_t> > res_sss;
    theSolve(s1,s2,s3,std::back_inserter(res_sss));
    if(res_sss.size() == 1) {
      assert(theSigh_at(s1,res_sss[0].first) == CGAL::ZERO);
      assert(theSigh_at(s2,res_sss[0].first) == CGAL::ZERO);
      assert(theSigh_at(s3,res_sss[0].first) == CGAL::ZERO);
    }
    if(res_sss.size() == 2) {
      assert(theSigh_at(s1,res_sss[0].first) == CGAL::ZERO);
      assert(theSigh_at(s2,res_sss[0].first) == CGAL::ZERO);
      assert(theSigh_at(s3,res_sss[0].first) == CGAL::ZERO);
      assert(theSigh_at(s1,res_sss[1].first) == CGAL::ZERO);
      assert(theSigh_at(s2,res_sss[1].first) == CGAL::ZERO);
      assert(theSigh_at(s3,res_sss[1].first) == CGAL::ZERO);
      if(res_sss[0].first.x() != res_sss[1].first.x()) {
        assert(res_sss[0].first.x() < res_sss[1].first.x());
      } else if(res_sss[0].first.y() != res_sss[1].first.y()) {
        assert(res_sss[0].first.y() < res_sss[1].first.y());
      } else {
        assert(res_sss[0].first.z() < res_sss[1].first.z());
      }
    }
  }

  // two spheres and a plane
  for(std::size_t i = 0; i < 500; i++){
    int a1,a2,a3,a4,b1,b2,b3,b4,c1,c2,c3,c4;
    Polynomial_for_spheres_2_3 s1, s2;
    Polynomial_1_3 p;
    do {
      a1 = theRandom.get_int(random_min,random_max);
      a2 = theRandom.get_int(random_min,random_max);
      a3 = theRandom.get_int(random_min,random_max);
      a4 = theRandom.get_int(random_min,random_max);
      b1 = theRandom.get_int(random_min,random_max);
      b2 = theRandom.get_int(random_min,random_max);
      b3 = theRandom.get_int(random_min,random_max);
      b4 = theRandom.get_int(random_min,random_max);
      c1 = theRandom.get_int(random_min,random_max);
      c2 = theRandom.get_int(random_min,random_max);
      c3 = theRandom.get_int(random_min,random_max);
      c4 = theRandom.get_int(random_min,random_max);
      s1 = Polynomial_for_spheres_2_3(a1,a2,a3,a4);
      s2 = Polynomial_for_spheres_2_3(b1,b2,b3,b4);
      p = Polynomial_1_3(c1,c2,c3,c4);
    } while((a4 <= 0) || (b4 <= 0) || (s1 == s2) ||
            (c1 == 0 && c2 == 0 && c3 == 0));

    if(CGAL::AlgebraicSphereFunctors::intersect<AK>(s1,s2)) {
      Polynomial_1_3 p1 = CGAL::AlgebraicSphereFunctors::plane_from_2_spheres<AK>(s1,s2);
      if(CGAL::same_solutions<FT>(p1,p)) {
        const FT sq_d1 = CGAL::square(p1.a()*s1.a() + p1.b()*s1.b() +p1.c()*s1.c() + p1.d()) /
                (CGAL::square(p1.a()) + CGAL::square(p1.b()) + CGAL::square(p1.c()));
        const FT r1_sqr = s1.r_sq() - sq_d1;
            if(r1_sqr != 0) { --i; continue; }
      }
    }


                std::cout << a1 << " " << a2 << " " << a3 << " " << a4 << std::endl;
    std::cout << b1 << " " << b2 << " " << b3 << " " << b4 << std::endl;
    std::cout << c1 << " " << c2 << " " << c3 << " " << c4 << std::endl;
    std::vector< std::pair<Root_for_spheres_2_3, size_t> > res_ssp;
    theSolve(p,s1,s2,std::back_inserter(res_ssp));
    if(res_ssp.size() == 1) {
      assert(theSigh_at(s1,res_ssp[0].first) == CGAL::ZERO);
      assert(theSigh_at(s2,res_ssp[0].first) == CGAL::ZERO);
      assert(theSigh_at(p,res_ssp[0].first) == CGAL::ZERO);
    }
    if(res_ssp.size() == 2) {
      assert(theSigh_at(s1,res_ssp[0].first) == CGAL::ZERO);
      assert(theSigh_at(s2,res_ssp[0].first) == CGAL::ZERO);
      assert(theSigh_at(p,res_ssp[0].first) == CGAL::ZERO);
      assert(theSigh_at(s1,res_ssp[1].first) == CGAL::ZERO);
      assert(theSigh_at(s2,res_ssp[1].first) == CGAL::ZERO);
      assert(theSigh_at(p,res_ssp[1].first) == CGAL::ZERO);
      if(res_ssp[0].first.x() != res_ssp[1].first.x()) {
        assert(res_ssp[0].first.x() < res_ssp[1].first.x());
      } else if(res_ssp[0].first.y() != res_ssp[1].first.y()) {
        assert(res_ssp[0].first.y() < res_ssp[1].first.y());
      } else {
        assert(res_ssp[0].first.z() < res_ssp[1].first.z());
      }
    }
  }

  // two planes and a sphere
  for(std::size_t i = 0; i < 500; i++){
    int a1,a2,a3,a4,b1,b2,b3,b4,c1,c2,c3,c4;
    Polynomial_for_spheres_2_3 s;
    Polynomial_1_3 p1, p2;
    do {
      a1 = theRandom.get_int(random_min,random_max);
      a2 = theRandom.get_int(random_min,random_max);
      a3 = theRandom.get_int(random_min,random_max);
      a4 = theRandom.get_int(random_min,random_max);
      b1 = theRandom.get_int(random_min,random_max);
      b2 = theRandom.get_int(random_min,random_max);
      b3 = theRandom.get_int(random_min,random_max);
      b4 = theRandom.get_int(random_min,random_max);
      c1 = theRandom.get_int(random_min,random_max);
      c2 = theRandom.get_int(random_min,random_max);
      c3 = theRandom.get_int(random_min,random_max);
      c4 = theRandom.get_int(random_min,random_max);
      p1 = Polynomial_1_3(a1,a2,a3,a4);
      p2 = Polynomial_1_3(b1,b2,b3,b4);
      s = Polynomial_for_spheres_2_3(c1,c2,c3,c4);
    } while((a1 == 0 && a2 == 0 && a3 == 0) ||
            (b1 == 0 && b2 == 0 && b3 == 0) ||
            (c4 <= 0) || (CGAL::same_solutions<FT>(p1,p2)));

    if(CGAL::same_solutions<FT>(p1,p2)) {
      const FT sq_d1 = CGAL::square(p1.a()*s.a() + p1.b()*s.b() +p1.c()*s.c() + p1.d()) /
              (CGAL::square(p1.a()) + CGAL::square(p1.b()) + CGAL::square(p1.c()));
      const FT r1_sqr = s.r_sq() - sq_d1;
      if(r1_sqr != 0) { --i; continue; }
    }

                std::cout << a1 << " " << a2 << " " << a3 << " " << a4 << std::endl;
                std::cout << b1 << " " << b2 << " " << b3 << " " << b4 << std::endl;
                std::cout << c1 << " " << c2 << " " << c3 << " " << c4 << std::endl;
    std::vector< std::pair<Root_for_spheres_2_3, size_t> > res_spp;
    theSolve(p1,p2,s,std::back_inserter(res_spp));
    if(res_spp.size() == 1) {
      assert(theSigh_at(p1,res_spp[0].first) == CGAL::ZERO);
      assert(theSigh_at(p2,res_spp[0].first) == CGAL::ZERO);
      assert(theSigh_at(s,res_spp[0].first) == CGAL::ZERO);
    }
    if(res_spp.size() == 2) {
      assert(theSigh_at(p1,res_spp[0].first) == CGAL::ZERO);
      assert(theSigh_at(p2,res_spp[0].first) == CGAL::ZERO);
      assert(theSigh_at(s,res_spp[0].first) == CGAL::ZERO);
      assert(theSigh_at(p1,res_spp[1].first) == CGAL::ZERO);
      assert(theSigh_at(p2,res_spp[1].first) == CGAL::ZERO);
      assert(theSigh_at(s,res_spp[1].first) == CGAL::ZERO);
      if(res_spp[0].first.x() != res_spp[1].first.x()) {
        assert(res_spp[0].first.x() < res_spp[1].first.x());
      } else if(res_spp[0].first.y() != res_spp[1].first.y()) {
        assert(res_spp[0].first.y() < res_spp[1].first.y());
      } else {
        assert(res_spp[0].first.z() < res_spp[1].first.z());
      }
    }
  }

  // A line and a sphere
  for(std::size_t i = 0; i < 500; i++){
    int a1,a2,a3,b1,b2,b3,c1,c2,c3,c4;
    Polynomials_for_line_3 l;
    Polynomial_for_spheres_2_3 s;
    do {
      a1 = theRandom.get_int(random_min,random_max);
      a2 = theRandom.get_int(random_min,random_max);
      a3 = theRandom.get_int(random_min,random_max);
      b1 = theRandom.get_int(random_min,random_max);
      b2 = theRandom.get_int(random_min,random_max);
      b3 = theRandom.get_int(random_min,random_max);
      c1 = theRandom.get_int(random_min,random_max);
      c2 = theRandom.get_int(random_min,random_max);
      c3 = theRandom.get_int(random_min,random_max);
      c4 = theRandom.get_int(random_min,random_max);
      l = Polynomials_for_line_3(a1,b1,a2,b2,a3,b3);
      s = Polynomial_for_spheres_2_3(c1,c2,c3,c4);
    } while((l.degenerated()) || (c4 <= 0));
                std::cout << a1 << " " << a2 << " " << a3 << std::endl;
                std::cout << b1 << " " << b2 << " " << b3 << std::endl;
                std::cout << c1 << " " << c2 << " " << c3 << std::endl;
    std::vector< std::pair<Root_for_spheres_2_3, size_t> > res_ls;
    theSolve(l,s,std::back_inserter(res_ls));
    if(res_ls.size() == 1) {
      assert(theSigh_at(s,res_ls[0].first) == CGAL::ZERO);
      assert(CGAL::AlgebraicSphereFunctors::contains<AK>(l,res_ls[0].first));
    }
    if(res_ls.size() == 2) {
      assert(theSigh_at(s,res_ls[0].first) == CGAL::ZERO);
      assert(CGAL::AlgebraicSphereFunctors::contains<AK>(l,res_ls[0].first));
      assert(theSigh_at(s,res_ls[1].first) == CGAL::ZERO);
      assert(CGAL::AlgebraicSphereFunctors::contains<AK>(l,res_ls[1].first));
      if(res_ls[0].first.x() != res_ls[1].first.x()) {
        assert(res_ls[0].first.x() < res_ls[1].first.x());
      } else if(res_ls[0].first.y() != res_ls[1].first.y()) {
        assert(res_ls[0].first.y() < res_ls[1].first.y());
      } else {
        assert(res_ls[0].first.z() < res_ls[1].first.z());
      }
    }
  }

}

template <class AK>
void _test_sign_at(AK ak)
{
  typedef typename AK::Root_for_spheres_2_3 Root_for_spheres_2_3;
  typename AK::Sign_at theSigh_at =
    ak.sign_at_object();

  //Polynomial_for_spheres_2_3
  typename AK::Construct_polynomial_for_spheres_2_3 theConstruct_2_3 =
    ak.construct_polynomial_for_spheres_2_3_object();

  assert(theSigh_at(theConstruct_2_3(5, 5, 5, 100),
                    Root_for_spheres_2_3(-5,5,5)) == CGAL::ZERO);
  assert(theSigh_at(theConstruct_2_3(5, 5, 5, 100),
                    Root_for_spheres_2_3(15,5,5)) == CGAL::ZERO);
  assert(theSigh_at(theConstruct_2_3(5, 5, 5, 100),
                    Root_for_spheres_2_3(5,15,5)) == CGAL::ZERO);
  assert(theSigh_at(theConstruct_2_3(5, 5, 5, 100),
                    Root_for_spheres_2_3(5,-5,5)) == CGAL::ZERO);
  assert(theSigh_at(theConstruct_2_3(5, 5, 5, 100),
                    Root_for_spheres_2_3(5,5,-5)) == CGAL::ZERO);
  assert(theSigh_at(theConstruct_2_3(5, 5, 5, 100),
                    Root_for_spheres_2_3(5,5,15)) == CGAL::ZERO);

  assert(theSigh_at(theConstruct_2_3(5, 5, 5, 100),
                    Root_for_spheres_2_3(5,5,5)) == CGAL::NEGATIVE);
  assert(theSigh_at(theConstruct_2_3(5, 5, 5, 100),
                    Root_for_spheres_2_3(5,16,5)) == CGAL::POSITIVE);
  assert(theSigh_at(theConstruct_2_3(5, 5, 5, 100),
                    Root_for_spheres_2_3(5,-6,5)) == CGAL::POSITIVE);

  //Polynomial_1_3
  typename AK::Construct_polynomial_1_3 theConstruct_1_3 =
    ak.construct_polynomial_1_3_object();

  assert(theSigh_at(theConstruct_1_3(1, 2, 3, 6),
                    Root_for_spheres_2_3(-1,-1,-1)) == CGAL::ZERO);
  assert(theSigh_at(theConstruct_1_3(1, 2, 3, 6),
                    Root_for_spheres_2_3(2,-1,-2)) == CGAL::ZERO);
  assert(theSigh_at(theConstruct_1_3(1, 2, 3, 6),
                    Root_for_spheres_2_3(1,1,1)) != CGAL::ZERO);
  assert(theSigh_at(theConstruct_1_3(1, 2, 3, 6),
                    Root_for_spheres_2_3(0,0,0)) != CGAL::ZERO);
}
template <class AK>
void _test_critical_points_circle( const std::pair<
      typename AK::Polynomial_for_spheres_2_3,
      typename AK::Polynomial_1_3> &c) {
  typedef typename AK::FT FT;

  typedef typename AK::Root_for_spheres_2_3 Root_for_spheres_2_3;
  typename AK::X_critical_points theX_critical_points =
    AK().x_critical_points_object();
  typename AK::Y_critical_points theY_critical_points =
    AK().y_critical_points_object();
  typename AK::Z_critical_points theZ_critical_points =
    AK().z_critical_points_object();
  typename AK::Solve theSolve = AK().solve_object();
  typedef typename AK::Polynomial_1_3
      Polynomial_1_3;



  typename AK::Sign_at theSigh_at =
    AK().sign_at_object();

  Root_for_spheres_2_3 rx_min, rx_max,
                       ry_min, ry_max,
                       rz_min, rz_max;

  // X
  if(!(is_zero(c.second.b()) && is_zero(c.second.c()))) {
    rx_min = theX_critical_points(c,true);
    rx_max = theX_critical_points(c,false);

    assert(theSigh_at(c.first, rx_min) == CGAL::ZERO);
    assert(theSigh_at(c.second, rx_min) == CGAL::ZERO);
    assert(theSigh_at(c.first, rx_max) == CGAL::ZERO);
    assert(theSigh_at(c.second, rx_max) == CGAL::ZERO);

    CGAL::Interval_nt<> ix_min = CGAL::to_interval(rx_min.x());
    CGAL::Interval_nt<> ix_max = CGAL::to_interval(rx_max.x());

    FT x1_min = FT(ix_min.inf());
    FT x1_max = FT(ix_min.sup());
    FT x2_min = FT(ix_max.inf());
    FT x2_max = FT(ix_max.sup());

    if(x1_min != x1_max) {
      Polynomial_1_3 pt_xneg_min = Polynomial_1_3(1,0,0,-(x1_min));
      Polynomial_1_3 pt_xneg_max = Polynomial_1_3(1,0,0,-(x1_max));
      std::vector< std::pair<Root_for_spheres_2_3, size_t> > intersection_test_x_1;
      std::vector< std::pair<Root_for_spheres_2_3, size_t> > intersection_test_x_2;
      theSolve(c, pt_xneg_min, std::back_inserter(intersection_test_x_1));
      theSolve(c, pt_xneg_max, std::back_inserter(intersection_test_x_2));
      assert(intersection_test_x_1.size() == 0);
      assert(intersection_test_x_2.size() > 0);
    }

    if(x2_min != x2_max) {
      Polynomial_1_3 pt_xpos_min = Polynomial_1_3(1,0,0,-(x2_min));
      Polynomial_1_3 pt_xpos_max = Polynomial_1_3(1,0,0,-(x2_max));
      std::vector< std::pair<Root_for_spheres_2_3, size_t> > intersection_test_x_1;
      std::vector< std::pair<Root_for_spheres_2_3, size_t> > intersection_test_x_2;
      theSolve(c, pt_xpos_min, std::back_inserter(intersection_test_x_1));
      theSolve(c, pt_xpos_max, std::back_inserter(intersection_test_x_2));
      assert(intersection_test_x_1.size() > 0);
      assert(intersection_test_x_2.size() == 0);
    }
  }

  // Y
  if(!(is_zero(c.second.a()) && is_zero(c.second.c()))) {
    ry_min = theY_critical_points(c,true);
    ry_max = theY_critical_points(c,false);

    assert(theSigh_at(c.first, ry_min) == CGAL::ZERO);
    assert(theSigh_at(c.second, ry_min) == CGAL::ZERO);
    assert(theSigh_at(c.first, ry_max) == CGAL::ZERO);
    assert(theSigh_at(c.second, ry_max) == CGAL::ZERO);

    CGAL::Interval_nt<> iy_min = CGAL::to_interval(ry_min.y());
    CGAL::Interval_nt<> iy_max = CGAL::to_interval(ry_max.y());

    if(iy_min.inf() > iy_max.sup()) {
      std::swap(iy_min, iy_max);
    }

    FT y1_min = FT(iy_min.inf());
    FT y1_max = FT(iy_min.sup());
    FT y2_min = FT(iy_max.inf());
    FT y2_max = FT(iy_max.sup());

    if(y1_min != y1_max) {
      Polynomial_1_3 pt_yneg_min = Polynomial_1_3(0,1,0,-(y1_min));
      Polynomial_1_3 pt_yneg_max = Polynomial_1_3(0,1,0,-(y1_max));
      std::vector< std::pair<Root_for_spheres_2_3, size_t> > intersection_test_y_1;
      std::vector< std::pair<Root_for_spheres_2_3, size_t> > intersection_test_y_2;
      theSolve(c, pt_yneg_min, std::back_inserter(intersection_test_y_1));
      theSolve(c, pt_yneg_max, std::back_inserter(intersection_test_y_2));
      assert(intersection_test_y_1.size() == 0);
      assert(intersection_test_y_2.size() > 0);
    }

    if(y2_min != y2_max) {
      Polynomial_1_3 pt_ypos_min = Polynomial_1_3(0,1,0,-(y2_min));
      Polynomial_1_3 pt_ypos_max = Polynomial_1_3(0,1,0,-(y2_max));
      std::vector< std::pair<Root_for_spheres_2_3, size_t> > intersection_test_y_1;
      std::vector< std::pair<Root_for_spheres_2_3, size_t> > intersection_test_y_2;
      theSolve(c, pt_ypos_min, std::back_inserter(intersection_test_y_1));
      theSolve(c, pt_ypos_max, std::back_inserter(intersection_test_y_2));
      assert(intersection_test_y_1.size() > 0);
      assert(intersection_test_y_2.size() == 0);
    }
  }

  // Z
  if(!(is_zero(c.second.a()) && is_zero(c.second.b()))) {
    rz_min = theZ_critical_points(c,true);
    rz_max = theZ_critical_points(c,false);

    assert(theSigh_at(c.first, rz_min) == CGAL::ZERO);
    assert(theSigh_at(c.second, rz_min) == CGAL::ZERO);
    assert(theSigh_at(c.first, rz_max) == CGAL::ZERO);
    assert(theSigh_at(c.second, rz_max) == CGAL::ZERO);

    CGAL::Interval_nt<> iz_min = CGAL::to_interval(rz_min.z());
    CGAL::Interval_nt<> iz_max = CGAL::to_interval(rz_max.z());

    if(iz_min.inf() > iz_max.sup()) {
      std::swap(iz_min, iz_max);
    }

    FT z1_min = FT(iz_min.inf());
    FT z1_max = FT(iz_min.sup());
    FT z2_min = FT(iz_max.inf());
    FT z2_max = FT(iz_max.sup());

    if(z1_min != z1_max) {
      Polynomial_1_3 pt_zneg_min = Polynomial_1_3(0,0,1,-(z1_min));
      Polynomial_1_3 pt_zneg_max = Polynomial_1_3(0,0,1,-(z1_max));
      std::vector< std::pair<Root_for_spheres_2_3, size_t> > intersection_test_z_1;
      std::vector< std::pair<Root_for_spheres_2_3, size_t> > intersection_test_z_2;
      theSolve(c, pt_zneg_min, std::back_inserter(intersection_test_z_1));
      theSolve(c, pt_zneg_max, std::back_inserter(intersection_test_z_2));
      assert(intersection_test_z_1.size() == 0);
      assert(intersection_test_z_2.size() > 0);
    }

    if(z2_min != z2_max) {
      Polynomial_1_3 pt_zpos_min = Polynomial_1_3(0,0,1,-(z2_min));
      Polynomial_1_3 pt_zpos_max = Polynomial_1_3(0,0,1,-(z2_max));
      std::vector< std::pair<Root_for_spheres_2_3, size_t> > intersection_test_z_1;
      std::vector< std::pair<Root_for_spheres_2_3, size_t> > intersection_test_z_2;
      theSolve(c, pt_zpos_min, std::back_inserter(intersection_test_z_1));
      theSolve(c, pt_zpos_max, std::back_inserter(intersection_test_z_2));
      assert(intersection_test_z_1.size() > 0);
      assert(intersection_test_z_2.size() == 0);
    }
  }
}

template <class AK>
void _test_critical_points(AK ak)
{
  CGAL::Random generatorOfgenerator;
  int random_seed = generatorOfgenerator.get_int(0, 123456);
  std::cout << "random_seed = " << random_seed << std::endl;
  CGAL::Random theRandom(random_seed);
  int random_max = 127;
  int random_min = -127;

  typedef typename AK::Root_for_spheres_2_3 Root_for_spheres_2_3;
  typename AK::Construct_polynomial_for_spheres_2_3 theConstruct_2_3 =
    ak.construct_polynomial_for_spheres_2_3_object();
  typename AK::X_critical_points theX_critical_points =
    ak.x_critical_points_object();
  typename AK::Y_critical_points theY_critical_points =
    ak.y_critical_points_object();
  typename AK::Z_critical_points theZ_critical_points =
    ak.z_critical_points_object();

  for(int i = 0; i < 20; i++){
    int x = theRandom.get_int(random_min,random_max);
    int y = theRandom.get_int(random_min,random_max);
    int z = theRandom.get_int(random_min,random_max);
    int r = theRandom.get_int(1,random_max);

    assert(theX_critical_points(theConstruct_2_3(x,y,z,r*r),true)
           == Root_for_spheres_2_3(x - r, y, z));
    assert(theX_critical_points(theConstruct_2_3(x,y,z,r*r),false)
           == Root_for_spheres_2_3(x + r, y, z));
    assert(theY_critical_points(theConstruct_2_3(x,y,z,r*r),true)
           == Root_for_spheres_2_3(x, y - r, z));
    assert(theY_critical_points(theConstruct_2_3(x,y,z,r*r),false)
           == Root_for_spheres_2_3(x, y + r, z));
    assert(theZ_critical_points(theConstruct_2_3(x,y,z,r*r),true)
           == Root_for_spheres_2_3(x, y, z - r));
    assert(theZ_critical_points(theConstruct_2_3(x,y,z,r*r),false)
           == Root_for_spheres_2_3(x, y, z + r));
  }

  typedef typename AK::Polynomial_for_spheres_2_3
      Polynomial_for_spheres_2_3;
  typedef typename AK::Polynomial_1_3
      Polynomial_1_3;
  typedef std::pair<
      Polynomial_for_spheres_2_3,
      Polynomial_1_3>       Equation_Circle;

  for(int i = 0; i < 200; i++){
    int a = theRandom.get_int(-3,3);
    int b = theRandom.get_int(-3,3);
    int c = theRandom.get_int(-3,3);
    if(a == 0 && b == 0 && c == 0) continue;
    Equation_Circle ec = std::make_pair(Polynomial_for_spheres_2_3(0,0,0,1),
                                       Polynomial_1_3(a,b,c,0));
    _test_critical_points_circle<AK>(ec);
  }

}

template <class AK>
void _test_compare_Root_for_spheres(AK ak)
{
  CGAL::Random generatorOfgenerator;
  int random_seed = generatorOfgenerator.get_int(0, 123456);
  std::cout << "random_seed = " << random_seed << std::endl;
  CGAL::Random theRandom(random_seed);
  int random_max = 127;
  int random_min = -127;

  typedef typename AK::Root_for_spheres_2_3 Root_for_spheres_2_3;

  typename AK::Compare_x theCompare_x =
    ak.compare_x_object();

  typename AK::Compare_y theCompare_y =
    ak.compare_y_object();

  typename AK::Compare_z theCompare_z =
    ak.compare_z_object();

  typename AK::Compare_xy theCompare_xy =
    ak.compare_xy_object();

  typename AK::Compare_xyz theCompare_xyz =
    ak.compare_xyz_object();

  for (int i = 0; i < 100; i++){
    Root_for_spheres_2_3 r1(theRandom.get_int(random_min,random_max),
                            theRandom.get_int(random_min,random_max),
                            theRandom.get_int(random_min,random_max));
    Root_for_spheres_2_3 r2(theRandom.get_int(random_min,random_max),
                            theRandom.get_int(random_min,random_max),
                            theRandom.get_int(random_min,random_max));
    if(r1.x() > r2.x()){
      assert(theCompare_x(r1, r2) == CGAL::LARGER);
      assert(theCompare_xy(r1, r2) == CGAL::LARGER);
    }
    else if(r1.x() == r2.x()){
      assert(theCompare_x(r1, r2) == CGAL::EQUAL);
      if(r1.y() < r2.y()){
        assert(theCompare_y(r1, r2) == CGAL::SMALLER);
        assert(theCompare_xy(r1, r2) == CGAL::SMALLER);
      }
      else if(r1.y() > r2.y()){
        assert(theCompare_y(r1, r2) == CGAL::LARGER);
        assert(theCompare_xy(r1, r2) == CGAL::LARGER);
      }
      else {
        assert(theCompare_y(r1, r2) == CGAL::EQUAL);
        assert(theCompare_xy(r1, r2) == CGAL::EQUAL);
        if(r1.z() < r2.z()){
          assert(theCompare_z(r1, r2) == CGAL::SMALLER);
          assert(theCompare_xyz(r1, r2) == CGAL::SMALLER);
        }
        else if(r1.z() > r2.z()){
          assert(theCompare_z(r1, r2) == CGAL::LARGER);
          assert(theCompare_xyz(r1, r2) == CGAL::LARGER);
        } else {
          assert(theCompare_z(r1, r2) == CGAL::EQUAL);
          assert(theCompare_xyz(r1, r2) == CGAL::EQUAL);
        }
      }
    }
    else {
      assert(theCompare_x(r1, r2) == CGAL::SMALLER);
      assert(theCompare_xy(r1, r2) == CGAL::SMALLER);
    }

    if(r1.y() > r2.y())
      assert(theCompare_y(r1, r2) == CGAL::LARGER);
    else if(r1.y() < r2.y())
      assert(theCompare_y(r1, r2) == CGAL::SMALLER);
    else
      assert(theCompare_y(r1, r2) == CGAL::EQUAL);

    if(r1.z() > r2.z())
      assert(theCompare_z(r1, r2) == CGAL::LARGER);
    else if(r1.z() < r2.z())
      assert(theCompare_z(r1, r2) == CGAL::SMALLER);
    else
      assert(theCompare_z(r1, r2) == CGAL::EQUAL);
  }
}

