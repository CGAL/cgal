// Copyright (c) 2017 GeometryFactory (France).
// All rights reserved.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is part of CGAL (www.cgal.org)
//
// SPDX-License-Identifier: LGPL-2.1-only
//
// Author(s)     : Mael Rouxel-Labbé

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Simple_homogeneous.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_kth_root.h>
#include <CGAL/Filtered_kernel.h>

#include <CGAL/Has_conversion.h>

#include <CGAL/internal/Exact_type_selector.h>
#include <CGAL/assertions.h>
#include <CGAL/use.h>

int main()
{
  typedef double                                                  coord_type;
  typedef CGAL::Interval_nt<true>                                 NT_approx;
  typedef CGAL::internal::Exact_field_selector<coord_type>::Type  NT_exact;

  typedef CGAL::Simple_cartesian<coord_type>                      SC;
  typedef CGAL::Simple_cartesian<NT_approx>                       ASC;
  typedef CGAL::Filtered_kernel<SC>                               FSC;

  typedef CGAL::Simple_homogeneous<NT_exact>                      SH;
  typedef CGAL::Filtered_kernel<SH>                               FSH;

  typedef CGAL::Exact_predicates_exact_constructions_kernel_with_kth_root EPECK;

  CGAL_USE_TYPE(ASC);
  CGAL_USE_TYPE(FSC);
  CGAL_USE_TYPE(SH);
  CGAL_USE_TYPE(FSH);
  CGAL_USE_TYPE(EPECK);

  CGAL_assertion((CGAL::Has_conversion<SC, SC, SC::Point_2, SC::Point_2>::value));
  CGAL_assertion((CGAL::Has_conversion<SC, SC, SC::Object_2, SC::Object_2>::value));

  CGAL_assertion(!(CGAL::Has_conversion<SC, SC, SC::Point_2, SC::Point_3>::value));
  CGAL_assertion(!(CGAL::Has_conversion<SC, SC, SC::Iso_cuboid_3, SC::Circle_2>::value));

  CGAL_assertion((CGAL::Has_conversion<SC, FSC, SC::Vector_2, FSC::Vector_2>::value));
  CGAL_assertion((CGAL::Has_conversion<FSC, SC, FSC::Vector_3, SC::Vector_3>::value));
  CGAL_assertion((CGAL::Has_conversion<SH, FSH, SH::Vector_3, FSH::Vector_3>::value));

  CGAL_assertion((CGAL::Has_conversion<SC, ASC, SC::Sphere_3, ASC::Sphere_3>::value));
  CGAL_assertion((CGAL::Has_conversion<SC, EPECK, SC::Triangle_2, EPECK::Triangle_2>::value));
  CGAL_assertion((CGAL::Has_conversion<EPECK, SC, EPECK::Circle_3, SC::Circle_3>::value));

  CGAL_assertion(!(CGAL::Has_conversion<SC, EPECK, SC::Weighted_point_2, EPECK::Weighted_point_3>::value));
  CGAL_assertion(!(CGAL::Has_conversion<SC, ASC, SC::Point_2, ASC::Weighted_point_2>::value));

  // below will produce static assert failures
//  CGAL_assertion((CGAL::Has_conversion<SC, SH, SC::Point_2, SH::Point_2>::value));
//  CGAL_assertion((CGAL::Has_conversion<FSH, EPECK, FSH::Point_3, EPECK::Point_2>::value));


}
