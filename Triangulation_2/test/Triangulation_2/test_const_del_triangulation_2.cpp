// ============================================================================
//
// Copyright (c) 1998 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : test/Triangulation/test_constrained_triangulation.C
// source        : $URL$
// revision      : $Id$
// revision_date : $Date$
// author(s)     : Francois Rebufat (Francois.Rebufat@sophia.inria.fr)
//               : Mariette Yvinec (Mariette.Yvinec@sophia.inria.fr)
// coordinator   : INRIA Sophia-Antipolis
// ============================================================================

#include <CGAL/_test_types.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/_test_cls_const_Del_triangulation_2.h>

// Explicit instantiation of the whole class :
template class CGAL::Constrained_Delaunay_triangulation_2<TestK>;

int main()
{
  std::cout << "Testing constrained_Delaunay_triangulation "<< std::endl;
  std::cout << " with No_constraint_intersection_requiring_constructions_tag : " << std::endl;
  typedef CGAL::Constrained_Delaunay_triangulation_2<TestK>        CDt2;
  typedef CGAL::Constrained_Delaunay_triangulation_2<CGAL::Epeck>  EPECK_CDt2;

  _test_cls_const_Del_triangulation(CDt2());
  _test_cls_const_Del_triangulation(EPECK_CDt2());

  //Testing insertion of a range of constraints
  std::vector<TestK::Point_2> points;
  points.push_back( TestK::Point_2(0,0) );
  points.push_back( TestK::Point_2(0,1) );
  points.push_back( TestK::Point_2(2,1) );
  points.push_back( TestK::Point_2(2,3) );
  {
  std::vector< std::pair<std::size_t,std::size_t> > csts;
  csts.push_back( std::make_pair(0,1) );
  csts.push_back( std::make_pair(1,2) );
  CDt2 cdt;
  cdt.insert_constraints(points.begin(), points.end(),
                         csts.begin(), csts.end() );
  }
  {
  std::vector< TestK::Segment_2 > csts;
  csts.push_back( TestK::Segment_2(points[0],points[1]) );
  csts.push_back( TestK::Segment_2(points[1],points[2]) );
  CDt2 cdt;
  cdt.insert_constraints(csts.begin(), csts.end() );
  }
  {
  std::vector< CDt2::Constraint > csts;
  csts.push_back( CDt2::Constraint(points[0],points[1]) );
  csts.push_back( CDt2::Constraint(points[1],points[2]) );
  CDt2 cdt;
  cdt.insert_constraints(csts.begin(), csts.end() );
  }
  return 0;
}
