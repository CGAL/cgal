// ============================================================================
//
// Copyright (c) 1998 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL
// release
// of the Computational Geometry Algorithms Library (CGAL). It is
// not
// intended for general use.
//
// ----------------------------------------------------------------------------
// 
// release       :
// release_date  :
// 
// file          : test/Triangulation/include/CGAL/_test_cls_const_Del_tr..
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Mariette Yvinec
//
// coordinator   : INRIA Sophia-Antipolis Mariette.Yvinec@sophia.inria.fr
// ============================================================================

#include <CGAL/_test_cls_constrained_triangulation_2.C>

template <class Triangul>
void 
_test_cls_const_Del_triangulation(const Triangul&)
{
  //typedef Triangulation                      Cls;
  typedef typename Triangul::Geom_traits          Gt;


  typedef typename Triangul::Point                Point;
  typedef typename Triangul::Segment              Segment;
  typedef typename Triangul::Triangle             Triangle;

  typedef typename Triangul::Vertex_handle         Vertex_handle;
  typedef typename Triangul::Face_handle           Face_handle;
  typedef std::pair<Face_handle,int>          Edge;

  typedef typename Triangul::Locate_type           Locate_type;

  typedef std::pair<Point,Point>              Constraint ;
  typedef std::list<Constraint>               list_constraints;
  typedef typename list_constraints::iterator list_iterator;

  /***********************/
  /***** SUBCLASS ******/
   _test_cls_constrained_triangulation( Triangul() );

   // build triangulation T2
   list_constraints l;

   Point lpt[20] = {
     Point(0,0), Point(1,0), Point(2,0), Point(3,0),Point(4,0),
     Point(4,1), Point(3,1), Point(2,1), Point(1,1),Point(0,1),
     Point(0,2), Point(1,2), Point(2,2), Point(3,2),Point(4,2),
     Point(4,3), Point(3,3), Point(2,3), Point(1,3),Point(0,3)
   };
   for (int m=0; m<19;  m++) 
     l.push_back(Constraint(lpt[m],lpt[m+1]));
   Triangul T2(l);
   assert( T2.dimension() == 2 );
   assert( T2.number_of_vertices() == 20);
   assert( T2.is_valid() );

   // test get_conflicts
  std:: cout << "    get conflicts" << std::endl;
  std::list<Face_handle> conflicts;
  std::list<Edge>  hole_bd;
  assert(T2.get_conflicts(Point(1,1,2), 
			   std::back_inserter(conflicts)));
  conflicts.clear();	 
  assert(T2.get_conflicts_and_boundary(Point(1,1,2), 
				       std::back_inserter(conflicts),
				       std::back_inserter(hole_bd)));
  assert(hole_bd.size() == conflicts.size() + 2);
  conflicts.clear();
  hole_bd.clear();
  assert(T2.get_conflicts(Point(0,1,2), 
			  std::back_inserter(conflicts)));
  assert(T2.get_boundary_of_conflicts(Point(0,1,2), 
				      std::back_inserter(hole_bd)));
  assert(hole_bd.size() == conflicts.size() + 2);
  conflicts.clear();
  unsigned int nch = hole_bd.size();
  hole_bd.clear();
  T2.find_conflicts(Point(0,1,2), hole_bd);
  assert(hole_bd.size() == nch);
  hole_bd.clear();
  assert(!T2.get_conflicts(Point(0,0,1),std::back_inserter(conflicts)));
  conflicts.clear();
  assert(T2.get_conflicts(Point(-1,-1,1),std::back_inserter(conflicts)));
  conflicts.clear();
   
  T2.get_conflicts_and_boundary(Point(0,1,2), 
				std::back_inserter(conflicts),
				std::back_inserter(hole_bd));
  T2.star_hole(Point(0,1,2), 
	       hole_bd.begin(),
	       hole_bd.end(),
	       conflicts.begin(),
	       conflicts.end());
}
