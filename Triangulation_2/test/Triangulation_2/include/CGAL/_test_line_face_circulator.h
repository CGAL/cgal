// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
// 
// release       :
// release_date  :
// 
// package       : Triangulation/test/Triangulation
// file          : include/CGAL/_test_cls_triangulation_2.
// source        : $URL$
// revision      : $Id$
// revision_date : $Date$

// author(s)     : Mariette Yvinec (Mariette.Yvinec@sophia.inria.fr)

// coordinator   : INRIA Sophia-Antipolis
// ============================================================================

#include <vector>

template <class Tri>
void
_test_line_face_circulator( const Tri & )
{
  //typedef Triangul                  Tri;
 typedef typename Tri::Point                Point;

 typedef typename Tri::Vertex_handle        Vertex_handle;
 typedef typename Tri::Face_handle          Face_handle;
 typedef typename Tri::Line_face_circulator Line_face_circulator;
 typedef typename Tri::Locate_type          Locate_type;

 // square - vertices of the triangulation
 Point p0(0,0,1);
 Point p1(2,0,1);
 Point p2(0,2,1);
 Point p3(2,2,1);

 // collinear outside points
 Point q1(-1,0,1);
 Point q2(0,-1,1);
 Point q3(2,-1,1);
 Point q4(3,0,1);
 Point q5(3,2,1);
 Point q6(2,3,1);
 Point q7(0,3,1);
 Point q8(-1,3,1);

 // middle of edges
 Point m1(1,0,1);
 Point m2(2,1,1);
 Point m3(1,2,1);
 Point m4(0,1,1);

 // middle of square and neighboring squares
 Point t0(1,1,1);
 Point t1(1,-1,1);
 Point t2(3,-1,1);
 Point t3(3,1,1);
 Point t4(3,3,1);
 Point t5(1,3,1);
 Point t6(-1,3,1);
 Point t7(-1,1,1);
 Point t8(-1,-1,1);

 Point zz(1,1,2);


 std::vector<Point> p;
 p.push_back(p0); p.push_back(p1); p.push_back(p2); p.push_back(p3);

 std::vector<Point> q;
 q.push_back(q1); q.push_back(q2); q.push_back(q3);
 q.push_back(q4); q.push_back(q5); q.push_back(q6);
 q.push_back(q7); q.push_back(q8);

 std::vector<Point> t;
 t.push_back(t0);
 t.push_back(t1); t.push_back(t2); t.push_back(t3);
 t.push_back(t4); t.push_back(t5); t.push_back(t6);
 t.push_back(t7); t.push_back(t8);

 std::vector<Point> m;
 m.push_back(m1); m.push_back(m2); m.push_back(m3); m.push_back(m4);

 Tri tr;
 typename std::vector<Point>::iterator pit;
 typename std::vector<Point>::iterator qit;
 typename std::vector<Point>::iterator mit;
 typename std::vector<Point>::iterator tit;
 int i; 
 Locate_type(lt);

 // insert points p - create Vertex_handle vector
 std::vector<Vertex_handle> v;
  for(pit=p.begin() ; pit != p.end(); pit++) {
    Vertex_handle vh = tr.insert(*pit);
    v.push_back(vh);
 }
 tr.is_valid();

 Face_handle f1,f2;
 assert(tr.is_face(v[0],v[1],v[2],f1));
 assert(tr.is_face(v[1],v[2],v[3],f2));

 //test line face circulator from vertex
 Line_face_circulator lfc;
 lfc = Line_face_circulator(v[0],&tr,t0);  assert(f1 == lfc);
 lfc = Line_face_circulator(v[0],&tr,m4);  assert(lfc == 0);
 lfc = Line_face_circulator(v[0],&tr,t7);  assert(lfc == 0);
 lfc = Line_face_circulator(v[0],&tr,q1);  assert(lfc == 0);
 lfc = Line_face_circulator(v[0],&tr,t8);  assert(lfc != 0);
 lfc = Line_face_circulator(v[0],&tr,q2);  assert(lfc != 0);
 lfc = Line_face_circulator(v[0],&tr,t1);  assert(lfc == 0);
 lfc = Line_face_circulator(v[0],&tr,m1);  assert(f1 == lfc);

 lfc = Line_face_circulator(v[1],&tr,p2); assert(f1 == lfc);
 lfc = Line_face_circulator(v[0],&tr,p1); assert(f1 == lfc);


 //Locate vertices, middle points and edges tested the
 //line_face_circulator creator from a vertex_handle and a point

//locate vertices
 for(pit=p.begin() ; pit != p.end(); pit++) {
   tr.locate(*pit,lt,i);
   assert(lt == Tri::VERTEX);
 }

// locate middle points
 for(mit=m.begin() ; mit != m.end(); mit++) {
   tr.locate(*mit, lt,i);
   assert(lt == Tri::EDGE);
 }

 // locate outside
 for(qit=q.begin() ; qit != q.end(); qit++) {
   tr.locate(*qit, lt,i);
   assert(lt == Tri::OUTSIDE_CONVEX_HULL);
 }

 //locate middle of squares
 tit = t.begin();
 tr.locate(*tit, lt,i);
 assert(lt == Tri::EDGE);
 ++tit;
 for(; tit != t.end(); tit++) {
   tr.locate(*tit, lt,i);
   assert(lt == Tri::OUTSIDE_CONVEX_HULL);
 }
 
 // test creator from two point
 lfc = tr.line_walk(p0,t0); assert(lfc != 0);
 lfc = tr.line_walk(p0,m4); assert(lfc == 0);
 lfc = tr.line_walk(p0,t5); assert(lfc != 0);
 lfc = tr.line_walk(p0,q1); assert(lfc == 0);
 lfc = tr.line_walk(p0,t8); assert(lfc != 0);
 lfc = tr.line_walk(p0,q2); assert(lfc != 0);
 lfc = tr.line_walk(p0,t1); assert(lfc == 0);
 lfc = tr.line_walk(p0,m1); assert(lfc != 0);

 lfc = tr.line_walk(p1,p2); assert(lfc != 0);
 assert( lfc->has_vertex(v[0]) && 
	 lfc->has_vertex(v[1]) &&
	 lfc->has_vertex(v[2]));
 lfc = tr.line_walk(p0,p1); assert(lfc != 0);
 assert (lfc->has_vertex(v[0]) && 
	 lfc->has_vertex(v[1]) &&
	 lfc->has_vertex(v[2]));

 lfc = tr.line_walk(t0,p0); assert(lfc != 0);
 lfc = tr.line_walk(m4,p0); assert(lfc != 0);
 lfc = tr.line_walk(t5,p0); assert(lfc != 0);
 lfc = tr.line_walk(q1,p0); assert(lfc != 0);
 lfc = tr.line_walk(t8,p0); assert(lfc != 0);
 lfc = tr.line_walk(q2,p0); assert(lfc == 0);
 lfc = tr.line_walk(t1,p0); assert(lfc == 0);
 lfc = tr.line_walk(m1,p0); assert(lfc == 0);

 // test creator from two points with a hint
 lfc = tr.line_walk(p0,t0,f1); assert(f1 == lfc);
 lfc = tr.line_walk(p0,m4,f1); assert(lfc == 0);
 lfc = tr.line_walk(p0,t5,f1); assert(f1 == lfc);
 lfc = tr.line_walk(p0,q1,f1); assert(lfc == 0);
 lfc = tr.line_walk(p0,t8,f1); assert(f1 == lfc);
 lfc = tr.line_walk(p0,q2,f1); assert(f1 == lfc);
 lfc = tr.line_walk(p0,t1,f1); assert(lfc == 0);
 lfc = tr.line_walk(p0,m1,f1); assert(f1 == lfc);
 lfc = tr.line_walk(t0,p0,f1); assert(f1 == lfc);
 lfc = tr.line_walk(t0,p0,f2); assert(f2 == lfc);
 lfc = tr.line_walk(m4,p0,f1); assert(f1 == lfc);
 lfc = tr.line_walk(m3,p0,f2); assert(f2 == lfc);
 lfc = tr.line_walk(m1,p0,f1); assert(lfc == 0);
 lfc = tr.line_walk(p1,p2,f2); assert(f1 == lfc);
 lfc = tr.line_walk(p1,p2,f1); assert(f1 == lfc);
 lfc = tr.line_walk(t0,p2,f1); assert(f1 == lfc);
 lfc = tr.line_walk(t0,p1,f1); assert(f2 == lfc);
 lfc = tr.line_walk(zz,p1,f1); assert(f1 == lfc);
 lfc = tr.line_walk(zz,m1,f1); assert(f1 == lfc);
 lfc = tr.line_walk(zz,m3,f1); assert(f1 == lfc);

 // a few more tests for difficult cases (collinear edges on convex_hull)
 v.push_back(tr.insert(m1));
 v.push_back(tr.insert(m2));
 v.push_back(tr.insert(m3));
 v.push_back(tr.insert(m4));
 assert(tr.is_face(v[0],v[4],v[7],f1));
 assert(tr.is_face(v[3],v[6],v[5],f2));

 lfc = tr.line_walk(m1,p1); assert(f1 == lfc);
 lfc = tr.line_walk(p1,m1); assert(lfc == 0);

 lfc = tr.line_walk(p0, p3, f1); assert(f1 == lfc);
 lfc = tr.line_walk(p0, p1, f1); assert(f1 == lfc);
 lfc = tr.line_walk(p0, p2, f1); assert(lfc == 0);
 lfc = tr.line_walk(p0, t8, f1); assert(lfc != 0);
 lfc = tr.line_walk(p0, t1, f1); assert(lfc == 0);
 return;
}
