#ifndef GL_DRAW_DELAUNAY_H
#define GL_DRAW_DELAUNAY_H

// CGAL
#include <CGAL/Reconstruction_triangulation_3.h>

// STL
#include <CGAL/gl.h>

/// OpenGL rendering of the vertices of a Reconstruction_triangulation_3 triangulation.
///
/// @heading Parameters:
/// @param BaseGt   Geometric traits class.
/// @param Gt       Geometric traits class / Point_3 is a typedef to Point_with_normal_3<BaseGt>.
/// @param Tds      Model of TriangulationDataStructure_3. The vertex class
///                 must derive from Reconstruction_vertex_base_3.
template <class BaseGt, class Gt, class Tds>
void gl_draw_delaunay_vertices(
  const CGAL::Reconstruction_triangulation_3<BaseGt,Gt,Tds>& triangulation,
  unsigned char r, unsigned char g, unsigned char b,
  float point_size)
{
  // Triangulation
  typedef CGAL::Reconstruction_triangulation_3<BaseGt,Gt,Tds> Triangulation;

  // Geometric types
  typedef typename Triangulation::Finite_vertices_iterator Finite_vertices_iterator;
  typedef typename Triangulation::Point Point;  ///< typedef to Point_with_normal_3<BaseGt>

  // Draw input points
  ::glPointSize(point_size);
  ::glColor3ub(r,g,b);
  ::glBegin(GL_POINTS);
  Finite_vertices_iterator v;
  for(v = triangulation.finite_vertices_begin(); 
      v != triangulation.finite_vertices_end(); 
      v++)
  {
    if(v->type() != Triangulation::INPUT)
      continue;
    const Point& p = v->point();
    ::glVertex3d(p.x(),p.y(),p.z());
  }
  ::glEnd();

  // Draw Steiner points
  ::glPointSize(1.0f);
  ::glColor3ub(0,0,125);
  ::glBegin(GL_POINTS);
  for(v = triangulation.finite_vertices_begin(); 
      v != triangulation.finite_vertices_end(); 
      v++)
  {
    if(v->type() != Triangulation::STEINER)
      continue;
    const Point& p = v->point();
    ::glVertex3d(p.x(),p.y(),p.z());
  }
  ::glEnd();
}

/// OpenGL rendering of the edges of a Reconstruction_triangulation_3 triangulation.
///
/// @heading Parameters:
/// @param BaseGt   Geometric traits class.
/// @param Gt       Geometric traits class / Point_3 is a typedef to Point_with_normal_3<BaseGt>.
/// @param Tds      Model of TriangulationDataStructure_3. The vertex class
///                 must derive from Reconstruction_vertex_base_3.
template <class BaseGt, class Gt, class Tds>
void gl_draw_delaunay_edges(
  const CGAL::Reconstruction_triangulation_3<BaseGt,Gt,Tds>& triangulation,
  unsigned char r, unsigned char g, unsigned char b,
  float width)
{
  // Triangulation
  typedef CGAL::Reconstruction_triangulation_3<BaseGt,Gt,Tds> Triangulation;

  // Geometric types
  typedef typename Triangulation::Segment Segment;
  typedef typename Triangulation::Finite_edges_iterator Finite_edges_iterator;
  typedef typename Triangulation::Point Point;  ///< typedef to Point_with_normal_3<BaseGt>

  ::glLineWidth(width);
  ::glColor3ub(r,g,b);
  ::glBegin(GL_LINES);
  Finite_edges_iterator it;
  for(it = triangulation.finite_edges_begin();
      it != triangulation.finite_edges_end();
      it++)
  {
    Segment s = triangulation.segment(*it);
    Point p1 = s.source();
    Point p2 = s.target();
    ::glVertex3d(p1.x(),p1.y(),p1.z());
    ::glVertex3d(p2.x(),p2.y(),p2.z());
  }
  ::glEnd();
}

/// OpenGL rendering of the vertex normals of a Reconstruction_triangulation_3 triangulation.
///
/// @heading Parameters:
/// @param BaseGt   Geometric traits class.
/// @param Gt       Geometric traits class / Point_3 is a typedef to Point_with_normal_3<BaseGt>.
/// @param Tds      Model of TriangulationDataStructure_3. The vertex class
///                 must derive from Reconstruction_vertex_base_3.
template <class BaseGt, class Gt, class Tds>
void gl_draw_delaunay_normals(
  const CGAL::Reconstruction_triangulation_3<BaseGt,Gt,Tds>& triangulation,
  unsigned char r, unsigned char g, unsigned char b,
  FT c = 1.0)
{
  // Triangulation
  typedef CGAL::Reconstruction_triangulation_3<BaseGt,Gt,Tds> Triangulation;

  // Geometric types
  typedef typename Triangulation::Finite_vertices_iterator Finite_vertices_iterator;
  typedef typename Triangulation::Point Point;  ///< typedef to Point_with_normal_3<BaseGt>

  ::glColor3ub(r,g,b);
  ::glBegin(GL_LINES);
  Finite_vertices_iterator v;
  for(v = triangulation.finite_vertices_begin(); 
      v != triangulation.finite_vertices_end(); 
      v++)
  {
    Vector n = v->normal();
    if (n != CGAL::NULL_VECTOR && v->type() == Triangulation::INPUT)
    {
      Point a = v->point();
      Point b = a + c * n;
      glVertex3d(a.x(),a.y(),a.z());
      glVertex3d(b.x(),b.y(),b.z());
    }
  }
  ::glEnd();
}


#endif // GL_DRAW_DELAUNAY_H
