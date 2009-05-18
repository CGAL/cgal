#ifndef POISSON_DT3_H
#define POISSON_DT3_H

// CGAL
#include <CGAL/Reconstruction_triangulation_3.h>

// STL
#include <vector>
#include <algorithm>
#include <CGAL/gl.h>


/// The Poisson_dt3 class is a 3D triangulation class that provides:
/// - the interface requested by the Poisson_reconstruction_function class
/// - OpenGL rendering
///
/// @heading Parameters:
/// @param BaseGt   Geometric traits class.
/// @param Gt       Geometric traits class / Point_3 == Point_with_normal_3<BaseGt>.
/// @param Tds      Model of TriangulationDataStructure_3. The vertex class
///                 must derive from Reconstruction_vertex_base_3.

template <class BaseGt,
          class Gt = CGAL::Reconstruction_triangulation_default_geom_traits_3<BaseGt>,
          class Tds = CGAL::Triangulation_data_structure_3<CGAL::Reconstruction_vertex_base_3<Gt> > >
class Poisson_dt3 : public CGAL::Reconstruction_triangulation_3<BaseGt,Gt,Tds>
{
// Private types
private:

  // Base class
  typedef CGAL::Reconstruction_triangulation_3<BaseGt,Gt,Tds> Base;

// Public types
public:

  /// Geometric traits class / Point_3 == Point_with_normal_3<BaseGt>.
  typedef Gt Geom_traits;

  // Repeat Reconstruction_triangulation_3 public types
  /// @cond SKIP_IN_MANUAL
  typedef Tds Triangulation_data_structure;
  typedef typename Base::Segment      Segment;
  typedef typename Base::Triangle     Triangle;
  typedef typename Base::Tetrahedron  Tetrahedron;
  typedef typename Base::Line         Line;
  typedef typename Base::Ray          Ray;
  typedef typename Base::Object Object;
  typedef typename Base::Cell_handle   Cell_handle;
  typedef typename Base::Vertex_handle Vertex_handle;
  typedef typename Base::Cell   Cell;
  typedef typename Base::Vertex Vertex;
  typedef typename Base::Facet  Facet;
  typedef typename Base::Edge   Edge;
  typedef typename Base::Cell_circulator  Cell_circulator;
  typedef typename Base::Facet_circulator Facet_circulator;
  typedef typename Base::Cell_iterator    Cell_iterator;
  typedef typename Base::Facet_iterator   Facet_iterator;
  typedef typename Base::Edge_iterator    Edge_iterator;
  typedef typename Base::Vertex_iterator  Vertex_iterator;
  typedef typename Base::Point_iterator Point_iterator;
  typedef typename Base::Finite_vertices_iterator Finite_vertices_iterator;
  typedef typename Base::Finite_cells_iterator    Finite_cells_iterator;
  typedef typename Base::Finite_facets_iterator   Finite_facets_iterator;
  typedef typename Base::Finite_edges_iterator    Finite_edges_iterator;
  typedef typename Base::All_cells_iterator       All_cells_iterator;
  typedef typename Base::Locate_type Locate_type;
  typedef typename Base::FT FT;
  typedef typename Base::Vector Vector; ///< == Vector_3<BaseGt>
  typedef typename Base::Point Point;  ///< == Point_with_normal_3<BaseGt>
  typedef typename Base::Point_with_normal Point_with_normal; ///< Point_with_normal_3<BaseGt>
  typedef typename Base::Sphere Sphere;
  typedef typename Base::Normal_iterator Normal_iterator;
  /// @endcond

// Public methods
public:

  // Default constructor, copy constructor and operator =() are fine

  void gl_draw_delaunay_vertices(unsigned char r, unsigned char g, unsigned char b,
                                 float point_size) const
  {
    // Draw input points
    ::glPointSize(point_size);
    ::glColor3ub(r,g,b);
    ::glBegin(GL_POINTS);
    Finite_vertices_iterator v;
    for(v = finite_vertices_begin(); v != finite_vertices_end(); v++)
    {
      if(v->type() != INPUT)
        continue;
      const Point& p = v->point();
      ::glVertex3d(p.x(),p.y(),p.z());
    }
    ::glEnd();

    // Draw Steiner points
    ::glPointSize(1.0f);
    ::glColor3ub(0,0,125);
    ::glBegin(GL_POINTS);
    for(v = finite_vertices_begin(); v != finite_vertices_end(); v++)
    {
      if(v->type() != STEINER)
        continue;
      const Point& p = v->point();
      ::glVertex3d(p.x(),p.y(),p.z());
    }
    ::glEnd();
  }

  // draw Delaunay edges
  void gl_draw_delaunay_edges(unsigned char r, unsigned char g, unsigned char b,
                              float width) const
  {
    ::glLineWidth(width);
    ::glColor3ub(r,g,b);
    ::glBegin(GL_LINES);
    Finite_edges_iterator it;
    for(it = finite_edges_begin();
      it != finite_edges_end();
      it++)
    {
      Segment s = segment(*it);
      Point p1 = s.source();
      Point p2 = s.target();
      ::glVertex3d(p1.x(),p1.y(),p1.z());
      ::glVertex3d(p2.x(),p2.y(),p2.z());
    }
    ::glEnd();
  }

  // Draw normals
  void gl_draw_normals(unsigned char r, unsigned char g, unsigned char b,
                       FT c = 1.0) const
  {
    ::glColor3ub(r,g,b);
    ::glBegin(GL_LINES);
    for(Finite_vertices_iterator v = finite_vertices_begin();
        v != finite_vertices_end();
        v++)
    {
      Vector n = v->normal();
      if ( n != CGAL::NULL_VECTOR && v->type() == 0)
      {
        Point a = v->point();
        Point b = a + c * n;
        glVertex3d(a.x(),a.y(),a.z());
        glVertex3d(b.x(),b.y(),b.z());
      }
    }
    ::glEnd();
  }

}; // end of class Poisson_dt3


#endif // POISSON_DT3_H
