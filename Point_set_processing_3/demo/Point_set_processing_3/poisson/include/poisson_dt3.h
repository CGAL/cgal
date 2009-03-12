#ifndef POISSON_DT3_H
#define POISSON_DT3_H

// CGAL
#include <CGAL/Reconstruction_triangulation_3.h>

// STL
#include <vector>
#include <algorithm>
#include <GL/gl.h>


/// The Poisson_dt3 class is a 3D triangulation class that provides:
/// - the interface requested by the Poisson_reconstruction_function class
/// - OpenGL rendering
///
/// @heading Is Model for the Concepts: 
/// Model of the ReconstructionTriangulation_3 concept.
///
/// @heading Parameters:
/// @param BaseGt   Kernel's geometric traits.
/// @param Gt       Geometric traits class / Point_3 is a model of PointWithNormal_3.
/// @param Tds      Model of TriangulationDataStructure_3. The cell base class must be
/// a model of ReconstructionCellBase_3 and the vertex base class
/// must be a model of ReconstructionVertexBase_3.

template <class BaseGt,
          class Gt = CGAL::Reconstruction_triangulation_default_geom_traits_3<BaseGt>,
          class Tds = CGAL::Triangulation_data_structure_3<CGAL::Reconstruction_vertex_base_3<Gt>,
                                                           CGAL::Reconstruction_cell_base_3<Gt> > >
class Poisson_dt3 : public CGAL::Reconstruction_triangulation_3<BaseGt,Gt,Tds>
{
// Private types
private:

  // Base class
  typedef CGAL::Reconstruction_triangulation_3<BaseGt,Gt,Tds> Base;

// Public types
public:

  // Repeat Reconstruction_triangulation_3 public types
  /// @cond SKIP_IN_MANUAL
  typedef Tds Triangulation_data_structure;
  typedef Gt  Geom_traits; ///< Geometric traits class / Point_3 is a model of PointWithNormal_3.
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
  typedef typename Base::Vector Vector;
  typedef typename Base::Iso_cuboid Iso_cuboid;
  typedef typename Base::Sphere Sphere;
  typedef typename Base::Normal_iterator Normal_iterator;
  /// @endcond
  
  /// The geometric traits class's Point_3 type is a model of PointWithNormal_3
  typedef typename Geom_traits::Point_3 Point;             ///< Model of PointWithNormal_3
  typedef typename Geom_traits::Point_3 Point_with_normal; ///< Model of PointWithNormal_3
  typedef typename Point_with_normal::Normal Normal; ///< Model of Kernel::Vector_3 concept.

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

  // Contour a triangular surface via a Marching Cubes style algorithm.
  // OutputIterator value_type must be Triangle.
  template <typename OutputIterator> 
  unsigned int marching_tet(OutputIterator output, ///< output triangles
                            const FT value)
  {
    unsigned int nb_tri = 0;
    for(Finite_cells_iterator c = finite_cells_begin();
        c != finite_cells_end();
        c++)
      nb_tri += contour(c,output,value);
    return nb_tri;
  }

private:

  // Contour a triangular surface in 'cell' via a Marching Cubes style algorithm.
  // OutputIterator value_type must be Triangle.
  template <typename OutputIterator> 
  unsigned int contour(Cell_handle cell,
                       OutputIterator& output, // notice the '&'
                       const FT value)
  {
    std::vector<Point> points;
    Point point;
    if(level_set(cell,value,0,1,point)) points.push_back(point);
    if(level_set(cell,value,0,2,point)) points.push_back(point);
    if(level_set(cell,value,0,3,point)) points.push_back(point);
    if(level_set(cell,value,1,2,point)) points.push_back(point);
    if(level_set(cell,value,1,3,point)) points.push_back(point);
    if(level_set(cell,value,2,3,point)) points.push_back(point);

    // only 3 or 4
    if(points.size() == 3)
    {
      std::vector<Point>::iterator it = points.begin();
      const Point& a = (*it); it++;
      const Point& b = (*it); it++;
      const Point& c = (*it);

      Triangle triangle = Triangle(a,b,c);
      *output++ = triangle;
      return 1;
    }
    else if(points.size() == 4)
    {
      std::vector<Point>::iterator it = points.begin();
      std::vector<Point> p(4);
      for(int i=0;i<4;i++)
      {
        p[i] = (*it);
        it++;
      }
      *output++ = Triangle(p[0],p[1],p[3]);
      *output++ = Triangle(p[0],p[3],p[2]);

      return 2;
    }
    return 0;
  }

  bool level_set(Cell_handle c,
                 const FT value,
                 const int i1,
                 const int i2,
                 Point& p)
  {
    const Point& p1 = c->vertex(i1)->point();
    const Point& p2 = c->vertex(i2)->point();
    double v1 = c->vertex(i1)->f();
    double v2 = c->vertex(i2)->f();

    if(v1 <= value && v2 >= value)
    {
      double ratio = (value - v1) / (v2 - v1);
      p = p1 + ratio * (p2-p1);
      return true;
    }
    else if(v2 <= value && v1 >= value)
    {
      double ratio = (value - v2) / (v1 - v2);
      p = p2 + ratio * (p1-p2);
      return true;
    }
    return false;
  }

public:

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
      Normal n = v->normal();
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
