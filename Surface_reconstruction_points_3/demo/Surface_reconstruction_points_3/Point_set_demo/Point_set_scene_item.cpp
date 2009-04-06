#include "Point_set_scene_item.h"
#include "Point_set_demo_types.h"
#include <CGAL/IO/read_off_point_set.h>
#include <CGAL/IO/write_off_point_set.h>
#include <CGAL/IO/read_xyz_point_set.h>
#include <CGAL/IO/write_xyz_point_set.h>

#include <QObject>

#include <set>
#include <stack>
#include <algorithm>
#include <boost/array.hpp>


// ----------------------------------------------------------------------------
// Private functions
// ----------------------------------------------------------------------------

// compute polyhedron facet's normal
static Vector compute_facet_normal(Polyhedron::Facet_const_handle f)
{
  Vector sum = CGAL::NULL_VECTOR;
  Polyhedron::Halfedge_around_facet_const_circulator h = f->facet_begin();
  do
  {
    Vector normal = CGAL::cross_product(
      h->next()->vertex()->point() - h->vertex()->point(),
      h->next()->next()->vertex()->point() - h->next()->vertex()->point());
    double sqnorm = normal * normal;
    if(sqnorm != 0)
      normal = normal / (double)std::sqrt(sqnorm);
    sum = sum + normal;
  }
  while(++h != f->facet_begin());
  
  // Normalize 'sum'
  double sqnorm = sum * sum; // dot product
  if(sqnorm != 0.0)
    sum = sum / std::sqrt(sqnorm);

  return sum;
}

// compute polyhedron vertex's normal from connectivity
static Vector compute_vertex_normal(Polyhedron::Vertex_const_handle v)
{
  Vector normal = CGAL::NULL_VECTOR;
  Polyhedron::Halfedge_around_vertex_const_circulator pHalfedge = v->vertex_begin(),
                                                      end = pHalfedge;
  CGAL_For_all(pHalfedge,end)
    if(!pHalfedge->is_border())
      normal = normal + compute_facet_normal(pHalfedge->facet());
      
  // Normalize 'normal'
  double sqnorm = normal * normal;
  if(sqnorm != 0.0)
    normal = normal / std::sqrt(sqnorm);
  
  return normal;
}


// ----------------------------------------------------------------------------
// Class Point_set_scene_item
// ----------------------------------------------------------------------------

Point_set_scene_item::Point_set_scene_item()
  : Scene_item_with_display_list(),
    m_points(new Point_set)
{
}

// Copy constructor
Point_set_scene_item::Point_set_scene_item(const Point_set_scene_item& toCopy)
  : Scene_item_with_display_list(), // do not call superclass' copy constructor
    m_points(new Point_set(*toCopy.m_points))
{
}

// Convert polyhedron to point set
Point_set_scene_item::Point_set_scene_item(const Polyhedron& input_mesh)
  : Scene_item_with_display_list(),
    m_points(new Point_set)
{
  // Convert Polyhedron vertices to point set.
  // Compute vertices' normals from connectivity.
  Polyhedron::Vertex_const_iterator v;
  for (v = input_mesh.vertices_begin(); v != input_mesh.vertices_end(); v++)
  {
    Point p = v->point();
    Vector n = compute_vertex_normal(v);
    m_points->push_back(UI_point(p,n));
  }

}

Point_set_scene_item::~Point_set_scene_item()
{
  Q_ASSERT(m_points != NULL);
  delete m_points; m_points = NULL;
}

// Duplicate scene item
Point_set_scene_item* 
Point_set_scene_item::clone() const 
{
  return new Point_set_scene_item(*this);
}

// Load point set from .OFF file
bool 
Point_set_scene_item::read_off_point_set(std::istream& in)
{
  Q_ASSERT(m_points != NULL);
  
  m_points->clear();
  return in &&
         CGAL::read_off_point_set(in, std::back_inserter(*m_points)) &&
         !isEmpty();
}

// Write point set to .OFF file
bool 
Point_set_scene_item::write_off_point_set(std::ostream& out) const
{
  Q_ASSERT(m_points != NULL);

  return out && 
         !isEmpty() && 
         CGAL::write_off_point_set(out, m_points->begin(), m_points->end());
}

// Load point set from .XYZ file
bool 
Point_set_scene_item::read_xyz_point_set(std::istream& in)
{
  Q_ASSERT(m_points != NULL);
  
  m_points->clear();
  return in &&
         CGAL::read_xyz_point_set(in, std::back_inserter(*m_points)) && 
         !isEmpty();
}

// Write point set to .XYZ file
bool 
Point_set_scene_item::write_xyz_point_set(std::ostream& out) const
{
  Q_ASSERT(m_points != NULL);
  
  return out && 
         !isEmpty() && 
         CGAL::write_xyz_point_set(out, m_points->begin(), m_points->end());
}

QString 
Point_set_scene_item::toolTip() const
{
  Q_ASSERT(m_points != NULL);

  return QObject::tr("<p><b>%1</b> (color: %4)<br />"
                     "<i>Point set</i></p>"
                     "<p>Number of vertices: %2</p>")
    .arg(name())
    .arg(m_points->size())
    .arg(color().name());
}

void
Point_set_scene_item::direct_draw() const 
{
  Q_ASSERT(m_points != NULL);

  Sphere region_of_interest = m_points->region_of_interest();

  // Draw points
  m_points->gl_draw_vertices(color().red(),color().blue(),color().green(), 
                             2.0f /*size*/);
                            
  // Draw normals
  bool points_have_normals = (m_points->begin() != m_points->end() &&
                              m_points->begin()->normal() != CGAL::NULL_VECTOR);
  if(points_have_normals)
  {
    float normal_length = (float)sqrt(region_of_interest.squared_radius() / 10000.);
    m_points->gl_draw_normals(0,255,0 /*green*/, 
                              normal_length);
  }
}

// Get wrapped point set
Point_set* Point_set_scene_item::point_set() 
{ 
  Q_ASSERT(m_points != NULL);
  return m_points; 
}
const Point_set* Point_set_scene_item::point_set() const
{ 
  Q_ASSERT(m_points != NULL);
  return m_points; 
}

bool
Point_set_scene_item::isEmpty() const 
{
  Q_ASSERT(m_points != NULL);
  return m_points->empty();
}

Point_set_scene_item::Bbox
Point_set_scene_item::bbox() const 
{
  Q_ASSERT(m_points != NULL);

  Iso_cuboid bbox = m_points->bounding_box();
  return Bbox(bbox.xmin(),bbox.ymin(),bbox.zmin(),
              bbox.xmax(),bbox.ymax(),bbox.zmax());
}

#include "Point_set_scene_item.moc"
