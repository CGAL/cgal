#include "Point_set_scene_item.h"
#include "Polyhedron_type.h"
#include <CGAL/compute_normal.h>

#include <CGAL/IO/read_off_point_set.h>
#include <CGAL/IO/write_off_point_set.h>
#include <CGAL/IO/read_xyz_point_set.h>
#include <CGAL/IO/write_xyz_point_set.h>

#include <QObject>

#include <set>
#include <stack>
#include <algorithm>
#include <boost/array.hpp>


Point_set_scene_item::Point_set_scene_item()
  : Scene_item_with_display_list(),
    m_points(new Point_set)
{
  setRenderingMode(PointsPlusNormals);
}

// Copy constructor
Point_set_scene_item::Point_set_scene_item(const Point_set_scene_item& toCopy)
  : Scene_item_with_display_list(), // do not call superclass' copy constructor
    m_points(new Point_set(*toCopy.m_points))
{
  setRenderingMode(PointsPlusNormals);
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
    const Point& p = v->point();
    Vector n = compute_vertex_normal<Polyhedron::Vertex,Kernel>(*v);
    m_points->push_back(UI_point(p,n));
  }

  setRenderingMode(PointsPlusNormals);
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
bool Point_set_scene_item::read_off_point_set(std::istream& in)
{
  Q_ASSERT(m_points != NULL);

  m_points->clear();
  bool success = in &&
                 CGAL::read_off_point_set(in, std::back_inserter(*m_points)) &&
                 !isEmpty();
    
  // Mark all normals as oriented
  m_points->unoriented_points_begin() = m_points->end();
  
  return success;
}

// Write point set to .OFF file
bool Point_set_scene_item::write_off_point_set(std::ostream& out) const
{
  Q_ASSERT(m_points != NULL);

  return out &&
         CGAL::write_off_point_set(out, m_points->begin(), m_points->end());
}

// Load point set from .XYZ file
bool Point_set_scene_item::read_xyz_point_set(std::istream& in)
{
  Q_ASSERT(m_points != NULL);

  m_points->clear();
  bool success = in &&
                 CGAL::read_xyz_point_set(in, std::back_inserter(*m_points)) &&
                 !isEmpty();
    
  // Mark all normals as oriented
  m_points->unoriented_points_begin() = m_points->end();
  
  return success;
}

// Write point set to .XYZ file
bool Point_set_scene_item::write_xyz_point_set(std::ostream& out) const
{
  Q_ASSERT(m_points != NULL);

  return out &&
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

bool Point_set_scene_item::supportsRenderingMode(RenderingMode m) const {
  return m==Points || m==PointsPlusNormals || m==Splatting;
}

// Points OpenGL drawing in a display list
void Point_set_scene_item::direct_draw() const
{
  Q_ASSERT(m_points != NULL);

  // Draw points
  m_points->gl_draw_vertices();
}

// Normals OpenGL drawing
void Point_set_scene_item::draw_normals() const
{
  Q_ASSERT(m_points != NULL);

  // Draw normals
  bool points_have_normals = (m_points->begin() != m_points->end() &&
                              m_points->begin()->normal() != CGAL::NULL_VECTOR);
  if(points_have_normals)
  {
    Sphere region_of_interest = m_points->region_of_interest();
    float normal_length = (float)sqrt(region_of_interest.squared_radius() / 1000.);

    m_points->gl_draw_normals(normal_length);
  }
}

void Point_set_scene_item::draw_splats() const
{
  Q_ASSERT(m_points != NULL);

  // Draw splats
  bool points_have_normals = (m_points->begin() != m_points->end() &&
                              m_points->begin()->normal() != CGAL::NULL_VECTOR);
  bool points_have_radii =   (m_points->begin() != m_points->end() &&
                              m_points->begin()->radius() != FT(0));
  if(points_have_normals && points_have_radii)
  {
    m_points->gl_draw_splats();
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
