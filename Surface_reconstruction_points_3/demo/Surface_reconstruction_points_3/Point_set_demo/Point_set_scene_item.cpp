#include "Point_set_scene_item.h"

#include <CGAL/IO/read_off_point_set.h>
#include <CGAL/IO/write_off_point_set.h>
#include <CGAL/IO/read_xyz_point_set.h>
#include <CGAL/IO/write_xyz_point_set.h>

#include <QObject>
#include <QtDebug>

#include <set>
#include <stack>
#include <algorithm>
#include <boost/array.hpp>



Point_set_scene_item::Point_set_scene_item()
  : Scene_item_with_display_list()
{
}

Point_set_scene_item::Point_set_scene_item(const Point_set_scene_item& toCopy)
  // do not call superclass' copy constructor
{
  this->m_points = toCopy.m_points;
}

Point_set_scene_item::~Point_set_scene_item()
{
}

Point_set_scene_item* 
Point_set_scene_item::clone() const 
{
  Point_set_scene_item* new_item = new Point_set_scene_item(*this);
  return new_item;
}

// Load point set from .OFF file
bool 
Point_set_scene_item::read_off_point_set(std::istream& in)
{
  m_points.clear();
  return in &&
         CGAL::read_off_point_set(in, std::back_inserter(m_points)) &&
         !isEmpty();
}

// Write point set to .OFF file
bool 
Point_set_scene_item::write_off_point_set(std::ostream& out) const
{
  return out && !isEmpty() && 
         CGAL::write_off_point_set(out, m_points.begin(), m_points.end());
}

// Load point set from .XYZ file
bool 
Point_set_scene_item::read_xyz_point_set(std::istream& in)
{
  m_points.clear();
  return in &&
         CGAL::read_xyz_point_set(in, std::back_inserter(m_points)) && 
         !isEmpty();
}

// Write point set to .XYZ file
bool 
Point_set_scene_item::write_xyz_point_set(std::ostream& out) const
{
  return out && !isEmpty() && 
         CGAL::write_xyz_point_set(out, m_points.begin(), m_points.end());
}

QString 
Point_set_scene_item::toolTip() const
{
  return QObject::tr("<p><b>%1</b> (color: %4)<br />"
                     "<i>Point set</i></p>"
                     "<p>Number of vertices: %2</p>")
    .arg(name())
    .arg(m_points.size())
    .arg(color().name());
}

void
Point_set_scene_item::direct_draw() const 
{
  Sphere region_of_interest = m_points.region_of_interest();

  // Draw points
  m_points.gl_draw_vertices(color().red(),color().blue(),color().green(), 
                            2.0f /*size*/);
                            
  // Draw normals
  bool points_have_normals = (m_points.begin() != m_points.end() &&
                              m_points.begin()->normal() != CGAL::NULL_VECTOR);
  if(points_have_normals)
  {
    float normal_length = (float)sqrt(region_of_interest.squared_radius() / 10000.);
    m_points.gl_draw_normals(0,255,0 /*green*/, 
                             normal_length);
  }
}

bool
Point_set_scene_item::isEmpty() const 
{
  return m_points.empty();
}

Point_set_scene_item::Bbox
Point_set_scene_item::bbox() const 
{
  Iso_cuboid bbox = m_points.bounding_box();
  return Bbox(bbox.xmin(),bbox.ymin(),bbox.zmin(),
              bbox.xmax(),bbox.ymax(),bbox.zmax());
}

#include "Point_set_scene_item.moc"
