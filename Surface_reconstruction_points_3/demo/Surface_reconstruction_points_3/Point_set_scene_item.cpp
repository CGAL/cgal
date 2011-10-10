#include "Point_set_scene_item.h"
#include "Polyhedron_type.h"
#include "CGAL/compute_normal.h"

#include <CGAL/IO/read_off_points.h>
#include <CGAL/IO/write_off_points.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/IO/write_xyz_points.h>
#include <CGAL/Timer.h>
#include <CGAL/Memory_sizer.h>

#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_3.h>

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

// Converts polyhedron to point set
Point_set_scene_item::Point_set_scene_item(const Polyhedron& input_mesh)
  : Scene_item_with_display_list(),
    m_points(new Point_set)
{
  // Converts Polyhedron vertices to point set.
  // Computes vertices normal from connectivity.
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

// Duplicates scene item
Point_set_scene_item*
Point_set_scene_item::clone() const
{
  return new Point_set_scene_item(*this);
}

// Is selection empty?
bool Point_set_scene_item::isSelectionEmpty() const
{
  return (m_points->nb_selected_points() == 0);
}

// Delete selection
void Point_set_scene_item::deleteSelection()
{
  CGAL::Timer task_timer; task_timer.start();
  std::cerr << "Delete " << m_points->nb_selected_points() << " points...";

  // Delete selected points
  m_points->delete_selection();

  long memory = CGAL::Memory_sizer().virtual_size();
  std::cerr << "done: " << task_timer.time() << " seconds, "
                        << (memory>>20) << " Mb allocated"
                        << std::endl;
}

// Reset selection mark
void Point_set_scene_item::resetSelection()
{
  // Un-select all points
  m_points->select(m_points->begin(), m_points->end(), false);
}

// Loads point set from .OFF file
bool Point_set_scene_item::read_off_point_set(std::istream& stream)
{
  Q_ASSERT(m_points != NULL);

  m_points->clear();
  bool ok = stream &&
            CGAL::read_off_points_and_normals(stream,
                                              std::back_inserter(*m_points),
                                              CGAL::make_normal_of_point_with_normal_pmap(std::back_inserter(*m_points))) &&
            !isEmpty();

  return ok;
}

// Write point set to .OFF file
bool Point_set_scene_item::write_off_point_set(std::ostream& stream) const
{
  Q_ASSERT(m_points != NULL);

  return stream &&
         CGAL::write_off_points_and_normals(stream,
                                            m_points->begin(), m_points->end(),
                                            CGAL::make_normal_of_point_with_normal_pmap(m_points->begin()));
}

// Loads point set from .XYZ file
bool Point_set_scene_item::read_xyz_point_set(std::istream& stream)
{
  Q_ASSERT(m_points != NULL);

  m_points->clear();
  bool ok = stream &&
            CGAL::read_xyz_points_and_normals(stream,
                                              std::back_inserter(*m_points),
                                              CGAL::make_normal_of_point_with_normal_pmap(std::back_inserter(*m_points))) &&
            !isEmpty();

  return ok;
}

// Write point set to .XYZ file
bool Point_set_scene_item::write_xyz_point_set(std::ostream& stream) const
{
  Q_ASSERT(m_points != NULL);

  return stream &&
         CGAL::write_xyz_points_and_normals(stream,
                                            m_points->begin(), m_points->end(),
                                            CGAL::make_normal_of_point_with_normal_pmap(m_points->begin()));
}

QString
Point_set_scene_item::toolTip() const
{
  Q_ASSERT(m_points != NULL);

  return QObject::tr("<p><b>%1</b> (color: %4)<br />"
                     "<i>Point set</i></p>"
                     "<p>Number of points: %2</p>")
    .arg(name())
    .arg(m_points->size())
    .arg(color().name());
}

bool Point_set_scene_item::supportsRenderingMode(RenderingMode m) const 
{
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
    float normal_length = (float)std::sqrt(region_of_interest.squared_radius() / 1000.);

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

// Gets wrapped point set
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

void Point_set_scene_item::computes_local_spacing(int k)
{
  typedef Kernel Geom_traits;
  typedef Geom_traits::FT FT;
  typedef CGAL::Search_traits_3<Geom_traits> TreeTraits;
  typedef CGAL::Orthogonal_k_neighbor_search<TreeTraits> Neighbor_search;
  typedef Neighbor_search::Tree Tree;

  Point_set::iterator end(m_points->end());

  // build kdtree
  Tree tree(m_points->begin(), end);

  // Compute the radius of each point = (distance max to k nearest neighbors)/2.
  {
    int i=0;
    for (Point_set::iterator it=m_points->begin(); it!=end; ++it, ++i)
    {
      Neighbor_search search(tree, *it, k+1);
      double maxdist2 = (--search.end())->second; // squared distance to furthest neighbor
      it->radius() = sqrt(maxdist2)/2.;
    }
  }

  m_points->set_radii_uptodate(true);
}

void Point_set_scene_item::setRenderingMode(RenderingMode m)
{
  Scene_item_with_display_list::setRenderingMode(m);
  if (rendering_mode==Splatting && (!m_points->are_radii_uptodate()))
  {
    computes_local_spacing(6); // default value = small
  }
}

#include "Point_set_scene_item.moc"
