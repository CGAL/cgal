#include "Scene_points_with_normal_item.h"
#include "Polyhedron_type.h"
#include <CGAL/Polygon_mesh_processing/compute_normal.h>

#include <CGAL/IO/read_ply_points.h>
#include <CGAL/IO/write_ply_points.h>
#include <CGAL/IO/read_off_points.h>
#include <CGAL/IO/write_off_points.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/IO/write_xyz_points.h>
#include <CGAL/Timer.h>
#include <CGAL/Memory_sizer.h>

#include <CGAL/Three/Viewer_interface.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_3.h>

#include <QObject>
#include <QMenu>
#include <QGLViewer/manipulatedCameraFrame.h>

#include <set>
#include <stack>
#include <algorithm>
#include <boost/array.hpp>


Scene_points_with_normal_item::Scene_points_with_normal_item()
    : Scene_item(NbOfVbos,NbOfVaos),
    m_points(new Point_set),
    m_has_normals(false)
{
  setRenderingMode(Points);
    is_selected = true;
    nb_points = 0;
    nb_selected_points = 0;
    nb_lines = 0;
}

// Copy constructor
Scene_points_with_normal_item::Scene_points_with_normal_item(const Scene_points_with_normal_item& toCopy)
    :Scene_item(NbOfVbos,NbOfVaos), // do not call superclass' copy constructor
    m_points(new Point_set(*toCopy.m_points)),
    m_has_normals(toCopy.m_has_normals)
{
  if (m_has_normals)
    {
    setRenderingMode(PointsPlusNormals);
        is_selected = true;
    }
  else
    {
    setRenderingMode(Points);
        is_selected = true;
    }
    nb_points = 0;
    nb_selected_points = 0;
    nb_lines = 0;
    invalidateOpenGLBuffers();
}

// Converts polyhedron to point set
Scene_points_with_normal_item::Scene_points_with_normal_item(const Polyhedron& input_mesh)
    : Scene_item(NbOfVbos,NbOfVaos),
    m_points(new Point_set),
    m_has_normals(true)
{
  // Converts Polyhedron vertices to point set.
  // Computes vertices normal from connectivity.

  Polyhedron::Vertex_iterator v;
  for (v = const_cast<Polyhedron&>(input_mesh).vertices_begin(); 
       v != const_cast<Polyhedron&>(input_mesh).vertices_end(); v++)
  {
    const Kernel::Point_3& p = v->point();
    Kernel::Vector_3 n =
      CGAL::Polygon_mesh_processing::compute_vertex_normal(v, input_mesh);
    m_points->push_back(UI_point(p,n));
  }

  setRenderingMode(PointsPlusNormals);
  is_selected = true;
  nb_points = 0;
  nb_selected_points = 0;
  nb_lines = 0;
  invalidateOpenGLBuffers();
}

Scene_points_with_normal_item::~Scene_points_with_normal_item()
{
  Q_ASSERT(m_points != NULL);
  delete m_points; m_points = NULL;
}



void Scene_points_with_normal_item::initialize_buffers(CGAL::Three::Viewer_interface *viewer) const
{
    compute_normals_and_vertices();
    //vao for the edges
    {
        program = getShaderProgram(PROGRAM_NO_SELECTION, viewer);
        program->bind();

        vaos[Edges]->bind();
        buffers[Edges_vertices].bind();
        buffers[Edges_vertices].allocate(positions_lines.data(),
                            static_cast<int>(positions_lines.size()*sizeof(double)));
        program->enableAttributeArray("vertex");
        program->setAttributeBuffer("vertex",GL_DOUBLE,0,3);
        buffers[Edges_vertices].release();

        vaos[Edges]->release();
        nb_lines = positions_lines.size();
        positions_lines.resize(0);
        std::vector<double>(positions_lines).swap(positions_lines);
        program->release();
    }
    //vao for the points
    {
        program = getShaderProgram(PROGRAM_NO_SELECTION, viewer);
        program->bind();

        vaos[ThePoints]->bind();
        buffers[Points_vertices].bind();
        buffers[Points_vertices].allocate(positions_points.data(),
                            static_cast<int>(positions_points.size()*sizeof(double)));
        program->enableAttributeArray("vertex");
        program->setAttributeBuffer("vertex",GL_DOUBLE,0,3);
        buffers[Points_vertices].release();
        vaos[ThePoints]->release();
        nb_points = positions_points.size();
        positions_points.resize(0);
        std::vector<double>(positions_points).swap(positions_points);
        program->release();
    }
    //vao for the selected points
    {
        program = getShaderProgram(PROGRAM_NO_SELECTION, viewer);
        program->bind();

        vaos[Selected_points]->bind();
        buffers[Selected_points_vertices].bind();
        buffers[Selected_points_vertices].allocate(positions_selected_points.data(),
                            static_cast<int>(positions_selected_points.size()*sizeof(double)));
        program->enableAttributeArray("vertex");
        program->setAttributeBuffer("vertex",GL_DOUBLE,0,3);
        buffers[Selected_points_vertices].release();

        vaos[Selected_points]->release();
        nb_selected_points = positions_selected_points.size();
        positions_selected_points.resize(0);
        std::vector<double>(positions_selected_points).swap(positions_selected_points);
        program->release();
    }
    are_buffers_filled = true;
}

void Scene_points_with_normal_item::compute_normals_and_vertices() const
{
    positions_points.resize(0);
    positions_lines.resize(0);
    positions_selected_points.resize(0);
    normals.resize(0);

    positions_points.reserve(m_points->size() * 3);
    positions_lines.reserve(m_points->size() * 3 * 2);


    //Shuffle container to allow quick display random points
    Point_set_3<Kernel> points = *m_points;
    std::random_shuffle (points.begin(), points.end() - m_points->nb_selected_points());
    std::random_shuffle (points.end() - m_points->nb_selected_points(), points.end());
    
    //The points
    {
        // The *non-selected* points
      for (Point_set_3<Kernel>::const_iterator it = points.begin(); it != points.first_selected(); it++)
	{
	  const UI_point& p = *it;
	  positions_points.push_back(p.x());
	  positions_points.push_back(p.y());
	  positions_points.push_back(p.z());
	}

        // Draw *selected* points
      for (Point_set_3<Kernel>::const_iterator it = points.first_selected(); it != points.end(); it++)
	{
	  const UI_point& p = *it;
	  positions_selected_points.push_back(p.x());
	  positions_selected_points.push_back(p.y());
	  positions_selected_points.push_back(p.z());
	}
    }

    //The lines
    {
        // Stock normals
        Kernel::Sphere_3 region_of_interest = points.region_of_interest();
        float normal_length = (float)std::sqrt(region_of_interest.squared_radius() / 1000.);

	for (Point_set_3<Kernel>::const_iterator it = points.begin(); it != points.end(); it++)
	  {
	    const UI_point& p = *it;
	    const Point_set_3<Kernel>::Vector& n = p.normal();
	    Point_set_3<Kernel>::Point q = p + normal_length * n;
	    positions_lines.push_back(p.x());
	    positions_lines.push_back(p.y());
	    positions_lines.push_back(p.z());

	    positions_lines.push_back(q.x());
	    positions_lines.push_back(q.y());
	    positions_lines.push_back(q.z());
	  }
    }
}

// Duplicates scene item
Scene_points_with_normal_item*
Scene_points_with_normal_item::clone() const
{
  return new Scene_points_with_normal_item(*this);
}

// Is selection empty?
bool Scene_points_with_normal_item::isSelectionEmpty() const
{
  return (m_points->nb_selected_points() == 0);
}

// Delete selection
void Scene_points_with_normal_item::deleteSelection()
{
  CGAL::Timer task_timer; task_timer.start();
  std::cerr << "Delete " << m_points->nb_selected_points() << " points...";

  // Delete selected points
  m_points->delete_selection();

  std::size_t memory = CGAL::Memory_sizer().virtual_size();
  std::cerr << "done: " << task_timer.time() << " seconds, "
                        << (memory>>20) << " Mb allocated"
                        << std::endl;
  invalidateOpenGLBuffers();
  Q_EMIT itemChanged();
}

// Invert selection
void Scene_points_with_normal_item::invertSelection()
{
  m_points->invert_selection();
  invalidateOpenGLBuffers();
  Q_EMIT itemChanged();
}

// Select everything
void Scene_points_with_normal_item::selectAll()
{
  m_points->select_all();
  invalidateOpenGLBuffers();
  Q_EMIT itemChanged();
}
// Reset selection mark
void Scene_points_with_normal_item::resetSelection()
{
  // Un-select all points
  m_points->unselect_all();
  invalidateOpenGLBuffers();
  Q_EMIT itemChanged();
}
  //Select duplicated points
void Scene_points_with_normal_item::selectDuplicates()
{
  std::set<Kernel::Point_3> unique_points;
  for (Point_set::iterator ptit=m_points->begin(); ptit!=m_points->end();++ptit )
    if ( !unique_points.insert(*ptit).second )
      m_points->select(ptit);
  invalidateOpenGLBuffers();
  Q_EMIT itemChanged();
}

// Loads point set from .PLY file
bool Scene_points_with_normal_item::read_ply_point_set(std::istream& stream)
{
  Q_ASSERT(m_points != NULL);

  m_points->clear();
  bool ok = stream &&
            CGAL::read_ply_points_and_normals(stream,
                                              std::back_inserter(*m_points),
                                              CGAL::make_normal_of_point_with_normal_pmap(Point_set::value_type())) &&
            !isEmpty();
  if (ok)
    {
      for (Point_set::iterator it=m_points->begin(),
             end=m_points->end();it!=end; ++it)
        {
          if (it->normal() != CGAL::NULL_VECTOR)
            {
              m_has_normals=true;
              setRenderingMode(PointsPlusNormals);
              break;
            }
        }
    }
  invalidateOpenGLBuffers();
  return ok;
}

// Write point set to .PLY file
bool Scene_points_with_normal_item::write_ply_point_set(std::ostream& stream) const
{
  Q_ASSERT(m_points != NULL);

  return stream &&
         CGAL::write_ply_points_and_normals(stream,
                                            m_points->begin(), m_points->end(),
                                            CGAL::make_normal_of_point_with_normal_pmap(Point_set::value_type()));
}

// Loads point set from .OFF file
bool Scene_points_with_normal_item::read_off_point_set(std::istream& stream)
{
  Q_ASSERT(m_points != NULL);

  m_points->clear();
  bool ok = stream &&
            CGAL::read_off_points_and_normals(stream,
                                              std::back_inserter(*m_points),
                                              CGAL::make_normal_of_point_with_normal_pmap(Point_set::value_type())) &&
            !isEmpty();
  invalidateOpenGLBuffers();
  return ok;
}

// Write point set to .OFF file
bool Scene_points_with_normal_item::write_off_point_set(std::ostream& stream) const
{
  Q_ASSERT(m_points != NULL);

  return stream &&
         CGAL::write_off_points_and_normals(stream,
                                            m_points->begin(), m_points->end(),
                                            CGAL::make_normal_of_point_with_normal_pmap(Point_set::value_type()));
}

// Loads point set from .XYZ file
bool Scene_points_with_normal_item::read_xyz_point_set(std::istream& stream)
{
  Q_ASSERT(m_points != NULL);

  m_points->clear();
  bool ok = stream &&
            CGAL::read_xyz_points_and_normals(stream,
                                              std::back_inserter(*m_points),
                                              CGAL::make_normal_of_point_with_normal_pmap(Point_set::value_type())) &&
            !isEmpty();

  if (ok)
  {
    for (Point_set::iterator it=m_points->begin(),
                             end=m_points->end();it!=end; ++it)
    {
      if (it->normal() != CGAL::NULL_VECTOR)
      {
        m_has_normals=true;
        setRenderingMode(PointsPlusNormals);
        break;
      }
    }
  }
  invalidateOpenGLBuffers();
  return ok;
}

// Write point set to .XYZ file
bool Scene_points_with_normal_item::write_xyz_point_set(std::ostream& stream) const
{
  Q_ASSERT(m_points != NULL);

  return stream &&
         CGAL::write_xyz_points_and_normals(stream,
                                            m_points->begin(), m_points->end(),
                                            CGAL::make_normal_of_point_with_normal_pmap(Point_set::value_type()));
}

QString
Scene_points_with_normal_item::toolTip() const
{
  Q_ASSERT(m_points != NULL);

  return QObject::tr("<p><b>%1</b> (color: %4)<br />"
                     "<i>Point set</i></p>"
                     "<p>Number of points: %2</p>")
    .arg(name())
    .arg(m_points->size())
    .arg(color().name());
}

bool Scene_points_with_normal_item::supportsRenderingMode(RenderingMode m) const 
{
    return m==Points ||
            ( has_normals() &&
              ( m==PointsPlusNormals || m==Splatting ) );
}

void Scene_points_with_normal_item::draw_splats(CGAL::Three::Viewer_interface* viewer) const
{
   // TODO add support for selection
   viewer->glBegin(GL_POINTS);
   for ( Point_set_3<Kernel>::const_iterator it = m_points->begin(); it != m_points->end(); it++)
   {
     const UI_point& p = *it;
     viewer->glNormal3dv(&p.normal().x());
     viewer->glMultiTexCoord1d(GL_TEXTURE2, p.radius());
     viewer->glVertex3dv(&p.x());

   }
   viewer->glEnd();



}

void Scene_points_with_normal_item::draw_edges(CGAL::Three::Viewer_interface* viewer) const
{
    double ratio_displayed = 1.0;
    if (viewer->inFastDrawing () &&
        (nb_lines/6 > 300000)) // arbitrary large value
      ratio_displayed = 6 * 300000. / (double)(nb_lines);

    if(!are_buffers_filled)
        initialize_buffers(viewer);
    vaos[Edges]->bind();
    program=getShaderProgram(PROGRAM_NO_SELECTION);
    attrib_buffers(viewer,PROGRAM_NO_SELECTION);
    program->bind();
    program->setAttributeValue("colors", this->color());
    viewer->glDrawArrays(GL_LINES, 0,
                         static_cast<GLsizei>(((std::size_t)(ratio_displayed * nb_lines)/3)));
    vaos[Edges]->release();
    program->release();
}
void Scene_points_with_normal_item::draw_points(CGAL::Three::Viewer_interface* viewer) const
{
    if(!are_buffers_filled)
        initialize_buffers(viewer);

    double ratio_displayed = 1.0;
    if (viewer->inFastDrawing () &&
        ((nb_points + nb_selected_points)/3 > 300000)) // arbitrary large value
      ratio_displayed = 3 * 300000. / (double)(nb_points + nb_selected_points);

    vaos[ThePoints]->bind();
    program=getShaderProgram(PROGRAM_NO_SELECTION);
    attrib_buffers(viewer,PROGRAM_NO_SELECTION);
    program->bind();
    program->setAttributeValue("colors", this->color());
    viewer->glDrawArrays(GL_POINTS, 0,
                         static_cast<GLsizei>(((std::size_t)(ratio_displayed * nb_points)/3)));
    vaos[ThePoints]->release();
    program->release();
    GLfloat point_size;
    viewer->glGetFloatv(GL_POINT_SIZE, &point_size);
    viewer->glPointSize(4.f);

    vaos[Selected_points]->bind();
    program=getShaderProgram(PROGRAM_NO_SELECTION);
    attrib_buffers(viewer,PROGRAM_NO_SELECTION);
    program->bind();
    program->setAttributeValue("colors", QColor(255,0,0));
    viewer->glDrawArrays(GL_POINTS, 0,
                         static_cast<GLsizei>(((std::size_t)(ratio_displayed * nb_selected_points)/3)));
    vaos[Selected_points]->release();
    program->release();
    viewer->glPointSize(point_size);
}
// Gets wrapped point set
Point_set* Scene_points_with_normal_item::point_set()
{
  Q_ASSERT(m_points != NULL);
  return m_points;
}
const Point_set* Scene_points_with_normal_item::point_set() const
{
  Q_ASSERT(m_points != NULL);
  return m_points;
}

bool
Scene_points_with_normal_item::isEmpty() const
{
  Q_ASSERT(m_points != NULL);
  return m_points->empty();
}

void
Scene_points_with_normal_item::compute_bbox() const
{
  Q_ASSERT(m_points != NULL);

  Kernel::Iso_cuboid_3 bbox = m_points->bounding_box();
  _bbox = Bbox(bbox.xmin(),bbox.ymin(),bbox.zmin(),
              bbox.xmax(),bbox.ymax(),bbox.zmax());
}

void Scene_points_with_normal_item::computes_local_spacing(int k)
{
  typedef Kernel Geom_traits;
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

QMenu* Scene_points_with_normal_item::contextMenu()
{
    const char* prop_name = "Menu modified by Scene_points_with_normal_item.";

    QMenu* menu = Scene_item::contextMenu();

    // Use dynamic properties:
    // http://doc.qt.io/qt-5/qobject.html#property
    bool menuChanged = menu->property(prop_name).toBool();

    if(!menuChanged) {
        actionDeleteSelection = menu->addAction(tr("Delete Selection"));
        actionDeleteSelection->setObjectName("actionDeleteSelection");
        connect(actionDeleteSelection, SIGNAL(triggered()),this, SLOT(deleteSelection()));

        actionResetSelection = menu->addAction(tr("Reset Selection"));
        actionResetSelection->setObjectName("actionResetSelection");
        connect(actionResetSelection, SIGNAL(triggered()),this, SLOT(resetSelection()));

        actionSelectDuplicatedPoints = menu->addAction(tr("Select duplicated points"));
        actionSelectDuplicatedPoints->setObjectName("actionSelectDuplicatedPoints");
        connect(actionSelectDuplicatedPoints, SIGNAL(triggered()),this, SLOT(selectDuplicates()));

        menu->setProperty(prop_name, true);
    }

    if (isSelectionEmpty())
    {
        actionDeleteSelection->setDisabled(true);
        actionResetSelection->setDisabled(true);
    }
    else
    {
        actionDeleteSelection->setDisabled(false);
        actionResetSelection->setDisabled(false);
    }

    return menu;
}

void Scene_points_with_normal_item::setRenderingMode(RenderingMode m)
{
    Scene_item::setRenderingMode(m);
    if (rendering_mode==Splatting && (!m_points->are_radii_uptodate()))
    {
        computes_local_spacing(6); // default value = small
    }
}

bool Scene_points_with_normal_item::has_normals() const { return m_has_normals; }

void Scene_points_with_normal_item::set_has_normals(bool b) {
  if (b!=m_has_normals){
    m_has_normals=b;
    //reset the context menu
    delete this->defaultContextMenu;
    this->defaultContextMenu = 0;
  }
}

void Scene_points_with_normal_item::invalidateOpenGLBuffers()
{
    are_buffers_filled = false;
    compute_bbox();
}

