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
#include <CGAL/compute_average_spacing.h>

#include <CGAL/Three/Viewer_interface.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_3.h>

#include <QObject>
#include <QApplication>
#include <QMenu>
#include <QSlider>
#include <QWidgetAction>
#include <QGLViewer/manipulatedCameraFrame.h>

#include <set>
#include <stack>
#include <algorithm>
#include <boost/array.hpp>

struct Scene_points_with_normal_item_priv
{
  Scene_points_with_normal_item_priv(Scene_points_with_normal_item* parent)
    :m_points(new Point_set),
      m_has_normals(false)
  {
    item = parent;
    nb_points = 0;
    nb_selected_points = 0;
    nb_lines = 0;
    normal_Slider = new QSlider(Qt::Horizontal);
    normal_Slider->setValue(20);
    point_Slider = new QSlider(Qt::Horizontal);
    point_Slider->setValue(5);
    point_Slider->setMinimum(1);
    point_Slider->setMaximum(25);
  }
  Scene_points_with_normal_item_priv(const Scene_points_with_normal_item& toCopy, Scene_points_with_normal_item* parent)
    : m_points(new Point_set(*toCopy.d->m_points)),
      m_has_normals(toCopy.d->m_has_normals)
  {
    item = parent;
    normal_Slider = new QSlider(Qt::Horizontal);
    normal_Slider->setValue(20);
    point_Slider = new QSlider(Qt::Horizontal);
    point_Slider->setValue(5);
    point_Slider->setMinimum(1);
    point_Slider->setMaximum(25);
  }
  Scene_points_with_normal_item_priv(const Polyhedron& input_mesh, Scene_points_with_normal_item* parent)
    : m_points(new Point_set),
      m_has_normals(true)
  {
    item = parent;
    nb_points = 0;
    nb_selected_points = 0;
    nb_lines = 0;
    Polyhedron::Vertex_iterator v;
    for (v = const_cast<Polyhedron&>(input_mesh).vertices_begin();
         v != const_cast<Polyhedron&>(input_mesh).vertices_end(); v++)
    {
      const Kernel::Point_3& p = v->point();
      Kernel::Vector_3 n =
        CGAL::Polygon_mesh_processing::compute_vertex_normal(v, input_mesh);
      m_points->push_back(UI_point(p,n));
    }
    normal_Slider = new QSlider(Qt::Horizontal);
    normal_Slider->setValue(20);
    point_Slider = new QSlider(Qt::Horizontal);
    point_Slider->setValue(5);
    point_Slider->setMinimum(1);
    point_Slider->setMaximum(25);
  }
  ~Scene_points_with_normal_item_priv()
  {
    Q_ASSERT(m_points != NULL);
    delete m_points; m_points = NULL;
    delete normal_Slider;
    delete point_Slider;
  }
  void initializeBuffers(CGAL::Three::Viewer_interface *viewer) const;
  void compute_normals_and_vertices() const;
  enum VAOs {
      Edges=0,
      ThePoints,
      Selected_points,
      NbOfVaos
  };
  enum VBOs {
      Edges_vertices = 0,
      Points_vertices,
      Points_normals,
      Selected_points_vertices,
      Selected_points_normals,
      NbOfVbos
  };
  Point_set* m_points;
  bool m_has_normals;
  QAction* actionDeleteSelection;
  QAction* actionResetSelection;
  QAction* actionSelectDuplicatedPoints;
  QSlider* normal_Slider;
  QSlider* point_Slider;
  mutable std::vector<double> positions_lines;
  mutable std::vector<double> positions_points;
  mutable std::vector<double> positions_selected_points;
  mutable std::vector<double> normals;
  mutable std::vector<double> positions_normals;
  mutable std::vector<double> positions_selected_normals;
  mutable std::size_t nb_points;
  mutable std::size_t nb_selected_points;
  mutable std::size_t nb_lines;
  mutable QOpenGLShaderProgram *program;

  Scene_points_with_normal_item* item;
};


Scene_points_with_normal_item::Scene_points_with_normal_item()
    : Scene_item(Scene_points_with_normal_item_priv::NbOfVbos,Scene_points_with_normal_item_priv::NbOfVaos)
{
    setRenderingMode(Points);
    is_selected = true;
    d = new Scene_points_with_normal_item_priv(this);
}

// Copy constructor
Scene_points_with_normal_item::Scene_points_with_normal_item(const Scene_points_with_normal_item& toCopy)
    :Scene_item(Scene_points_with_normal_item_priv::NbOfVbos,Scene_points_with_normal_item_priv::NbOfVaos)
{

  d = new Scene_points_with_normal_item_priv(toCopy, this);
  if (d->m_has_normals)
    {
        setRenderingMode(PointsPlusNormals);
        is_selected = true;
    }
  else
    {
        setRenderingMode(Points);
        is_selected = true;
    }

    invalidateOpenGLBuffers();
}

// Converts polyhedron to point set
Scene_points_with_normal_item::Scene_points_with_normal_item(const Polyhedron& input_mesh)
    : Scene_item(Scene_points_with_normal_item_priv::NbOfVbos,Scene_points_with_normal_item_priv::NbOfVaos)
{
  // Converts Polyhedron vertices to point set.
  // Computes vertices normal from connectivity.
  d = new Scene_points_with_normal_item_priv(input_mesh, this);
  setRenderingMode(PointsPlusNormals);
  is_selected = true;
  invalidateOpenGLBuffers();
}

Scene_points_with_normal_item::~Scene_points_with_normal_item()
{
  delete d;
}



void Scene_points_with_normal_item_priv::initializeBuffers(CGAL::Three::Viewer_interface *viewer) const
{
    compute_normals_and_vertices();
    //vao for the edges
    {
        program = item->getShaderProgram(Scene_points_with_normal_item::PROGRAM_NO_SELECTION, viewer);
        program->bind();

        item->vaos[Edges]->bind();
        item->buffers[Edges_vertices].bind();
        item->buffers[Edges_vertices].allocate(positions_lines.data(),
                            static_cast<int>(positions_lines.size()*sizeof(double)));
        program->enableAttributeArray("vertex");
        program->setAttributeBuffer("vertex",GL_DOUBLE,0,3);
        item->buffers[Edges_vertices].release();

        item->vaos[Edges]->release();
        nb_lines = positions_lines.size();
        positions_lines.resize(0);
        std::vector<double>(positions_lines).swap(positions_lines);
        program->release();
    }
    //vao for the points
    {
        if(item->has_normals())
          program = item->getShaderProgram(Scene_points_with_normal_item::PROGRAM_WITH_LIGHT, viewer);
        else
          program = item->getShaderProgram(Scene_points_with_normal_item::PROGRAM_NO_SELECTION, viewer);
        program->bind();

        item->vaos[ThePoints]->bind();
        item->buffers[Points_vertices].bind();
        item->buffers[Points_vertices].allocate(positions_points.data(),
                            static_cast<int>(positions_points.size()*sizeof(double)));
        program->enableAttributeArray("vertex");
        program->setAttributeBuffer("vertex",GL_DOUBLE,0,3);
        item->buffers[Points_vertices].release();

        if(item->has_normals())
        {
          item->buffers[Points_normals].bind();
          item->buffers[Points_normals].allocate(positions_normals.data(),
                                          static_cast<int>(positions_normals.size()*sizeof(double)));
          program->enableAttributeArray("normals");
          program->setAttributeBuffer("normals",GL_DOUBLE,0,3);
          item->buffers[Points_normals].release();
          positions_normals.resize(0);
          std::vector<double>(positions_normals).swap(positions_normals);
        }
        item->vaos[ThePoints]->release();
        nb_points = positions_points.size();
        positions_points.resize(0);
        std::vector<double>(positions_points).swap(positions_points);
        program->release();
    }
    //vao for the selected points
    {
        if(item->has_normals())
          program = item->getShaderProgram(Scene_points_with_normal_item::PROGRAM_WITH_LIGHT, viewer);
        else
          program = item->getShaderProgram(Scene_points_with_normal_item::PROGRAM_NO_SELECTION, viewer);

        program->bind();

        item->vaos[Selected_points]->bind();
        item->buffers[Selected_points_vertices].bind();
        item->buffers[Selected_points_vertices].allocate(positions_selected_points.data(),
                            static_cast<int>(positions_selected_points.size()*sizeof(double)));
        program->enableAttributeArray("vertex");
        program->setAttributeBuffer("vertex",GL_DOUBLE,0,3);
        item->buffers[Selected_points_vertices].release();

        if(item->has_normals())
        {
          item->buffers[Selected_points_normals].bind();
          item->buffers[Selected_points_normals].allocate(positions_selected_normals.data(),
                                          static_cast<int>(positions_selected_normals.size()*sizeof(double)));
          program->enableAttributeArray("normals");
          program->setAttributeBuffer("normals",GL_DOUBLE,0,3);
          item->buffers[Selected_points_normals].release();
          positions_selected_normals.resize(0);
          std::vector<double>(positions_selected_normals).swap(positions_selected_normals);
        }
        item->vaos[Selected_points]->release();
        nb_selected_points = positions_selected_points.size();
        positions_selected_points.resize(0);
        std::vector<double>(positions_selected_points).swap(positions_selected_points);
        program->release();
    }
    item->are_buffers_filled = true;
}

void Scene_points_with_normal_item_priv::compute_normals_and_vertices() const
{
    QApplication::setOverrideCursor(Qt::WaitCursor);
    positions_points.resize(0);
    positions_lines.resize(0);
    positions_selected_points.resize(0);
    normals.resize(0);
    positions_normals.resize(0);
    positions_selected_normals.resize(0);
    normals.resize(0);

    positions_points.reserve(m_points->size() * 3);
    positions_lines.reserve(m_points->size() * 3 * 2);
    if(item->has_normals())
    {
      positions_normals.reserve((m_points->size() - m_points->nb_selected_points()) * 3);
      positions_selected_normals.reserve(m_points->nb_selected_points() * 3);
    }
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

#ifdef LINK_WITH_TBB
       typedef CGAL::Parallel_tag Concurrency_tag;
#else
        typedef CGAL::Sequential_tag Concurrency_tag;
#endif

        double average_spacing = CGAL::compute_average_spacing<Concurrency_tag>(
              points.begin(), points.end(),
              6);

        double normal_length = (std::min)(average_spacing, std::sqrt(region_of_interest.squared_radius() / 1000.));
        double length_factor = 5.0/100*normal_Slider->value();
        for (Point_set_3<Kernel>::const_iterator it = points.begin(); it != points.first_selected(); it++)
	  {
	    const UI_point& p = *it;
	    const Point_set_3<Kernel>::Vector& n = p.normal();
            Point_set_3<Kernel>::Point q = p + normal_length * length_factor* n;
	    positions_lines.push_back(p.x());
	    positions_lines.push_back(p.y());
	    positions_lines.push_back(p.z());

	    positions_lines.push_back(q.x());
	    positions_lines.push_back(q.y());
	    positions_lines.push_back(q.z());


            positions_normals.push_back(n.x());
            positions_normals.push_back(n.y());
            positions_normals.push_back(n.z());
	  }
        for (Point_set_3<Kernel>::const_iterator it = points.first_selected(); it != points.end(); it++)
          {
            const UI_point& p = *it;
            const Point_set_3<Kernel>::Vector& n = p.normal();
            Point_set_3<Kernel>::Point q = p + normal_length * length_factor* n;
            positions_lines.push_back(p.x());
            positions_lines.push_back(p.y());
            positions_lines.push_back(p.z());

            positions_lines.push_back(q.x());
            positions_lines.push_back(q.y());
            positions_lines.push_back(q.z());


            positions_selected_normals.push_back(n.x());
            positions_selected_normals.push_back(n.y());
            positions_selected_normals.push_back(n.z());
          }
    }
    QApplication::restoreOverrideCursor();
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
  return (d->m_points->nb_selected_points() == 0);
}

// Delete selection
void Scene_points_with_normal_item::deleteSelection()
{
  CGAL::Timer task_timer; task_timer.start();
  std::cerr << "Delete " << d->m_points->nb_selected_points() << " points...";

  // Delete selected points
  d->m_points->delete_selection();

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
  d->m_points->invert_selection();
  invalidateOpenGLBuffers();
  Q_EMIT itemChanged();
}

// Select everything
void Scene_points_with_normal_item::selectAll()
{
  d->m_points->select_all();
  invalidateOpenGLBuffers();
  Q_EMIT itemChanged();
}
// Reset selection mark
void Scene_points_with_normal_item::resetSelection()
{
  // Un-select all points
  d->m_points->unselect_all();
  invalidateOpenGLBuffers();
  Q_EMIT itemChanged();
}
  //Select duplicated points
void Scene_points_with_normal_item::selectDuplicates()
{
  std::set<Kernel::Point_3> unique_points;
  for (Point_set::iterator ptit=d->m_points->begin(); ptit!=d->m_points->end();++ptit )
    if ( !unique_points.insert(*ptit).second )
      d->m_points->select(ptit);
  invalidateOpenGLBuffers();
  Q_EMIT itemChanged();
}

// Loads point set from .PLY file
bool Scene_points_with_normal_item::read_ply_point_set(std::istream& stream)
{
  Q_ASSERT(d->m_points != NULL);

  d->m_points->clear();
  bool ok = stream &&
            CGAL::read_ply_points_and_normals(stream,
                                              std::back_inserter(*d->m_points),
                                              CGAL::make_normal_of_point_with_normal_pmap(Point_set::value_type())) &&
            !isEmpty();
  if (ok)
    {
      for (Point_set::iterator it=d->m_points->begin(),
             end=d->m_points->end();it!=end; ++it)
        {
          if (it->normal() != CGAL::NULL_VECTOR)
            {
              d->m_has_normals=true;
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
  Q_ASSERT(d->m_points != NULL);

  return stream &&
         CGAL::write_ply_points_and_normals(stream,
                                            d->m_points->begin(), d->m_points->end(),
                                            CGAL::make_normal_of_point_with_normal_pmap(Point_set::value_type()));
}

// Loads point set from .OFF file
bool Scene_points_with_normal_item::read_off_point_set(std::istream& stream)
{
  Q_ASSERT(d->m_points != NULL);

  d->m_points->clear();
  bool ok = stream &&
            CGAL::read_off_points_and_normals(stream,
                                              std::back_inserter(*d->m_points),
                                              CGAL::make_normal_of_point_with_normal_pmap(Point_set::value_type())) &&
            !isEmpty();
  invalidateOpenGLBuffers();
  return ok;
}

// Write point set to .OFF file
bool Scene_points_with_normal_item::write_off_point_set(std::ostream& stream) const
{
  Q_ASSERT(d->m_points != NULL);

  return stream &&
         CGAL::write_off_points_and_normals(stream,
                                            d->m_points->begin(), d->m_points->end(),
                                            CGAL::make_normal_of_point_with_normal_pmap(Point_set::value_type()));
}

// Loads point set from .XYZ file
bool Scene_points_with_normal_item::read_xyz_point_set(std::istream& stream)
{
  Q_ASSERT(d->m_points != NULL);

  d->m_points->clear();
  bool ok = stream &&
            CGAL::read_xyz_points_and_normals(stream,
                                              std::back_inserter(*d->m_points),
                                              CGAL::make_normal_of_point_with_normal_pmap(Point_set::value_type())) &&
            !isEmpty();

  if (ok)
  {
    for (Point_set::iterator it=d->m_points->begin(),
                             end=d->m_points->end();it!=end; ++it)
    {
      if (it->normal() != CGAL::NULL_VECTOR)
      {
        d->m_has_normals=true;
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
  Q_ASSERT(d->m_points != NULL);

  return stream &&
         CGAL::write_xyz_points_and_normals(stream,
                                            d->m_points->begin(), d->m_points->end(),
                                            CGAL::make_normal_of_point_with_normal_pmap(Point_set::value_type()));
}

QString
Scene_points_with_normal_item::toolTip() const
{
  Q_ASSERT(d->m_points != NULL);

  return QObject::tr("<p><b>%1</b> (color: %4)<br />"
                     "<i>Point set</i></p>"
                     "<p>Number of points: %2</p>")
    .arg(name())
    .arg(d->m_points->size())
    .arg(color().name());
}

bool Scene_points_with_normal_item::supportsRenderingMode(RenderingMode m) const 
{
    return m==Points ||
            ( has_normals() &&
              ( m==PointsPlusNormals || m==Splatting ) );
}

void Scene_points_with_normal_item::drawSplats(CGAL::Three::Viewer_interface* viewer) const
{
   // TODO add support for selection
   viewer->glBegin(GL_POINTS);
   for ( Point_set_3<Kernel>::const_iterator it = d->m_points->begin(); it != d->m_points->end(); it++)
   {
     const UI_point& p = *it;
     viewer->glNormal3dv(&p.normal().x());
     viewer->glMultiTexCoord1d(GL_TEXTURE2, p.radius());
     viewer->glVertex3dv(&p.x());

   }
   viewer->glEnd();

}

void Scene_points_with_normal_item::drawEdges(CGAL::Three::Viewer_interface* viewer) const
{
    double ratio_displayed = 1.0;
    if (viewer->inFastDrawing () &&
        (d->nb_lines/6 > 300000)) // arbitrary large value
      ratio_displayed = 6 * 300000. / (double)(d->nb_lines);

    if(!are_buffers_filled)
        d->initializeBuffers(viewer);
    vaos[Scene_points_with_normal_item_priv::Edges]->bind();
    d->program=getShaderProgram(PROGRAM_NO_SELECTION);
    attribBuffers(viewer,PROGRAM_NO_SELECTION);
    d->program->bind();
    d->program->setAttributeValue("colors", this->color());
    viewer->glDrawArrays(GL_LINES, 0,
                         static_cast<GLsizei>(((std::size_t)(ratio_displayed * d->nb_lines)/3)));
    vaos[Scene_points_with_normal_item_priv::Edges]->release();
    d->program->release();
}
void Scene_points_with_normal_item::drawPoints(CGAL::Three::Viewer_interface* viewer) const
{
    if(!are_buffers_filled)
        d->initializeBuffers(viewer);

    GLfloat point_size;
    viewer->glGetFloatv(GL_POINT_SIZE, &point_size);
    viewer->glPointSize(d->point_Slider->value());
    double ratio_displayed = 1.0;
    if (viewer->inFastDrawing () &&
        ((d->nb_points + d->nb_selected_points)/3 > 300000)) // arbitrary large value
      ratio_displayed = 3 * 300000. / (double)(d->nb_points + d->nb_selected_points);

    vaos[Scene_points_with_normal_item_priv::ThePoints]->bind();
    if(has_normals())
    {
      d->program=getShaderProgram(PROGRAM_WITH_LIGHT);
      attribBuffers(viewer,PROGRAM_WITH_LIGHT);
    }
    else
    {
      d->program=getShaderProgram(PROGRAM_NO_SELECTION);
      attribBuffers(viewer,PROGRAM_NO_SELECTION);
    }
    d->program->bind();
    d->program->setAttributeValue("colors", this->color());
    viewer->glDrawArrays(GL_POINTS, 0,
                         static_cast<GLsizei>(((std::size_t)(ratio_displayed * d->nb_points)/3)));
    vaos[Scene_points_with_normal_item_priv::ThePoints]->release();
    d->program->release();

    vaos[Scene_points_with_normal_item_priv::Selected_points]->bind();
    if(has_normals())
    {
      d->program=getShaderProgram(PROGRAM_WITH_LIGHT);
      attribBuffers(viewer,PROGRAM_WITH_LIGHT);
    }
    else
    {
      d->program=getShaderProgram(PROGRAM_NO_SELECTION);
      attribBuffers(viewer,PROGRAM_NO_SELECTION);
    }
    d->program->bind();
    d->program->setAttributeValue("colors", QColor(255,0,0));
    viewer->glDrawArrays(GL_POINTS, 0,
                         static_cast<GLsizei>(((std::size_t)(ratio_displayed * d->nb_selected_points)/3)));
    vaos[Scene_points_with_normal_item_priv::Selected_points]->release();
    d->program->release();
    viewer->glPointSize(point_size);
}
// Gets wrapped point set
Point_set* Scene_points_with_normal_item::point_set()
{
  Q_ASSERT(d->m_points != NULL);
  return d->m_points;
}
const Point_set* Scene_points_with_normal_item::point_set() const
{
  Q_ASSERT(d->m_points != NULL);
  return d->m_points;
}

bool
Scene_points_with_normal_item::isEmpty() const
{
  Q_ASSERT(d->m_points != NULL);
  return d->m_points->empty();
}

void
Scene_points_with_normal_item::compute_bbox() const
{
  Q_ASSERT(d->m_points != NULL);

  Kernel::Iso_cuboid_3 bbox = d->m_points->bounding_box();
  _bbox = Bbox(bbox.xmin(),bbox.ymin(),bbox.zmin(),
              bbox.xmax(),bbox.ymax(),bbox.zmax());
}

void Scene_points_with_normal_item::computes_local_spacing(int k)
{
  typedef Kernel Geom_traits;
  typedef CGAL::Search_traits_3<Geom_traits> TreeTraits;
  typedef CGAL::Orthogonal_k_neighbor_search<TreeTraits> Neighbor_search;
  typedef Neighbor_search::Tree Tree;

  Point_set::iterator end(d->m_points->end());

  // build kdtree
  Tree tree(d->m_points->begin(), end);

  // Compute the radius of each point = (distance max to k nearest neighbors)/2.
  {
    int i=0;
    for (Point_set::iterator it=d->m_points->begin(); it!=end; ++it, ++i)
    {
      Neighbor_search search(tree, *it, k+1);
      double maxdist2 = (--search.end())->second; // squared distance to furthest neighbor
      it->radius() = sqrt(maxdist2)/2.;
    }
  }

  d->m_points->set_radii_uptodate(true);
}

QMenu* Scene_points_with_normal_item::contextMenu()
{
    const char* prop_name = "Menu modified by Scene_points_with_normal_item.";

    QMenu* menu = Scene_item::contextMenu();

    //add a slider to modify the normals length
    // Use dynamic properties:
    // http://doc.qt.io/qt-5/qobject.html#property
    bool menuChanged = menu->property(prop_name).toBool();

    if(!menuChanged) {
      if(has_normals())
      {
        QMenu *container = new QMenu(tr("Normals Length"));
        QWidgetAction *sliderAction = new QWidgetAction(0);
        connect(d->normal_Slider, &QSlider::valueChanged, this, &Scene_points_with_normal_item::invalidateOpenGLBuffers);
        connect(d->normal_Slider, &QSlider::valueChanged, this, &Scene_points_with_normal_item::itemChanged);

        sliderAction->setDefaultWidget(d->normal_Slider);

        container->addAction(sliderAction);
        menu->addMenu(container);
      }
        QMenu *container = new QMenu(tr("Points Size"));
        QWidgetAction *sliderAction = new QWidgetAction(0);
        connect(d->point_Slider, &QSlider::valueChanged, this, &Scene_points_with_normal_item::invalidateOpenGLBuffers);
        connect(d->point_Slider, &QSlider::valueChanged, this, &Scene_points_with_normal_item::itemChanged);

        sliderAction->setDefaultWidget(d->point_Slider);

        container->addAction(sliderAction);
        menu->addMenu(container);

        d->actionDeleteSelection = menu->addAction(tr("Delete Selection"));
        d->actionDeleteSelection->setObjectName("actionDeleteSelection");
        connect(d->actionDeleteSelection, SIGNAL(triggered()),this, SLOT(deleteSelection()));

        d->actionResetSelection = menu->addAction(tr("Reset Selection"));
        d->actionResetSelection->setObjectName("actionResetSelection");
        connect(d->actionResetSelection, SIGNAL(triggered()),this, SLOT(resetSelection()));

        d->actionSelectDuplicatedPoints = menu->addAction(tr("Select duplicated points"));
        d->actionSelectDuplicatedPoints->setObjectName("actionSelectDuplicatedPoints");
        connect(d->actionSelectDuplicatedPoints, SIGNAL(triggered()),this, SLOT(selectDuplicates()));

        menu->setProperty(prop_name, true);
    }

    if (isSelectionEmpty())
    {
        d->actionDeleteSelection->setDisabled(true);
        d->actionResetSelection->setDisabled(true);
    }
    else
    {
        d->actionDeleteSelection->setDisabled(false);
        d->actionResetSelection->setDisabled(false);
    }

    return menu;
}

void Scene_points_with_normal_item::setRenderingMode(RenderingMode m)
{
    Scene_item::setRenderingMode(m);
    if (rendering_mode==Splatting && (!d->m_points->are_radii_uptodate()))
    {
        computes_local_spacing(6); // default value = small
    }
}

bool Scene_points_with_normal_item::has_normals() const { return d->m_has_normals; }

void Scene_points_with_normal_item::set_has_normals(bool b) {
  if (b!=d->m_has_normals){
    d->m_has_normals=b;
    //reset the context menu
    if (defaultContextMenu)
      defaultContextMenu->deleteLater();
    this->defaultContextMenu = 0;
  }
}

void Scene_points_with_normal_item::invalidateOpenGLBuffers()
{
    are_buffers_filled = false;
    compute_bbox();
}

