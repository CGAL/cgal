#include "Scene_points_with_normal_item.h"
#include "Polyhedron_type.h"
#include <CGAL/Polygon_mesh_processing/compute_normal.h>

#include <CGAL/IO/read_ply_points.h>
#include <CGAL/IO/read_ply_point_set_3.h>
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
#include <CGAL/Search_traits_adapter.h>

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
    :m_points(new Point_set)
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
    : m_points(new Point_set(*toCopy.d->m_points))
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
    : m_points(new Point_set)
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
      m_points->push_back(p,n);
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
      Points_colors,
      Selected_points_vertices,
      Selected_points_normals,
      NbOfVbos
  };
  Point_set* m_points;
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
  mutable std::vector<double> colors_points;
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
  if (has_normals())
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

        if (!(colors_points.empty()))
          {
            item->buffers[Points_colors].bind();
            item->buffers[Points_colors].allocate (colors_points.data(),
                                                   static_cast<int>(colors_points.size()*sizeof(double)));
            program->enableAttributeArray("colors");
            program->setAttributeBuffer("colors",GL_DOUBLE,0,3);
            item->buffers[Points_colors].release();
          }

        item->vaos[Edges]->release();
        nb_lines = positions_lines.size();
        positions_lines.resize(0);
        std::vector<double>(positions_lines).swap(positions_lines);
        program->release();
    }
    //vao for the points
    {
        if(item->has_normals() && !(m_points->has_colors()))
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

        if(item->has_normals() && !(m_points->has_colors()))
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
        
        if (!(colors_points.empty()))
          {
            item->buffers[Points_colors].bind();
            item->buffers[Points_colors].allocate (colors_points.data(),
                                                   static_cast<int>(colors_points.size()*sizeof(double)));
            program->enableAttributeArray("colors");
            program->setAttributeBuffer("colors",GL_DOUBLE,0,3,6*sizeof(double));
            item->buffers[Points_colors].release();
            colors_points.resize(0);
            std::vector<double>(colors_points).swap(colors_points);
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
    colors_points.resize(0);

    positions_points.reserve(m_points->size() * 3);
    positions_lines.reserve(m_points->size() * 3 * 2);

    //Shuffle container to allow quick display random points
    std::random_shuffle (m_points->begin(), m_points->first_selected());
    if (m_points->nb_selected_points() != 0)
      std::random_shuffle (m_points->first_selected(), m_points->end());
    
    //The points
    {
        // The *non-selected* points
      std::size_t i = 0;
      for (; i < m_points->size () - m_points->nb_selected_points(); ++ i)
	{
	  positions_points.push_back(m_points->point(i).x());
	  positions_points.push_back(m_points->point(i).y());
	  positions_points.push_back(m_points->point(i).z());
	}

        // Draw *selected* points
      for (; i < m_points->size (); ++ i)
	{
	  positions_selected_points.push_back(m_points->point(i).x());
	  positions_selected_points.push_back(m_points->point(i).y());
	  positions_selected_points.push_back(m_points->point(i).z());
	}
    }

    //The lines
    if (item->has_normals())
    {
        positions_normals.reserve((m_points->size() - m_points->nb_selected_points()) * 3);
        positions_selected_normals.reserve(m_points->nb_selected_points() * 3);
        
        // Store normals
        Kernel::Sphere_3 region_of_interest = m_points->region_of_interest();

#ifdef LINK_WITH_TBB
       typedef CGAL::Parallel_tag Concurrency_tag;
#else
        typedef CGAL::Sequential_tag Concurrency_tag;
#endif

        double average_spacing = CGAL::compute_average_spacing<Concurrency_tag>(
              m_points->begin(), m_points->end(), m_points->point_pmap(),
              6);

        double normal_length = (std::min)(average_spacing, std::sqrt(region_of_interest.squared_radius() / 1000.));
        double length_factor = 5.0/100*normal_Slider->value();
        for (Point_set_3<Kernel>::const_iterator it = m_points->begin(); it != m_points->first_selected(); it++)
	  {
	    const Kernel::Point_3& p = m_points->point(it);
	    const Kernel::Vector_3& n = m_points->normal(it);
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
        for (Point_set_3<Kernel>::const_iterator it = m_points->first_selected(); it != m_points->end(); it++)
          {
	    const Kernel::Point_3& p = m_points->point(it);
	    const Kernel::Vector_3& n = m_points->normal(it);
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
    //The colors
    if (m_points->has_colors())
    {
        colors_points.reserve((m_points->size() - m_points->nb_selected_points()) * 3);

        for (Point_set_3<Kernel>::const_iterator it = m_points->begin(); it != m_points->end(); it++)
	  {
            colors_points.push_back ((double)(m_points->red(it) / 255.));
            colors_points.push_back ((double)(m_points->green(it) / 255.));
            colors_points.push_back ((double)(m_points->blue(it) / 255.));
            colors_points.push_back ((double)(m_points->red(it) / 255.));
            colors_points.push_back ((double)(m_points->green(it) / 255.));
            colors_points.push_back ((double)(m_points->blue(it) / 255.));
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
  for (Point_set::iterator ptit = d->m_points->begin(); ptit!= d->m_points->end();++ptit )
    if ( !unique_points.insert(d->m_points->point(ptit)).second)
      d->m_points->select(ptit);
  invalidateOpenGLBuffers();
  Q_EMIT itemChanged();
}

// Loads point set from .PLY file
bool Scene_points_with_normal_item::read_ply_point_set(std::istream& stream)
{
  Q_ASSERT(d->m_points != NULL);

  d->m_points->clear();

  CGAL::Ply_interpreter_point_set_3<Kernel> interpreter (*(d->m_points));
  
  bool ok = stream &&
            CGAL::read_ply_custom_points (stream,
                                          interpreter,
                                          Kernel()) &&
            !isEmpty();

  std::cerr << d->m_points->properties();

  if (d->m_points->has_normals())
    setRenderingMode(PointsPlusNormals);
  if (d->m_points->check_colors())
    std::cerr << "-> Point set has colors" << std::endl;
  
  invalidateOpenGLBuffers();
  return ok;
}

// Write point set to .PLY file
bool Scene_points_with_normal_item::write_ply_point_set(std::ostream& stream) const
{
  Q_ASSERT(d->m_points != NULL);

  if (d->m_points->has_normals())
    return stream &&
      CGAL::write_ply_points_and_normals(stream,
                                         d->m_points->begin(), d->m_points->end(),
                                         d->m_points->point_pmap(), d->m_points->normal_pmap(),
                                         Kernel());

  return stream &&
    CGAL::write_ply_points(stream,
                           d->m_points->begin(), d->m_points->end(),
                           d->m_points->point_pmap(),
                           Kernel());
}

// Loads point set from .OFF file
bool Scene_points_with_normal_item::read_off_point_set(std::istream& stream)
{
  Q_ASSERT(d->m_points != NULL);

  d->m_points->clear();
  bool ok = stream &&
            CGAL::read_off_points_and_normals(stream,
                                              d->m_points->index_back_inserter(),
                                              d->m_points->point_push_pmap(),
                                              d->m_points->normal_push_pmap(),
                                              Kernel()) &&
            !isEmpty();
  invalidateOpenGLBuffers();
  return ok;
}

// Write point set to .OFF file
bool Scene_points_with_normal_item::write_off_point_set(std::ostream& stream) const
{
  Q_ASSERT(d->m_points != NULL);

  if (d->m_points->has_normals())
    return stream &&
      CGAL::write_off_points_and_normals(stream,
                                         d->m_points->begin(), d->m_points->end(),
                                         d->m_points->point_pmap(), d->m_points->normal_pmap(),
                                         Kernel());
  return stream &&
    CGAL::write_off_points (stream,
                            d->m_points->begin(), d->m_points->end(),
                            d->m_points->point_pmap(),
                            Kernel());
}

// Loads point set from .XYZ file
bool Scene_points_with_normal_item::read_xyz_point_set(std::istream& stream)
{
  Q_ASSERT(d->m_points != NULL);

  d->m_points->clear();
  d->m_points->add_normal_property();
  bool ok = stream &&
            CGAL::read_xyz_points_and_normals(stream,
                                              d->m_points->index_back_inserter(),
                                              d->m_points->point_push_pmap(),
                                              d->m_points->normal_push_pmap(),
                                              Kernel()) &&
            !isEmpty();

  if (ok)
  {
    bool has_normals = false;
    for (Point_set::iterator it=d->m_points->begin(),
           end=d->m_points->end();it!=end; ++it)
      {
        if (d->m_points->normal(*it) != CGAL::NULL_VECTOR)
          {
            has_normals=true;
            setRenderingMode(PointsPlusNormals);
            break;
          }
      }
    if (!has_normals)
      d->m_points->remove_normal_property();
  }
  invalidateOpenGLBuffers();
  return ok;
}

// Write point set to .XYZ file
bool Scene_points_with_normal_item::write_xyz_point_set(std::ostream& stream) const
{
  Q_ASSERT(d->m_points != NULL);

  if (d->m_points->has_normals())
    return stream &&
      CGAL::write_xyz_points_and_normals(stream,
                                         d->m_points->begin(), d->m_points->end(),
                                         d->m_points->point_pmap(),
                                         d->m_points->normal_pmap(),
                                         Kernel());
  return stream &&
    CGAL::write_xyz_points (stream,
                            d->m_points->begin(), d->m_points->end(),
                            d->m_points->point_pmap(),
                            Kernel());
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
   if (d->m_points->has_colors())
     for ( Point_set_3<Kernel>::const_iterator it = d->m_points->begin(); it != d->m_points->end(); it++)
       {
         const Point_set::Point& p = d->m_points->point (it);
         const Point_set::Vector& n = d->m_points->normal (it);
         viewer->glColor4d((double)(d->m_points->red(it)) / 255.,
                           (double)(d->m_points->green(it)) / 255.,
                           (double)(d->m_points->blue(it)) / 255.,
                           1.0);
         viewer->glNormal3dv(&n.x());
         viewer->glMultiTexCoord1d(GL_TEXTURE2, d->m_points->radius(*it));
         viewer->glVertex3dv(&p.x());

       }
   else
     for ( Point_set_3<Kernel>::const_iterator it = d->m_points->begin(); it != d->m_points->end(); it++)
       {
         const Point_set::Point& p = d->m_points->point (it);
         const Point_set::Vector& n = d->m_points->normal (it);
         viewer->glNormal3dv(&n.x());
         viewer->glMultiTexCoord1d(GL_TEXTURE2, d->m_points->radius(*it));
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
    if(has_normals() && !(d->m_points->has_colors()))
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
    if (!(d->m_points->has_colors()))
      d->program->setAttributeValue("colors", this->color());
    viewer->glDrawArrays(GL_POINTS, 0,
                         static_cast<GLsizei>(((std::size_t)(ratio_displayed * d->nb_points)/3)));
    vaos[Scene_points_with_normal_item_priv::ThePoints]->release();
    d->program->release();

    vaos[Scene_points_with_normal_item_priv::Selected_points]->bind();
    if(has_normals() && !(d->m_points->has_colors()))
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

  typedef CGAL::Search_traits_3<Geom_traits> SearchTraits_3;
  typedef CGAL::Search_traits_adapter <std::size_t, Point_set::Point_pmap, SearchTraits_3> Search_traits;
  typedef CGAL::Orthogonal_k_neighbor_search<Search_traits> Neighbor_search;
  typedef Neighbor_search::Tree Tree;
  typedef Neighbor_search::Distance Distance;

  // build kdtree
  Tree tree(d->m_points->begin(),
            d->m_points->end(),
            typename Tree::Splitter(),
            Search_traits (d->m_points->point_pmap())
            );
  Distance tr_dist(d->m_points->point_pmap());

  if (!(d->m_points->has_property<double> ("radius")))
    d->m_points->add_radius();

  // Compute the radius of each point = (distance max to k nearest neighbors)/2.
  {
    d->m_points->test();
    int i=0;
    for (Point_set::iterator it=d->m_points->begin(); it!=d->m_points->end(); ++it, ++i)
    {
      Neighbor_search search(tree, d->m_points->point(it), k+1, 0, true, tr_dist);
      double maxdist2 = (--search.end())->second; // squared distance to furthest neighbor
      d->m_points->radius(it) = sqrt(maxdist2)/2.;
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

bool Scene_points_with_normal_item::has_normals() const { return d->m_points->has_normals(); }

void Scene_points_with_normal_item::invalidateOpenGLBuffers()
{
    are_buffers_filled = false;
    compute_bbox();
}

