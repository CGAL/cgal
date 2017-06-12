#define CGAL_data_type float
#define CGAL_GL_data_type GL_FLOAT
#include "Scene_points_with_normal_item.h"
#include <CGAL/Polygon_mesh_processing/compute_normal.h>

#include <CGAL/Point_set_3/IO.h>
#include <CGAL/Timer.h>
#include <CGAL/Memory_sizer.h>
#include <CGAL/compute_average_spacing.h>

#include <CGAL/Three/Viewer_interface.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/linear_least_squares_fitting_3.h>

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

#include <CGAL/boost/graph/properties_Surface_mesh.h>
#include "Polyhedron_type.h"


const std::size_t limit_fast_drawing = 300000; //arbitraty large value

struct Scene_points_with_normal_item_priv
{
  void init_values(Scene_points_with_normal_item* parent)
  {
    item = parent;
    nb_points = 0;
    nb_selected_points = 0;
    nb_lines = 0;
    is_point_slider_moving = false;
    normal_Slider = new QSlider(Qt::Horizontal);
    normal_Slider->setValue(20);
    point_Slider = new QSlider(Qt::Horizontal);
    point_Slider->setValue(2);
    point_Slider->setMinimum(1);
    point_Slider->setMaximum(25);
  }
  Scene_points_with_normal_item_priv(Scene_points_with_normal_item* parent)
    :m_points(new Point_set)
  {
    init_values(parent);
  }
  Scene_points_with_normal_item_priv(const Scene_points_with_normal_item& toCopy, Scene_points_with_normal_item* parent)
    : m_points(new Point_set(*toCopy.d->m_points))
  {
    init_values(parent);
  }

  Scene_points_with_normal_item_priv(const SMesh& input_mesh, Scene_points_with_normal_item* parent)
    : m_points(new Point_set)
  {
   init_values(parent);
   boost::graph_traits<SMesh>::vertex_iterator v;
    m_points->add_normal_map();
    for (v = const_cast<SMesh&>(input_mesh).vertices_begin();
         v != const_cast<SMesh&>(input_mesh).vertices_end(); v++)
    {
      boost::graph_traits<SMesh>::vertex_descriptor vd(*v);
      const Kernel::Point_3& p = input_mesh.point(vd);
      Kernel::Vector_3 n =
        CGAL::Polygon_mesh_processing::compute_vertex_normal(vd, input_mesh);
      m_points->insert(p,n);
    }
  }

  Scene_points_with_normal_item_priv(const Polyhedron& input_mesh, Scene_points_with_normal_item* parent)
     : m_points(new Point_set)
   {
    init_values(parent);
     Polyhedron::Vertex_iterator v;
     m_points->add_normal_map();
     for (v = const_cast<Polyhedron&>(input_mesh).vertices_begin();
          v != const_cast<Polyhedron&>(input_mesh).vertices_end(); v++)
     {
       const Kernel::Point_3& p = v->point();
       Kernel::Vector_3 n =
         CGAL::Polygon_mesh_processing::compute_vertex_normal(v, input_mesh);
       m_points->insert(p,n);
     }
   }

  ~Scene_points_with_normal_item_priv()
  {
    if(m_points)
    {
      delete m_points;
      m_points = NULL;
    }
    delete normal_Slider;
    delete point_Slider;
  }
  bool isPointSliderMoving() { return is_point_slider_moving; }
  void initializeBuffers(CGAL::Three::Viewer_interface *viewer) const;
  void compute_normals_and_vertices() const;
  enum VAOs {
      Edges=0,
      ThePoints,
      TheShadedPoints,
      Selected_points,
      Selected_shaded_points,
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
  mutable bool is_point_slider_moving;
  mutable std::vector<CGAL_data_type> positions_lines;
  mutable std::vector<CGAL_data_type> normals;
  mutable std::vector<CGAL_data_type> positions_normals;
  mutable std::vector<CGAL_data_type> positions_selected_normals;
  mutable std::vector<CGAL_data_type> colors_points;
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

Scene_points_with_normal_item::Scene_points_with_normal_item(const SMesh& input_mesh)
    : Scene_item(Scene_points_with_normal_item_priv::NbOfVbos,Scene_points_with_normal_item_priv::NbOfVaos)
{
  // Converts Polyhedron vertices to point set.
  // Computes vertices normal from connectivity.
  d = new Scene_points_with_normal_item_priv(input_mesh, this);
  setRenderingMode(PointsPlusNormals);
  is_selected = true;
  invalidateOpenGLBuffers();
}

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
    if(item->has_normals())
    {
        program = item->getShaderProgram(Scene_points_with_normal_item::PROGRAM_NO_SELECTION, viewer);
        program->bind();

        item->vaos[Edges]->bind();
        item->buffers[Edges_vertices].bind();
        item->buffers[Edges_vertices].allocate(positions_lines.data(),
                            static_cast<int>(positions_lines.size()*sizeof(CGAL_data_type)));
        program->enableAttributeArray("vertex");
        program->setAttributeBuffer("vertex",CGAL_GL_data_type,0,3);
        item->buffers[Edges_vertices].release();

        if (!(colors_points.empty()))
          {
            item->buffers[Points_colors].bind();
            item->buffers[Points_colors].allocate (colors_points.data(),
                                                   static_cast<int>(colors_points.size()*sizeof(CGAL_data_type)));
            program->enableAttributeArray("colors");
            program->setAttributeBuffer("colors",CGAL_GL_data_type,0,3);
            item->buffers[Points_colors].release();
          }

        item->vaos[Edges]->release();

        nb_lines = positions_lines.size();
        program->release();
    }
    //vao for the points
    {
        program = item->getShaderProgram(Scene_points_with_normal_item::PROGRAM_NO_SELECTION, viewer);
        program->bind();

        item->vaos[ThePoints]->bind();
        item->buffers[Edges_vertices].bind();
        if(!item->has_normals()) {
          item->buffers[Edges_vertices].allocate(positions_lines.data(),
                              static_cast<int>(positions_lines.size()*sizeof(CGAL_data_type)));
        }
        program->enableAttributeArray("vertex");
        if(item->has_normals())
          program->setAttributeBuffer("vertex",CGAL_GL_data_type,0,3,
                                      static_cast<int>(6*sizeof(CGAL_GL_data_type)));
        else
          program->setAttributeBuffer("vertex",CGAL_GL_data_type,0,3);
        item->buffers[Edges_vertices].release();
        if (!(colors_points.empty()))
          {
            item->buffers[Points_colors].bind();
            item->buffers[Points_colors].allocate (colors_points.data(),
                                                   static_cast<int>(colors_points.size()*sizeof(CGAL_data_type)));
            program->enableAttributeArray("colors");
            program->setAttributeBuffer("colors",CGAL_GL_data_type,0,3,6*sizeof(CGAL_data_type));
            item->buffers[Points_colors].release();
            colors_points.resize(0);
            std::vector<CGAL_data_type>(colors_points).swap(colors_points);
          }

        item->vaos[ThePoints]->release();
        program->release();
        if(item->has_normals())
        {
          program = item->getShaderProgram(Scene_points_with_normal_item::PROGRAM_WITH_LIGHT, viewer);
          item->vaos[TheShadedPoints]->bind();
          item->buffers[Edges_vertices].bind();
          program->enableAttributeArray("vertex");
          program->setAttributeBuffer("vertex",CGAL_GL_data_type,0,3,
                                      static_cast<int>(6*sizeof(CGAL_GL_data_type)));
          item->buffers[Edges_vertices].release();
          item->buffers[Points_normals].bind();
          item->buffers[Points_normals].allocate(positions_normals.data(),
                                                 static_cast<int>(positions_normals.size()*sizeof(CGAL_data_type)));
          program->enableAttributeArray("normals");
          program->setAttributeBuffer("normals",CGAL_GL_data_type,0,3);
          item->buffers[Points_normals].release();
          positions_normals.resize(0);
          std::vector<CGAL_data_type>(positions_normals).swap(positions_normals);
          item->vaos[TheShadedPoints]->release();
          program->release();
        }


        if(item->has_normals())
          nb_points = positions_lines.size()/2;
        else
          nb_points = positions_lines.size();
    }
    //vao for the selected points
    {
      program = item->getShaderProgram(Scene_points_with_normal_item::PROGRAM_NO_SELECTION, viewer);
      program->bind();

      item->vaos[Selected_points]->bind();
      item->buffers[Edges_vertices].bind();
      program->enableAttributeArray("vertex");
      if(!item->has_normals())
      {
        program->setAttributeBuffer("vertex",CGAL_GL_data_type,
                                    static_cast<int>( 3*(m_points->size()-m_points->nb_selected_points())*sizeof(CGAL_data_type) ),
                                    3,
                                    0);
      }
      else
      {
        program->setAttributeBuffer("vertex",CGAL_GL_data_type,
                                    static_cast<int>( 6*(m_points->size()-m_points->nb_selected_points())*sizeof(CGAL_data_type) ),
                                    3,
                                    static_cast<int>(6*sizeof(CGAL_GL_data_type)));
      }
      item->buffers[Edges_vertices].release();
      item->vaos[Selected_points]->release();
      program->release();
      if(item->has_normals())
      {
        program = item->getShaderProgram(Scene_points_with_normal_item::PROGRAM_WITH_LIGHT, viewer);
        item->vaos[Selected_shaded_points]->bind();
        item->buffers[Edges_vertices].bind();
        program->enableAttributeArray("vertex");
        program->setAttributeBuffer("vertex",CGAL_GL_data_type,
                                    static_cast<int>( 6*(m_points->size()-m_points->nb_selected_points())*sizeof(CGAL_data_type) ),
                                    3,
                                    static_cast<int>(6*sizeof(CGAL_GL_data_type)));

        item->buffers[Edges_vertices].release();

        item->buffers[Selected_points_normals].bind();
        item->buffers[Selected_points_normals].allocate(positions_selected_normals.data(),
                                                        static_cast<int>(positions_selected_normals.size()*sizeof(CGAL_data_type)));
        program->enableAttributeArray("normals");
        program->setAttributeBuffer("normals",CGAL_GL_data_type,0,3);
        item->buffers[Selected_points_normals].release();
        positions_selected_normals.resize(0);
        std::vector<CGAL_data_type>(positions_selected_normals).swap(positions_selected_normals);

        item->vaos[Selected_shaded_points]->release();
        program->release();
      }
      nb_selected_points = 3*m_points->nb_selected_points();
    }
    positions_lines.resize(0);
    positions_lines.shrink_to_fit();
    item->are_buffers_filled = true;
}

void Scene_points_with_normal_item_priv::compute_normals_and_vertices() const
{
    const qglviewer::Vec offset = static_cast<CGAL::Three::Viewer_interface*>(QGLViewer::QGLViewerPool().first())->offset();
    QApplication::setOverrideCursor(Qt::WaitCursor);
    positions_lines.resize(0);
    normals.resize(0);
    positions_normals.resize(0);
    positions_selected_normals.resize(0);
    normals.resize(0);
    colors_points.resize(0);

    //Shuffle container to allow quick display random points
    std::random_shuffle (m_points->begin(), m_points->first_selected());
    if (m_points->nb_selected_points() != 0)
      std::random_shuffle (m_points->first_selected(), m_points->end());
    //if item has normals, points will be one point out of two in the lines data.
    //else points will be lines and lines discarded.
    double average_spacing = 0;
    double normal_length =0;
    double length_factor =0;
    if (item->has_normals())
    {

#ifdef LINK_WITH_TBB
      typedef CGAL::Parallel_tag Concurrency_tag;
#else
      typedef CGAL::Sequential_tag Concurrency_tag;
#endif
    // Store normals
    Kernel::Sphere_3 region_of_interest = m_points->region_of_interest();
      positions_lines.reserve(m_points->size() * 6);
      positions_normals.reserve((m_points->size() - m_points->nb_selected_points()) * 3);
      positions_selected_normals.reserve(m_points->nb_selected_points() * 3);
      average_spacing = CGAL::compute_average_spacing<Concurrency_tag>(
            m_points->begin(), m_points->end(), m_points->point_map(),
            6);
      normal_length = (std::min)(average_spacing, std::sqrt(region_of_interest.squared_radius() / 1000.));
      length_factor = 5.0/100*normal_Slider->value();
    }
    else
    {
      positions_lines.reserve(m_points->size() * 3);
    }


    for (Point_set_3<Kernel>::const_iterator it = m_points->begin(); it != m_points->first_selected(); ++it)
    {
      const Kernel::Point_3& p = m_points->point(*it);
      positions_lines.push_back(p.x()+offset.x);
      positions_lines.push_back(p.y()+offset.y);
      positions_lines.push_back(p.z()+offset.z);
      if(item->has_normals())
      {
        const Kernel::Vector_3& n = m_points->normal(*it);
        Point_set_3<Kernel>::Point q = p + normal_length * length_factor* n;
        positions_lines.push_back(q.x()+offset.x);
        positions_lines.push_back(q.y()+offset.y);
        positions_lines.push_back(q.z()+offset.z);


        positions_normals.push_back(n.x());
        positions_normals.push_back(n.y());
        positions_normals.push_back(n.z());
      }
    }
    for (Point_set_3<Kernel>::const_iterator it = m_points->first_selected(); it != m_points->end(); ++it)
    {
      const Kernel::Point_3& p = m_points->point(*it);
      positions_lines.push_back(p.x()+offset.x);
      positions_lines.push_back(p.y()+offset.y);
      positions_lines.push_back(p.z()+offset.z);
      if(item->has_normals())
      {
        const Kernel::Vector_3& n = m_points->normal(*it);
        Point_set_3<Kernel>::Point q = p + normal_length * length_factor* n;
        positions_lines.push_back(q.x()+offset.x);
        positions_lines.push_back(q.y()+offset.y);
        positions_lines.push_back(q.z()+offset.z);


        positions_selected_normals.push_back(n.x());
        positions_selected_normals.push_back(n.y());
        positions_selected_normals.push_back(n.z());
      }
    }
    //The colors
    if (m_points->has_colors())
    {
        colors_points.reserve((m_points->size() - m_points->nb_selected_points()) * 6);

        for (Point_set_3<Kernel>::const_iterator it = m_points->begin(); it != m_points->end(); ++it)
	  {
            colors_points.push_back (m_points->red(*it));
            colors_points.push_back (m_points->green(*it));
            colors_points.push_back (m_points->blue(*it));
            colors_points.push_back (m_points->red(*it));
            colors_points.push_back (m_points->green(*it));
            colors_points.push_back (m_points->blue(*it));
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
  std::vector<Point_set::Index> unselected, selected;
  for (Point_set::iterator ptit = d->m_points->begin(); ptit!= d->m_points->end(); ++ ptit)
    if ( !unique_points.insert(d->m_points->point(*ptit)).second)
      selected.push_back (*ptit);
    else
      unselected.push_back (*ptit);
  
  for (std::size_t i = 0; i < unselected.size(); ++ i)
    *(d->m_points->begin() + i) = unselected[i];
  for (std::size_t i = 0; i < selected.size(); ++ i)
    *(d->m_points->begin() + (unselected.size() + i)) = selected[i];

  if (selected.empty ())
  {
    d->m_points->unselect_all();
  }
  else
  {
    d->m_points->set_first_selected
      (d->m_points->begin() + unselected.size());
  } 

  invalidateOpenGLBuffers();
  Q_EMIT itemChanged();
}

#if !defined(CGAL_CFG_NO_CPP0X_RVALUE_REFERENCE) && !defined(CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES)
#ifdef CGAL_LINKED_WITH_LASLIB
// Loads point set from .LAS file
bool Scene_points_with_normal_item::read_las_point_set(std::istream& stream)
{
  Q_ASSERT(d->m_points != NULL);

  d->m_points->clear();

  bool ok = stream &&
    CGAL::read_las_point_set (stream, *(d->m_points)) &&
            !isEmpty();

  std::cerr << d->m_points->info();

  if (d->m_points->has_normal_map())
    setRenderingMode(PointsPlusNormals);
  if (d->m_points->check_colors())
    std::cerr << "-> Point set has colors" << std::endl;
  
  invalidateOpenGLBuffers();
  return ok;
}

// Write point set to .LAS file
bool Scene_points_with_normal_item::write_las_point_set(std::ostream& stream) const
{
  Q_ASSERT(d->m_points != NULL);

  return stream &&
    CGAL::write_las_point_set (stream, *(d->m_points));
}

#endif // LAS

// Loads point set from .PLY file
bool Scene_points_with_normal_item::read_ply_point_set(std::istream& stream)
{
  Q_ASSERT(d->m_points != NULL);

  d->m_points->clear();

  bool ok = stream &&
    CGAL::read_ply_point_set (stream, *(d->m_points)) &&
            !isEmpty();

  std::cerr << d->m_points->info();

  if (d->m_points->has_normal_map())
    setRenderingMode(PointsPlusNormals);
  if (d->m_points->check_colors())
    std::cerr << "-> Point set has colors" << std::endl;
  
  invalidateOpenGLBuffers();
  return ok;
}

// Write point set to .PLY file
bool Scene_points_with_normal_item::write_ply_point_set(std::ostream& stream, bool binary) const
{
  Q_ASSERT(d->m_points != NULL);

  if (!stream)
    return false;

  if (binary)
    CGAL::set_binary_mode (stream);
  stream << *(d->m_points);

  return true;
}

#endif // CXX11

// Loads point set from .OFF file
bool Scene_points_with_normal_item::read_off_point_set(std::istream& stream)
{
  Q_ASSERT(d->m_points != NULL);

  d->m_points->clear();
  bool ok = stream &&
    CGAL::read_off_point_set(stream, *(d->m_points)) &&
            !isEmpty();
  invalidateOpenGLBuffers();
  return ok;
}

// Write point set to .OFF file
bool Scene_points_with_normal_item::write_off_point_set(std::ostream& stream) const
{
  Q_ASSERT(d->m_points != NULL);

  return stream &&
    CGAL::write_off_point_set (stream, *(d->m_points));
}

// Loads point set from .XYZ file
bool Scene_points_with_normal_item::read_xyz_point_set(std::istream& stream)
{
  Q_ASSERT(d->m_points != NULL);

  d->m_points->clear();

  bool ok = stream &&
    CGAL::read_xyz_point_set (stream, *(d->m_points)) &&
    !isEmpty();

  invalidateOpenGLBuffers();
  return ok;
}

// Write point set to .XYZ file
bool Scene_points_with_normal_item::write_xyz_point_set(std::ostream& stream) const
{
  Q_ASSERT(d->m_points != NULL);

  return stream &&
    CGAL::write_xyz_point_set (stream, *(d->m_points));
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
  switch ( m )
  {
  case Points:
    return true;
  case ShadedPoints:
  case PointsPlusNormals:
  case Splatting:
    return has_normals();

  default:
    return false;
  }
}

void Scene_points_with_normal_item::drawSplats(CGAL::Three::Viewer_interface* viewer) const
{
  const qglviewer::Vec v_offset = static_cast<CGAL::Three::Viewer_interface*>(QGLViewer::QGLViewerPool().first())->offset();
 Kernel::Vector_3 offset(v_offset.x, v_offset.y, v_offset.z);

   // TODO add support for selection
   viewer->glBegin(GL_POINTS);
   if (d->m_points->has_colors())
     for ( Point_set_3<Kernel>::const_iterator it = d->m_points->begin(); it != d->m_points->end(); it++)
       {
         Point_set::Point p = d->m_points->point (*it) + offset;
         const Point_set::Vector& n = d->m_points->normal (*it);
         viewer->glColor4d(d->m_points->red(*it),
                           d->m_points->green(*it),
                           d->m_points->blue(*it),
                           1.0);
         viewer->glNormal3dv(&n.x());
         viewer->glMultiTexCoord1d(GL_TEXTURE2, d->m_points->radius(*it));
         viewer->glVertex3dv(&p.x());

       }
   else
     for ( Point_set_3<Kernel>::const_iterator it = d->m_points->begin(); it != d->m_points->end(); it++)
       {
         const Point_set::Point p = d->m_points->point (*it) + offset;
         const Point_set::Vector& n = d->m_points->normal (*it);
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
        (d->nb_lines/6 > limit_fast_drawing)) // arbitrary large value
      ratio_displayed = 6 * limit_fast_drawing / (double)(d->nb_lines);

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
    if ((viewer->inFastDrawing () || d->isPointSliderMoving())
        &&((d->nb_points + d->nb_selected_points)/3 > limit_fast_drawing)) // arbitrary large value
      ratio_displayed = 3 * limit_fast_drawing / (double)(d->nb_points + d->nb_selected_points);

    // POINTS
    if(has_normals() && renderingMode() == ShadedPoints)
    {
      vaos[Scene_points_with_normal_item_priv::TheShadedPoints]->bind();
      d->program=getShaderProgram(PROGRAM_WITH_LIGHT);
      attribBuffers(viewer,PROGRAM_WITH_LIGHT);
    }
    else
    {
      vaos[Scene_points_with_normal_item_priv::ThePoints]->bind();
      d->program=getShaderProgram(PROGRAM_NO_SELECTION);
      attribBuffers(viewer,PROGRAM_NO_SELECTION);
    }
    d->program->bind();
    if (!(d->m_points->has_colors()) || renderingMode() == ShadedPoints)
      d->program->setAttributeValue("colors", this->color());
    viewer->glDrawArrays(GL_POINTS, 0,
                         static_cast<GLsizei>(((std::size_t)(ratio_displayed * d->nb_points)/3)));

    if(has_normals() && renderingMode() == ShadedPoints)
      vaos[Scene_points_with_normal_item_priv::TheShadedPoints]->release();
    else
      vaos[Scene_points_with_normal_item_priv::ThePoints]->release();
    d->program->release();


    // SELECTED POINTS
    vaos[Scene_points_with_normal_item_priv::Selected_points]->bind();
    if(has_normals() && renderingMode() == ShadedPoints)
    {
      vaos[Scene_points_with_normal_item_priv::Selected_shaded_points]->bind();
      d->program=getShaderProgram(PROGRAM_WITH_LIGHT);
      attribBuffers(viewer,PROGRAM_WITH_LIGHT);
    }
    else
    {
      vaos[Scene_points_with_normal_item_priv::Selected_points]->bind();
      d->program=getShaderProgram(PROGRAM_NO_SELECTION);
      attribBuffers(viewer,PROGRAM_NO_SELECTION);
    }
    d->program->bind();
    d->program->setAttributeValue("colors", QColor(255,0,0));
    viewer->glDrawArrays(GL_POINTS, 0,
                         static_cast<GLsizei>(((std::size_t)(ratio_displayed * d->nb_selected_points)/3)));

    if(has_normals() && renderingMode() == ShadedPoints)
      vaos[Scene_points_with_normal_item_priv::Selected_shaded_points]->bind();
    else
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
  typedef CGAL::Search_traits_adapter <Point_set::Index, Point_set::Point_map, SearchTraits_3> Search_traits;
  typedef CGAL::Orthogonal_k_neighbor_search<Search_traits> Neighbor_search;
  typedef Neighbor_search::Tree Tree;
  typedef Neighbor_search::Distance Distance;

  // build kdtree
  Tree tree(d->m_points->begin(),
            d->m_points->end(),
            Tree::Splitter(),
            Search_traits (d->m_points->point_map())
            );
  Distance tr_dist(d->m_points->point_map());

  if (!(d->m_points->has_property_map<double> ("radius")))
    d->m_points->add_radius();

  // Compute the radius of each point = (distance max to k nearest neighbors)/2.
  {
    int i=0;
    for (Point_set::iterator it=d->m_points->begin(); it!=d->m_points->end(); ++it, ++i)
    {
      Neighbor_search search(tree, d->m_points->point(*it), k+1, 0, true, tr_dist);
      double maxdist2 = (--search.end())->second; // squared distance to furthest neighbor
      d->m_points->radius(*it) = sqrt(maxdist2)/2.;
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
        if((d->nb_points + d->nb_selected_points)/3 <= limit_fast_drawing)
        {
          connect(d->normal_Slider, &QSlider::valueChanged, this, &Scene_points_with_normal_item::invalidateOpenGLBuffers);
          connect(d->normal_Slider, &QSlider::valueChanged, this, &Scene_points_with_normal_item::itemChanged);
        }
        else
        {
          connect(d->normal_Slider, &QSlider::sliderReleased, this, &Scene_points_with_normal_item::invalidateOpenGLBuffers);
          connect(d->normal_Slider, &QSlider::sliderReleased, this, &Scene_points_with_normal_item::itemChanged);
        }
        sliderAction->setDefaultWidget(d->normal_Slider);

        container->addAction(sliderAction);
        menu->addMenu(container);
      }
        QMenu *container = new QMenu(tr("Points Size"));
        QWidgetAction *sliderAction = new QWidgetAction(0);
        connect(d->point_Slider, &QSlider::sliderPressed, this, &Scene_points_with_normal_item::pointSliderPressed);
        connect(d->point_Slider, &QSlider::sliderReleased, this, &Scene_points_with_normal_item::pointSliderReleased);
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

bool Scene_points_with_normal_item::has_normals() const { return d->m_points->has_normal_map(); }

void Scene_points_with_normal_item::invalidateOpenGLBuffers()
{
    are_buffers_filled = false;
    compute_bbox();
}

void Scene_points_with_normal_item::pointSliderPressed()
{
  d->is_point_slider_moving = true;
}

void Scene_points_with_normal_item::pointSliderReleased()
{
  d->is_point_slider_moving = false;
}

void Scene_points_with_normal_item::copyProperties(Scene_item *item)
{
  Scene_points_with_normal_item* point_item = qobject_cast<Scene_points_with_normal_item*>(item);
  if(!point_item)
    return;
  d->point_Slider->setValue(point_item->getPointSliderValue());
  if(has_normals())
    d->normal_Slider->setValue(point_item->getNormalSliderValue());
}

int Scene_points_with_normal_item::getNormalSliderValue()
{
  return d->normal_Slider->value();
}

int Scene_points_with_normal_item::getPointSliderValue()
{
  return d->point_Slider->value();
}

void Scene_points_with_normal_item::itemAboutToBeDestroyed(Scene_item *item)
{
  Scene_item::itemAboutToBeDestroyed(item);
  if(d && d->m_points && item == this)
  {
    delete d->m_points;
    d->m_points = NULL;
  }
}

void Scene_points_with_normal_item::
zoomToPosition(const QPoint &, CGAL::Three::Viewer_interface *viewer) const
{
  if (point_set()->nb_selected_points() == 0)
    return;
  const qglviewer::Vec offset = static_cast<CGAL::Three::Viewer_interface*>(QGLViewer::QGLViewerPool().first())->offset();
  // fit plane to triangles
  Point_set points;
  Bbox selected_points_bbox;
  for(Point_set::const_iterator it = point_set()->first_selected();
      it != point_set()->end();
      ++it)
  {
    points.insert(point_set()->point(*it));
    selected_points_bbox += point_set()->point(*it).bbox();
  }
  Kernel::Plane_3 plane;
  Kernel::Point_3 center_of_mass;
  CGAL::linear_least_squares_fitting_3
      (points.points().begin(),points.points().end(),plane, center_of_mass,
       CGAL::Dimension_tag<0>());

  Kernel::Vector_3 plane_normal= plane.orthogonal_vector();
  plane_normal = plane_normal/(CGAL::sqrt(plane_normal.squared_length()));
  Kernel::Point_3 centroid(center_of_mass.x() + offset.x,
                           center_of_mass.y() + offset.y,
                           center_of_mass.z() + offset.z);

  qglviewer::Quaternion new_orientation(qglviewer::Vec(0,0,-1),
                                        qglviewer::Vec(-plane_normal.x(), -plane_normal.y(), -plane_normal.z()));
  double max_side = (std::max)((std::max)(selected_points_bbox.xmax() - selected_points_bbox.xmin(),
                                          selected_points_bbox.ymax() - selected_points_bbox.ymin()),
                               selected_points_bbox.zmax() - selected_points_bbox.zmin());
  //put the camera in way we are sure the longest side is entirely visible on the screen
  //See openGL's frustum definition
  double factor = max_side/(tan(viewer->camera()->aspectRatio()/
                                  (viewer->camera()->fieldOfView()/2)));

  Kernel::Point_3 new_pos = centroid + factor*plane_normal ;
  viewer->camera()->setSceneCenter(qglviewer::Vec(centroid.x(),
                                                  centroid.y(),
                                                  centroid.z()));
  viewer->moveCameraToCoordinates(QString("%1 %2 %3 %4 %5 %6 %7").arg(new_pos.x())
                                                                 .arg(new_pos.y())
                                                                 .arg(new_pos.z())
                                                                 .arg(new_orientation[0])
                                                                 .arg(new_orientation[1])
                                                                 .arg(new_orientation[2])
                                                                 .arg(new_orientation[3]));

}
