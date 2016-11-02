#include "Scene_point_set_classification_item.h"
#include "Color_ramp.h"

#include <CGAL/Timer.h>
#include <CGAL/Memory_sizer.h>

#include <CGAL/Three/Viewer_interface.h>

#include <QObject>
#include <QMenu>
#include <QGLViewer/manipulatedCameraFrame.h>
 
#include <set>
#include <stack>
#include <algorithm>
#include <boost/array.hpp>

Scene_point_set_classification_item::Scene_point_set_classification_item(PSC* psc)
: Scene_item(NbOfVbos,NbOfVaos),
  m_points (NULL),
  m_psc (psc),
  m_helper (NULL)
{
  setRenderingMode(PointsPlusNormals);
  m_nb_scales = 5;
  m_nb_trials = 300;
  m_smoothing = 0.5;
  m_index_color = 1;
  is_selected = true;
  nb_points = 0;
  nb_lines = 0;
}

Scene_point_set_classification_item::Scene_point_set_classification_item(Scene_points_with_normal_item* points)
  : Scene_item(NbOfVbos,NbOfVaos),
    m_points (points),
    m_psc (NULL),
    m_helper (NULL)
{
  setRenderingMode(PointsPlusNormals);
  m_nb_scales = 5;
  m_index_color = 1;
  m_nb_trials = 300;
  m_smoothing = 0.5;

  reset_indices();
  
  m_psc = new PSC(m_points->point_set()->begin(), m_points->point_set()->end(), m_points->point_set()->point_map());

  Type_handle ground = m_psc->add_classification_type("ground");
  Type_handle vegetation = m_psc->add_classification_type("vegetation");
  Type_handle roof = m_psc->add_classification_type("roof");
  Type_handle facade = m_psc->add_classification_type("facade");
  m_types.push_back (std::make_pair(ground, QColor(245, 180, 0)));
  m_types.push_back (std::make_pair(vegetation, QColor(0, 255, 27)));
  m_types.push_back (std::make_pair(roof, QColor(255, 0, 170)));
  m_types.push_back (std::make_pair(facade, QColor(100, 0, 255)));
  
  is_selected = true;
  nb_points = 0;
  nb_lines = 0;
}


// Copy constructor
Scene_point_set_classification_item::Scene_point_set_classification_item(const Scene_point_set_classification_item&)
  :Scene_item(NbOfVbos,NbOfVaos), // do not call superclass' copy constructor
   m_points (NULL),
   m_psc (NULL),
   m_helper (NULL)
{
  setRenderingMode(PointsPlusNormals);
  m_nb_scales = 5;
  m_index_color = 1;
  m_nb_trials = 300;
  m_smoothing = 0.5;
  
  nb_points = 0;
  invalidateOpenGLBuffers();
}

Scene_point_set_classification_item::~Scene_point_set_classification_item()
{
  if (m_psc != NULL)
    delete m_psc;
  if (m_helper != NULL)
    delete m_helper;
}

void Scene_point_set_classification_item::initializeBuffers(CGAL::Three::Viewer_interface *viewer) const
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
    if (!(colors_points.empty()))
      {
        buffers[Points_colors].bind();
        buffers[Points_colors].allocate (colors_points.data(),
                                         static_cast<int>(colors_points.size()*sizeof(double)));
        program->enableAttributeArray("colors");
        program->setAttributeBuffer("colors",GL_DOUBLE,0,3);
        buffers[Points_colors].release();
        std::vector<double>(colors_points).swap(colors_points);
      }

    vaos[ThePoints]->release();
    nb_points = positions_points.size();
    positions_points.resize(0);
    std::vector<double>(positions_points).swap(positions_points);
    program->release();
  }
  are_buffers_filled = true;
}

void Scene_point_set_classification_item::compute_normals_and_vertices() const
{
  int index_color = real_index_color();
    
  positions_lines.resize(0);
  positions_points.resize(0);
  colors_points.resize(0);
  if (index_color == 9) // Show clusters centroids
    {
      // positions_points.reserve(m_psc->clusters().size() * 3);
      // colors_points.reserve(m_psc->clusters().size() * 3);

      // for (std::size_t i = 0; i < m_psc->clusters().size(); ++ i)
      //   {
      //     const Kernel::Point_3& p = m_psc->clusters()[i].centroid;
      //     positions_points.push_back(p.x());
      //     positions_points.push_back(p.y());
      //     positions_points.push_back(p.z());

      //     for (std::set<std::size_t>::iterator it = m_psc->clusters()[i].neighbors.begin ();
      //          it != m_psc->clusters()[i].neighbors.end (); ++ it)
      //       {
      //         const Kernel::Point_3& q = m_psc->clusters()[*it].centroid;
      //         positions_lines.push_back(p.x());
      //         positions_lines.push_back(p.y());
      //         positions_lines.push_back(p.z());

      //         positions_lines.push_back(q.x());
      //         positions_lines.push_back(q.y());
      //         positions_lines.push_back(q.z());
      //       }
      //   }
    }
  else // Show points
    {
      positions_points.reserve(m_points->point_set()->size() * 3);
      colors_points.reserve(m_points->point_set()->size() * 3);

      for (Point_set::const_iterator it = m_points->point_set()->begin();
           it != m_points->point_set()->end(); ++ it)
        {
          const Kernel::Point_3& p = m_points->point_set()->point(*it);
          positions_points.push_back(p.x());
          positions_points.push_back(p.y());
          positions_points.push_back(p.z());
        }
    }

  // Colors
  static Color_ramp ramp;
  ramp.build_red();

  if (index_color == -1) // item color
    {

      for (Point_set::const_iterator it = m_points->point_set()->begin();
           it != m_points->point_set()->first_selected(); ++ it)
        {
          colors_points.push_back (color().redF());
          colors_points.push_back (color().greenF());
          colors_points.push_back (color().blueF());
        }
    }
  else if (index_color == 0) // real colors
    {

      for (Point_set::const_iterator it = m_points->point_set()->begin();
           it != m_points->point_set()->first_selected(); ++ it)
        {
          colors_points.push_back ((double)(m_points->point_set()->red(*it)) / 255.);
          colors_points.push_back ((double)(m_points->point_set()->green(*it)) / 255.);
          colors_points.push_back ((double)(m_points->point_set()->blue(*it)) / 255.);
        }
    }
  else if (index_color == 1) // classif
    {
      std::map<Type_handle, QColor> map_colors;
      for (std::size_t i = 0; i < m_types.size(); ++ i)
        map_colors.insert (m_types[i]);
          
      for (Point_set::const_iterator it = m_points->point_set()->begin();
           it != m_points->point_set()->first_selected(); ++ it)
        {
          QColor color (0, 0, 0);
          Type_handle c = m_psc->classification_type_of(*it);
          
          if (c != Type_handle())
            color = map_colors[c];

          colors_points.push_back ((double)(color.red()) / 255.);
          colors_points.push_back ((double)(color.green()) / 255.);
          colors_points.push_back ((double)(color.blue()) / 255.);
        }
    }
  else if (index_color == 2) // training
    {
      std::map<Type_handle, QColor> map_colors;
      for (std::size_t i = 0; i < m_types.size(); ++ i)
        map_colors.insert (m_types[i]);
          
      for (Point_set::const_iterator it = m_points->point_set()->begin();
           it != m_points->point_set()->first_selected(); ++ it)
        {
          QColor color (0, 0, 0);
          Type_handle c = m_psc->training_type_of(*it);
          Type_handle c2 = m_psc->classification_type_of(*it);
          
          if (c != Type_handle())
            color = map_colors[c];
          double div = 255.;
          if (c != c2)
            div = 400.;
          
          colors_points.push_back ((double)(color.red()) / div);
          colors_points.push_back ((double)(color.green()) / div);
          colors_points.push_back ((double)(color.blue()) / div);
        }
    }
  else
    {
      Attribute_handle att = m_psc->get_attribute(index_color - 3);
      double weight = att->weight;
      att->weight = att->max;
      for (Point_set::const_iterator it = m_points->point_set()->begin();
           it != m_points->point_set()->first_selected(); ++ it)
        {
          colors_points.push_back (ramp.r(att->normalized(*it)));
          colors_points.push_back (ramp.g(att->normalized(*it)));
          colors_points.push_back (ramp.b(att->normalized(*it)));
        }
      att->weight = weight;
    }

  for (Point_set::const_iterator it = m_points->point_set()->first_selected();
       it != m_points->point_set()->end(); ++ it)
    {
      colors_points.push_back (1.);
      colors_points.push_back (0.);
      colors_points.push_back (0.);
    }
}


// Duplicates scene item
Scene_point_set_classification_item*
Scene_point_set_classification_item::clone() const
{
  return new Scene_point_set_classification_item(*this);
}

// Write point set to .PLY file
bool Scene_point_set_classification_item::write_ply_point_set(std::ostream& stream)
{
  if (m_helper == NULL)
    return false;
    
  stream.precision (std::numeric_limits<double>::digits10 + 2);

  reset_indices();

  std::vector<Color> colors;
  for (std::size_t i = 0; i < m_types.size(); ++ i)
    {
      Color c = {{ (unsigned char)(m_types[i].second.red()),
                   (unsigned char)(m_types[i].second.green()),
                   (unsigned char)(m_types[i].second.blue()) }};
      colors.push_back (c);
    }
  
  m_helper->write_ply (stream,
                       m_points->point_set()->begin(),
                       m_points->point_set()->end(),
                       m_points->point_set()->point_map(),
                       *m_psc,
                       &colors);
  return true;
}

QString
Scene_point_set_classification_item::toolTip() const
{
  return QObject::tr("<p><b>%1</b><br />"
                     "<i>Point set classification</i></p>"
                     "<p>Number of points: %2</p>")
    .arg(name())
    .arg(m_points->point_set()->size());
}

bool Scene_point_set_classification_item::supportsRenderingMode(RenderingMode m) const 
{
  return m==PointsPlusNormals;

}

void Scene_point_set_classification_item::drawSplats(CGAL::Three::Viewer_interface*) const
{
}

void Scene_point_set_classification_item::drawEdges(CGAL::Three::Viewer_interface* viewer) const
{
  if(!are_buffers_filled)
    initializeBuffers(viewer);
  vaos[Edges]->bind();
  program=getShaderProgram(PROGRAM_NO_SELECTION);
  attribBuffers(viewer,PROGRAM_NO_SELECTION);
  program->bind();
  program->setAttributeValue("colors", this->color());
  viewer->glDrawArrays(GL_LINES, 0,
                       static_cast<GLsizei>((nb_lines)/3));
  vaos[Edges]->release();
  program->release();
}

void Scene_point_set_classification_item::drawPoints(CGAL::Three::Viewer_interface* viewer) const
{
  if(!are_buffers_filled)
    initializeBuffers(viewer);
  GLfloat point_size;
  viewer->glGetFloatv(GL_POINT_SIZE, &point_size);
  if (m_index_color == 9)
    viewer->glPointSize(15.f);
  else
    viewer->glPointSize(3.f);

  vaos[ThePoints]->bind();
  program=getShaderProgram(PROGRAM_NO_SELECTION);
  attribBuffers(viewer,PROGRAM_NO_SELECTION);
  program->bind();
  if (colors_points.empty())
    program->setAttributeValue("colors", this->color());
  viewer->glDrawArrays(GL_POINTS, 0,
                       static_cast<GLsizei>(((std::size_t)(nb_points)/3)));
  vaos[ThePoints]->release();
  program->release();
  viewer->glPointSize(point_size);
}

bool
Scene_point_set_classification_item::isEmpty() const
{
  return false;
}

void
Scene_point_set_classification_item::compute_bbox() const
{
  if (m_points->point_set()->empty())
    return;
  
  double xmin = std::numeric_limits<double>::max();
  double ymin = std::numeric_limits<double>::max();
  double zmin = std::numeric_limits<double>::max();
  double xmax = -std::numeric_limits<double>::max();
  double ymax = -std::numeric_limits<double>::max();
  double zmax = -std::numeric_limits<double>::max();

  for (Point_set::const_iterator it = m_points->point_set()->begin();
       it != m_points->point_set()->end(); ++ it)
    {
      const Kernel::Point_3& p = m_points->point_set()->point(*it);
      xmin = (std::min)(xmin, p.x());
      ymin = (std::min)(ymin, p.y());
      zmin = (std::min)(zmin, p.z());
      xmax = (std::max)(xmax, p.x());
      ymax = (std::max)(ymax, p.y());
      zmax = (std::max)(zmax, p.z());
    }
  
  _bbox = Bbox(xmin,ymin,zmin,
               xmax,ymax,zmax);
}

QMenu* Scene_point_set_classification_item::contextMenu()
{
  QMenu* menu = Scene_item::contextMenu();

  return menu;
}

void Scene_point_set_classification_item::setRenderingMode(RenderingMode m)
{
  Scene_item::setRenderingMode(m);
}

void Scene_point_set_classification_item::invalidateOpenGLBuffers()
{
  are_buffers_filled = false;
  compute_bbox();
}

void Scene_point_set_classification_item::change_color (int index)
{
  m_index_color = index;
  invalidateOpenGLBuffers();
}

int Scene_point_set_classification_item::real_index_color() const
{
  int out = m_index_color;
  
  if (out == 0 && !(m_points->point_set()->has_colors()))
    out = -1;
  return out;
}

void Scene_point_set_classification_item::reset_indices ()
{
  Point_set::Property_map<Point_set::Index> indices;

  boost::tie (indices, boost::tuples::ignore)
    = m_points->point_set()->property_map<Point_set::Index>("index");

  m_points->point_set()->unselect_all();
  Point_set::Index idx;
  ++ idx;
  for (std::size_t i = 0; i < m_points->point_set()->size(); ++ i)
    *(indices.begin() + i) = idx ++;
}

void Scene_point_set_classification_item::compute_features ()
{
  Q_ASSERT (!(m_points->point_set()->empty()));
  Q_ASSERT (m_psc != NULL);
  m_psc->clear_attributes();
  reset_indices();
  
  std::cerr << "Computing features with " << m_nb_scales << " scale(s)" << std::endl;
  compute_bbox();
  if (m_helper != NULL) delete m_helper;


  m_helper = new Helper (m_points->point_set()->begin(),
                         m_points->point_set()->end(),
                         m_points->point_set()->point_map(),
                         m_nb_scales);
  
  m_helper->generate_point_based_attributes (*m_psc,
                                             m_points->point_set()->begin(),
                                             m_points->point_set()->end(),
                                             m_points->point_set()->point_map());

  if (m_points->point_set()->has_normal_map())
    m_helper->generate_normal_based_attributes (*m_psc,
                                                m_points->point_set()->begin(),
                                                m_points->point_set()->end(),
                                                m_points->point_set()->normal_map());
  else
    m_helper->generate_normal_based_attributes (*m_psc,
                                                m_points->point_set()->begin(),
                                                m_points->point_set()->end());

  Point_set::Property_map<boost::uint8_t> echo_map;
  bool okay;
  boost::tie (echo_map, okay) = m_points->point_set()->template property_map<boost::uint8_t>("echo");
  if (okay)
    m_helper->generate_echo_based_attributes (*m_psc,
                                              m_points->point_set()->begin(),
                                              m_points->point_set()->end(),
                                              echo_map);

  if (m_points->point_set()->has_colors())
    m_helper->generate_color_based_attributes (*m_psc,
                                               m_points->point_set()->begin(),
                                               m_points->point_set()->end(),
                                               Color_map(m_points->point_set()));

}



void Scene_point_set_classification_item::train()
{
  if (m_helper == NULL)
    {
      std::cerr << "Error: features not computed" << std::endl;
      return;
    }

  m_psc->training(m_nb_trials);
  m_psc->run();
  m_helper->info();

}

bool Scene_point_set_classification_item::run (int method)
{
  if (m_helper == NULL)
    {
      std::cerr << "Error: features not computed" << std::endl;
      return false;
    }
  reset_indices();
  
  if (method == 0)
    m_psc->run();
  else if (method == 1)
    m_psc->run_with_local_smoothing (m_helper->neighborhood(),
                                     3. * m_helper->radius_neighbors());
  else if (method == 2)
    m_psc->run_with_graphcut (m_helper->neighborhood(), m_smoothing);
  invalidateOpenGLBuffers();
  Q_EMIT itemChanged();
  
  return true;
}
