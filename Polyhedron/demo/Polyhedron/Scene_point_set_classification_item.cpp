#include "Scene_point_set_classification_item.h"
#include "Color_ramp.h"

#include <CGAL/IO/read_ply_points.h>
#include <CGAL/IO/write_ply_points.h>
#include <CGAL/Timer.h>
#include <CGAL/Memory_sizer.h>

#include <CGAL/Three/Viewer_interface.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/jet_estimate_normals.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Shape_detection_3.h>
#include <CGAL/regularize_planes.h>
#include <CGAL/compute_average_spacing.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Alpha_shape_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>

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
  m_grid_resolution = 0.1;
  m_radius_neighbors = 0.5;
  m_radius_dtm = 2.5;
  m_nb_scales = 1;
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
  m_grid_resolution = 0.1;
  m_radius_neighbors = 0.5;
  m_radius_dtm = 2.5;
  m_nb_scales = 1;
  m_index_color = 1;
  m_nb_trials = 300;
  m_smoothing = 0.5;

  m_points->point_set()->unselect_all();
  m_points->point_set()->collect_garbage();
  
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
  m_grid_resolution = 0.1;
  m_radius_neighbors = 0.5;
  m_radius_dtm = 2.5;
  m_nb_scales = 1;
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
  else if (index_color == 2) // confidence
    {
      std::map<Type_handle, QColor> map_colors;
      for (std::size_t i = 0; i < m_types.size(); ++ i)
        map_colors.insert (m_types[i]);
          
      for (Point_set::const_iterator it = m_points->point_set()->begin();
           it != m_points->point_set()->first_selected(); ++ it)
        {
          QColor color (0, 0, 0);
          Type_handle c = m_psc->training_type_of(*it);
          
          if (c != Type_handle())
            color = map_colors[c];

          double confidence = m_psc->confidence_of(*it);
          colors_points.push_back (1. - confidence * (1. - (double)(color.red()) / 255.));
          colors_points.push_back (1. - confidence * (1. - (double)(color.green()) / 255.));
          colors_points.push_back (1. - confidence * (1. - (double)(color.blue()) / 255.));
        }
    }
  else if (index_color == 3) // RANSAC
    {
      int seed = time(NULL);
        
      for (Point_set::const_iterator it = m_points->point_set()->begin();
           it != m_points->point_set()->first_selected(); ++ it)
        {
          if (m_psc->group_of(*it) == (std::size_t)(-1))
            {
              colors_points.push_back (0.);
              colors_points.push_back (0.);
              colors_points.push_back (0.);
            }
          else
            {
              srand (m_psc->group_of(*it) + seed);
              colors_points.push_back (0.25 + 0.6 * (rand() / (double)RAND_MAX));
              colors_points.push_back (0.25 + 0.6 * (rand() / (double)RAND_MAX));
              colors_points.push_back (0.25 + 0.6 * (rand() / (double)RAND_MAX));
            }
        }
    }
  else
    {
      Attribute_handle att = m_psc->get_attribute(index_color - 4);
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
  if (m_psc->classification_type_of(0) == Type_handle())
    {
      std::cerr << "Error: classification was not performed." << std::endl;
      return false;
    }

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
  typename Point_set::Property_map<Point_set::Index> indices;

  boost::tie (indices, boost::tuples::ignore)
    = m_points->point_set()->property_map<Point_set::Index>("index");

  m_points->point_set()->unselect_all();
  Point_set::Index idx;
  ++ idx;
  for (std::size_t i = 0; i < m_points->point_set()->size(); ++ i)
    *(indices.begin() + i) = idx ++;
}

void Scene_point_set_classification_item::estimate_parameters ()
{
  double average_spacing = CGAL::compute_average_spacing<CGAL::Sequential_tag> (m_points->point_set()->begin(),
                                                                                m_points->point_set()->end(),
                                                                                m_points->point_set()->point_map(), 6);
  m_grid_resolution = average_spacing;
  m_radius_neighbors = 5 * m_grid_resolution;
  m_radius_dtm = 5 * m_radius_neighbors;
}

void Scene_point_set_classification_item::compute_features ()
{
  Q_ASSERT (!(m_points->point_set()->empty()));
  Q_ASSERT (m_psc != NULL);
  m_psc->clear_attributes();
  reset_indices();
  
  std::cerr << "Computing features with parameters "
            << m_grid_resolution << " " << m_radius_neighbors << " and " << m_radius_dtm << std::endl;
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

  typename Point_set::Property_map<boost::uint8_t> echo_map;
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


void Scene_point_set_classification_item::compute_ransac ()
{
  Q_ASSERT (m_psc != NULL);

  if (!(m_points->point_set()->has_normal_map()))
    {
      std::cerr << "Computing normals..." << std::endl;
      m_points->point_set()->add_normal_map();

      CGAL::jet_estimate_normals<CGAL::Sequential_tag>(m_points->point_set()->begin(),
                                                       m_points->point_set()->end(),
                                                       m_points->point_set()->point_map(),
                                                       m_points->point_set()->normal_map(),
                                                       12);
    }
  
  typedef CGAL::Shape_detection_3::Efficient_RANSAC_traits<Kernel, Point_set,
                                                           Point_set::Point_map,
                                                           Point_set::Vector_map> Traits;
  typedef CGAL::Shape_detection_3::Efficient_RANSAC<Traits> Shape_detection;

  Shape_detection shape_detection;

  shape_detection.set_input(*m_points->point_set(), m_points->point_set()->point_map(), m_points->point_set()->normal_map());

  shape_detection.add_shape_factory<My_plane<Traits> >();
  Shape_detection::Parameters op;
  op.probability = 0.05;
  op.min_points = std::max ((std::size_t)10, m_points->point_set()->size() / 250000);
  
  op.epsilon = m_radius_neighbors;
  op.cluster_epsilon = m_radius_neighbors;
  op.normal_threshold = 0.9;

  std::cerr << "Computing RANSAC..." << std::endl;
  shape_detection.detect(op);
  
  CGAL::regularize_planes (shape_detection, true, true, true, false,
                           10., m_radius_neighbors);
  
  m_psc->reset_groups();

  std::cerr << "Storing groups..." << std::endl;

  std::size_t nb_planes = 0;
  BOOST_FOREACH(boost::shared_ptr<Shape_detection::Shape> shape, shape_detection.shapes())
    {
      m_psc->add_plane ((Kernel::Plane_3)(*(dynamic_cast<CGAL::Shape_detection_3::Plane<Traits>*>(shape.get ()))));
      BOOST_FOREACH(std::size_t i, shape->indices_of_assigned_points())
        m_psc->set_group_of(shape->indices_of_assigned_points()[i], nb_planes);
      ++ nb_planes;
    }
  std::cerr << "Found " << nb_planes << " group(s)" << std::endl;
}


void Scene_point_set_classification_item::compute_clusters ()
{
  Q_ASSERT (m_psc != NULL);
  if (m_helper == NULL)
    {
      std::cerr << "Error: no neighborhood" << std::endl;
      return;
    }

  m_psc->cluster_points (m_helper->neighborhood(), m_radius_neighbors);

  invalidateOpenGLBuffers();
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
    m_psc->run_with_graphcut (m_helper->neighborhood(), m_smoothing);
  else if (method == 2)
    m_psc->run_with_groups (m_helper->neighborhood(), m_radius_neighbors);
  invalidateOpenGLBuffers();
  Q_EMIT itemChanged();
  
  return true;
}
