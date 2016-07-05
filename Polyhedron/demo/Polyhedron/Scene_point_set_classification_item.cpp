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
  m_psc (psc),
  m_scat (NULL),
  m_elev (NULL),
  m_hori (NULL),
  m_plan (NULL),
  m_colo (NULL)
{
  setRenderingMode(PointsPlusNormals);
  m_index_color = 1;
  is_selected = true;
  nb_points = 0;
  nb_lines = 0;
}

Scene_point_set_classification_item::Scene_point_set_classification_item(const Scene_points_with_normal_item* points,
                                                                         double grid_resolution)
  : Scene_item(NbOfVbos,NbOfVaos),
    m_psc (NULL),
    m_scat (NULL),
    m_elev (NULL),
    m_hori (NULL),
    m_plan (NULL),
    m_colo (NULL)
{
  setRenderingMode(PointsPlusNormals);
  m_index_color = 1;
  m_psc = new PSC(points->point_set()->begin(), points->point_set()->end(), grid_resolution);
  
  is_selected = true;
  nb_points = 0;
  nb_lines = 0;
}


// Copy constructor
Scene_point_set_classification_item::Scene_point_set_classification_item(const Scene_point_set_classification_item&)
  :Scene_item(NbOfVbos,NbOfVaos), // do not call superclass' copy constructor
   m_psc (NULL),
   m_scat (NULL),
   m_elev (NULL),
   m_hori (NULL),
   m_plan (NULL),
   m_colo (NULL)
{
  setRenderingMode(PointsPlusNormals);
  m_index_color = 1;
  
  nb_points = 0;
  invalidateOpenGLBuffers();
}

Scene_point_set_classification_item::~Scene_point_set_classification_item()
{
  if (m_psc != NULL)
    delete m_psc;
  if (m_scat != NULL)
    delete m_scat;
  if (m_elev != NULL)
    delete m_elev;
  if (m_hori != NULL)
    delete m_hori;
  if (m_plan != NULL)
    delete m_plan;
  if (m_colo != NULL)
    delete m_colo;
}

void Scene_point_set_classification_item::initialize_buffers(CGAL::Three::Viewer_interface *viewer) const
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
      positions_points.reserve(m_psc->clusters.size() * 3);
      colors_points.reserve(m_psc->clusters.size() * 3);

      for (std::size_t i = 0; i < m_psc->clusters.size(); ++ i)
        {
          const Kernel::Point_3& p = m_psc->clusters[i].centroid;
          positions_points.push_back(p.x());
          positions_points.push_back(p.y());
          positions_points.push_back(p.z());

          for (std::set<std::size_t>::iterator it = m_psc->clusters[i].neighbors.begin ();
               it != m_psc->clusters[i].neighbors.end (); ++ it)
            {
              const Kernel::Point_3& q = m_psc->clusters[*it].centroid;
              positions_lines.push_back(p.x());
              positions_lines.push_back(p.y());
              positions_lines.push_back(p.z());

              positions_lines.push_back(q.x());
              positions_lines.push_back(q.y());
              positions_lines.push_back(q.z());
            }
        }
    }
  else // Show points
    {
      positions_points.reserve(m_psc->HPS.size() * 3);
      colors_points.reserve(m_psc->HPS.size() * 3);

      for (std::size_t i = 0; i < m_psc->HPS.size(); ++ i)
        {
          const Kernel::Point_3& p = m_psc->HPS[i].position;
          positions_points.push_back(p.x());
          positions_points.push_back(p.y());
          positions_points.push_back(p.z());
        }
    }

  // Colors
  static Color_ramp ramp;
  ramp.build_red();
    

  if (index_color == 0) // real colors
    {
      for (std::size_t i = 0; i < m_psc->HPS.size(); ++ i)
        {
          Color color = m_psc->HPS[i].color;

          colors_points.push_back ((double)(color[0]) / 255.);
          colors_points.push_back ((double)(color[1]) / 255.);
          colors_points.push_back ((double)(color[2]) / 255.);
        }
    }
  else if (index_color == 1) // classif
    {
      for (std::size_t i = 0; i < m_psc->HPS.size(); ++ i)
        {
          QColor color (0, 0, 0);
          int c = m_psc->HPS[i].AE_label;
          
          if (c != -1)
            {        
              if (m_psc->segmentation_classes[c]->id() == "vegetation")
                color = QColor(0, 255, 27);
              else if (m_psc->segmentation_classes[c]->id() == "ground")
                color = QColor(245, 180, 0);
              else if (m_psc->segmentation_classes[c]->id() == "road")
                color = QColor(114, 114, 130);
              else if (m_psc->segmentation_classes[c]->id() == "roof")
                color = QColor(255, 0, 170);
              else if (m_psc->segmentation_classes[c]->id() == "facade")
                color = QColor(100, 0, 255);
              else if (m_psc->segmentation_classes[c]->id() == "building")
                color = QColor(0, 114, 225);
            }


          colors_points.push_back ((double)(color.red()) / 255.);
          colors_points.push_back ((double)(color.green()) / 255.);
          colors_points.push_back ((double)(color.blue()) / 255.);
        }
    }
  else if (index_color == 2) // confidence
    {
      for (std::size_t i = 0; i < m_psc->HPS.size(); ++ i)
        {
          QColor color (0, 0, 0);
          int c = m_psc->HPS[i].AE_label;
          
          if (c != -1)
            {        
              if (m_psc->segmentation_classes[c]->id() == "vegetation")
                color = QColor(0, 255, 27);
              else if (m_psc->segmentation_classes[c]->id() == "ground")
                color = QColor(245, 180, 0);
              else if (m_psc->segmentation_classes[c]->id() == "road")
                color = QColor(114, 114, 130);
              else if (m_psc->segmentation_classes[c]->id() == "roof")
                color = QColor(255, 0, 170);
              else if (m_psc->segmentation_classes[c]->id() == "facade")
                color = QColor(100, 0, 255);
              else if (m_psc->segmentation_classes[c]->id() == "building")
                color = QColor(0, 114, 225);
            }

          double confidence = m_psc->HPS[i].confidence;
          colors_points.push_back (1. - confidence * (1. - (double)(color.red()) / 255.));
          colors_points.push_back (1. - confidence * (1. - (double)(color.green()) / 255.));
          colors_points.push_back (1. - confidence * (1. - (double)(color.blue()) / 255.));
        }
    }
  else if (index_color == 3) // scatter
    {
      double weight = m_scat->weight;
      m_scat->weight = m_scat->max;
      for (std::size_t i = 0; i < m_psc->HPS.size(); ++ i)
        {
          colors_points.push_back (ramp.r(m_scat->value(i)));
          colors_points.push_back (ramp.g(m_scat->value(i)));
          colors_points.push_back (ramp.b(m_scat->value(i)));
        }
      m_scat->weight = weight;
    }
  else if (index_color == 4) // planarity
    {
      double weight = m_plan->weight;
      m_plan->weight = m_plan->max;
      for (std::size_t i = 0; i < m_psc->HPS.size(); ++ i)
        {
          colors_points.push_back (ramp.r(m_plan->value(i)));
          colors_points.push_back (ramp.g(m_plan->value(i)));
          colors_points.push_back (ramp.b(m_plan->value(i)));
        }
      m_plan->weight = weight;
    }
  else if (index_color == 5) // horizontality
    {
      double weight = m_hori->weight;
      m_hori->weight = m_hori->max;
      for (std::size_t i = 0; i < m_psc->HPS.size(); ++ i)
        {
          colors_points.push_back (ramp.r(m_hori->value(i)));
          colors_points.push_back (ramp.g(m_hori->value(i)));
          colors_points.push_back (ramp.b(m_hori->value(i)));
        }
      m_hori->weight = weight;
    }
  else if (index_color == 6) // elevation
    {
      double weight = m_elev->weight;
      m_elev->weight = m_elev->max;
      for (std::size_t i = 0; i < m_psc->HPS.size(); ++ i)
        {
          colors_points.push_back (ramp.r(m_elev->value(i)));
          colors_points.push_back (ramp.g(m_elev->value(i)));
          colors_points.push_back (ramp.b(m_elev->value(i)));
        }
      m_elev->weight = weight;
    }
  else if (index_color == 7) // color
    {
      double weight = m_colo->weight;
      m_colo->weight = m_colo->max;
      for (std::size_t i = 0; i < m_psc->HPS.size(); ++ i)
        {
          colors_points.push_back (ramp.r(m_colo->value(i)));
          colors_points.push_back (ramp.g(m_colo->value(i)));
          colors_points.push_back (ramp.b(m_colo->value(i)));
        }
      m_colo->weight = weight;
    }
  else if (index_color == 8) // RANSAC
    {
      int seed = time(NULL);
        
      for (std::size_t i = 0; i < m_psc->HPS.size(); ++ i)
        {
          if (m_psc->HPS[i].group == (std::size_t)(-1))
            {
              colors_points.push_back (0.);
              colors_points.push_back (0.);
              colors_points.push_back (0.);
            }
          else
            {
              srand (m_psc->HPS[i].group + seed);
              colors_points.push_back (0.25 + 0.6 * (rand() / (double)RAND_MAX));
              colors_points.push_back (0.25 + 0.6 * (rand() / (double)RAND_MAX));
              colors_points.push_back (0.25 + 0.6 * (rand() / (double)RAND_MAX));
            }
        }
    }
  else if (index_color == 9) // Clusters
    {
      for (std::size_t i = 0; i < m_psc->clusters.size(); ++ i)
        {
          QColor color (0, 0, 0);
          int c = m_psc->HPS[m_psc->clusters[i].indices[0]].AE_label;
          
          if (c != -1)
            {        
              if (m_psc->segmentation_classes[c]->id() == "vegetation")
                color = QColor(0, 255, 27);
              else if (m_psc->segmentation_classes[c]->id() == "ground")
                color = QColor(245, 180, 0);
              else if (m_psc->segmentation_classes[c]->id() == "road")
                color = QColor(114, 114, 130);
              else if (m_psc->segmentation_classes[c]->id() == "roof")
                color = QColor(255, 0, 170);
              else if (m_psc->segmentation_classes[c]->id() == "facade")
                color = QColor(100, 0, 255);
              else if (m_psc->segmentation_classes[c]->id() == "building")
                color = QColor(0, 114, 225);
            }


          colors_points.push_back ((double)(color.red()) / 255.);
          colors_points.push_back ((double)(color.green()) / 255.);
          colors_points.push_back ((double)(color.blue()) / 255.);
        }
    }
}


// Duplicates scene item
Scene_point_set_classification_item*
Scene_point_set_classification_item::clone() const
{
  return new Scene_point_set_classification_item(*this);
}

// Loads point set from .PLY file
bool Scene_point_set_classification_item::read_ply_point_set(std::istream& stream,
                                                             double grid_resolution)
{
  CGAL::Timer timer;
  timer.start();
  std::vector<Kernel::Point_3> points;
  std::vector<Color> colors;
  std::vector<unsigned char> echo;

  My_ply_interpreter<float> interpreter_f (points, colors, echo);

  if (!CGAL::read_ply_custom_points (stream, interpreter_f, Kernel()))
    {
      std::cerr << "PLY reader with float not applicable." << std::endl;
      stream.seekg(0);
      My_ply_interpreter<double> interpreter_d (points, colors, echo);
      if (!CGAL::read_ply_custom_points (stream, interpreter_d, Kernel()))
        {
          std::cerr << "PLY reader with double not applicable." << std::endl;
          stream.seekg(0);
          if (!CGAL::read_ply_points (stream, std::back_inserter (points)))
            {
              std::cerr << "Error: cannot read file " << std::endl;
              return false;
            }
          m_psc = new PSC(points.begin(), points.end(), grid_resolution);
          invalidateOpenGLBuffers();
          return true;
        }
    }

  m_psc = new PSC(points.begin(), points.end(), grid_resolution);

  Color black = {{ 0, 0, 0 }};
  
  for (std::size_t i = 0; i < m_psc->HPS.size(); ++ i)
    {
      m_psc->HPS[i].color = colors[i];
      if (colors[i] != black)
        m_psc->has_colors = true;
      m_psc->HPS[i].echo = echo[i];
      if (echo[i] != 0)
        m_psc->is_echo_given = true;
    }

  timer.stop();
  std::cerr << "Reading took "<< timer.time() << "sec."<<std::endl;
  if (m_psc->has_colors)
    std::cerr << "Has colors" << std::endl;

  invalidateOpenGLBuffers();
  return true;
}

// Write point set to .PLY file
bool Scene_point_set_classification_item::write_ply_point_set(std::ostream& stream) const
{
  if (m_psc->HPS[0].AE_label == (unsigned char)(-1))
    {
      std::cerr << "Error: classification was not performed." << std::endl;
      return false;
    }
    
  stream.precision (std::numeric_limits<double>::digits10 + 2);
  
  stream << "ply" << std::endl
         << "format ascii 1.0" << std::endl
         << "element vertex " << m_psc->HPS.size() << std::endl
         << "property float x" << std::endl
         << "property float y" << std::endl
         << "property float z" << std::endl
         << "property uchar red" << std::endl
         << "property uchar green" << std::endl
         << "property uchar blue" << std::endl
         << "end_header" << std::endl;

  for (std::size_t i = 0; i < m_psc->HPS.size(); ++ i)
    {
      QColor color (0, 0, 0);
      int c = m_psc->HPS[i].AE_label;
          
      if (c != -1)
        {        
          if (m_psc->segmentation_classes[c]->id() == "vegetation")
            color = QColor(0, 255, 27);
          else if (m_psc->segmentation_classes[c]->id() == "ground")
            color = QColor(245, 180, 0);
          else if (m_psc->segmentation_classes[c]->id() == "road")
            color = QColor(114, 114, 130);
          else if (m_psc->segmentation_classes[c]->id() == "roof")
            color = QColor(255, 0, 170);
          else if (m_psc->segmentation_classes[c]->id() == "facade")
            color = QColor(100, 0, 255);
          else if (m_psc->segmentation_classes[c]->id() == "building")
                      color = QColor(0, 114, 225);
        }
      
      stream << m_psc->HPS[i].position << " "
             << color.red() << " "
             << color.green() << " "
             << color.blue() << std::endl;
    }
    
  return true;
}

QString
Scene_point_set_classification_item::toolTip() const
{
  return QObject::tr("<p><b>%1</b><br />"
                     "<i>Point set classification</i></p>"
                     "<p>Number of points: %2</p>")
    .arg(name())
    .arg(m_psc->HPS.size());
}

bool Scene_point_set_classification_item::supportsRenderingMode(RenderingMode m) const 
{
  return m==PointsPlusNormals;

}

void Scene_point_set_classification_item::draw_splats(CGAL::Three::Viewer_interface*) const
{
}

void Scene_point_set_classification_item::draw_edges(CGAL::Three::Viewer_interface* viewer) const
{
  if(!are_buffers_filled)
    initialize_buffers(viewer);
  vaos[Edges]->bind();
  program=getShaderProgram(PROGRAM_NO_SELECTION);
  attrib_buffers(viewer,PROGRAM_NO_SELECTION);
  program->bind();
  program->setAttributeValue("colors", this->color());
  viewer->glDrawArrays(GL_LINES, 0,
                       static_cast<GLsizei>((nb_lines)/3));
  vaos[Edges]->release();
  program->release();
}

void Scene_point_set_classification_item::draw_points(CGAL::Three::Viewer_interface* viewer) const
{
  if(!are_buffers_filled)
    initialize_buffers(viewer);
  GLfloat point_size;
  viewer->glGetFloatv(GL_POINT_SIZE, &point_size);
  if (m_index_color == 9)
    viewer->glPointSize(15.f);
  else
    viewer->glPointSize(3.f);

  vaos[ThePoints]->bind();
  program=getShaderProgram(PROGRAM_NO_SELECTION);
  attrib_buffers(viewer,PROGRAM_NO_SELECTION);
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
  if (m_psc == NULL)
    return;
  
  double xmin = std::numeric_limits<double>::max();
  double ymin = std::numeric_limits<double>::max();
  double zmin = std::numeric_limits<double>::max();
  double xmax = -std::numeric_limits<double>::max();
  double ymax = -std::numeric_limits<double>::max();
  double zmax = -std::numeric_limits<double>::max();

  for (std::size_t i = 0; i < m_psc->HPS.size(); ++ i)
    {
      const Kernel::Point_3& p = m_psc->HPS[i].position;
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
  
  if (out == 1 && m_psc->HPS[0].AE_label == (unsigned char)(-1))
    out = 0;

  if (out == 3 && m_scat == NULL)
    out = 0;
  if (out == 4 && m_plan == NULL)
    out = 0;
  if (out == 5 && m_hori == NULL)
    out = 0;
  if (out == 6 && m_elev == NULL)
    out = 0;
  if (out == 7 && m_colo == NULL)
    out = 0;

  if (out == 0 && !(m_psc->has_colors))
    out = -1;
  return out;
}

void Scene_point_set_classification_item::estimate_parameters (double& grid_resolution,
                                                               double& radius_neighbors,
                                                               double& radius_dtm)
{
  std::vector<Kernel::Point_3> pts;
  for (std::size_t i = 0; i < m_psc->HPS.size(); ++ i)
    pts.push_back(m_psc->HPS[i].position);

  double average_spacing = CGAL::compute_average_spacing<CGAL::Sequential_tag> (pts.begin(), pts.end(), 6);

  grid_resolution = average_spacing;
  radius_neighbors = 5 * grid_resolution;
  radius_dtm = 5 * radius_neighbors;
}

void Scene_point_set_classification_item::compute_features (const double& grid_resolution,
                                                            const double& radius_neighbors,
                                                            const double& radius_dtm,
                                                            const QColor& c)
{
  Q_ASSERT (m_psc != NULL);

  bool on_groups = !(m_psc->groups.empty());
  
  m_psc->set_parameters (grid_resolution, radius_neighbors, radius_dtm);
  m_psc->initialization();

  if (m_scat != NULL) delete m_scat;
  m_scat = new Scatter (*m_psc, 1., on_groups);
  m_psc->segmentation_attributes.push_back (m_scat);
  
  if (m_plan != NULL) delete m_plan;
  m_plan = new NonPlanarity (*m_psc, 1., on_groups);
  m_psc->segmentation_attributes.push_back (m_plan);
  
  if (m_hori != NULL) delete m_hori;
  m_hori = new Horizontality (*m_psc, 1., on_groups);
  m_psc->segmentation_attributes.push_back (m_hori);

  if (m_elev != NULL) delete m_elev;
  m_elev = new Elevation (*m_psc, 1., on_groups);
  m_psc->segmentation_attributes.push_back (m_elev);

  if (m_colo != NULL) delete m_colo;
  m_colo = new ColorSeg (*m_psc, 1., c.hue(), c.saturation(), c.value());
  m_psc->segmentation_attributes.push_back (m_colo);

}


void Scene_point_set_classification_item::compute_ransac (const double& radius_neighbors)
{
  Q_ASSERT (m_psc != NULL);

  std::cerr << "Computing normals..." << std::endl;

  // Estimates normals direction.
  std::vector<std::size_t> indices (m_psc->HPS.size());
  for (std::size_t i = 0; i < m_psc->HPS.size(); ++ i)
    indices[i] = i;

  std::vector<Kernel::Vector_3> normals (m_psc->HPS.size(), CGAL::NULL_VECTOR);

  HPS_property_map<PSC::HPoint> hps_pmap (&(m_psc->HPS));
  CGAL::jet_estimate_normals<CGAL::Sequential_tag>(indices.begin(), indices.end(),
                                                   hps_pmap,
                                                   &normals[0],
                                                   12);

  typedef CGAL::Shape_detection_3::Efficient_RANSAC_traits<Kernel, std::vector<std::size_t>,
                                                           HPS_property_map<PSC::HPoint>,
                                                           const Kernel::Vector_3* > Traits;
  typedef CGAL::Shape_detection_3::Efficient_RANSAC<Traits> Shape_detection;

  Shape_detection shape_detection;


  shape_detection.set_input(indices, hps_pmap, &(normals[0]));
  //  shape_detection.add_shape_factory<CGAL::Shape_detection_3::Plane<Traits> >();
  shape_detection.add_shape_factory<My_plane<Traits> >();
  Shape_detection::Parameters op;
  op.probability = 0.05;
  op.min_points = std::max ((std::size_t)10, m_psc->HPS.size() / 250000);
  
  op.epsilon = radius_neighbors;
  op.cluster_epsilon = radius_neighbors;
  op.normal_threshold = 0.9;

  std::cerr << "Computing RANSAC..." << std::endl;
  shape_detection.detect(op);

  
  CGAL::regularize_planes (shape_detection, true, true, true, false,
                           10., radius_neighbors);
  
  m_psc->reset_groups();

  std::cerr << "Storing groups..." << std::endl;
  
  BOOST_FOREACH(boost::shared_ptr<Shape_detection::Shape> shape, shape_detection.shapes())
    {
      m_psc->groups.push_back ((Kernel::Plane_3)(*(dynamic_cast<CGAL::Shape_detection_3::Plane<Traits>*>(shape.get ()))));
      BOOST_FOREACH(std::size_t i, shape->indices_of_assigned_points())
        m_psc->HPS[indices[i]].group = m_psc->groups.size() - 1;
      
    }
  std::cerr << "Found " << m_psc->groups.size() << " group(s)" << std::endl;
}


void Scene_point_set_classification_item::compute_clusters (const double& radius_neighbors)
{
  Q_ASSERT (m_psc != NULL);

  m_psc->cluster_points (radius_neighbors);

  std::size_t index_ground = 0;
  std::size_t index_roof = 0;
  std::size_t index_facade = 0;
  std::size_t index_vege = 0;
  
  for (std::size_t i = 0; i < m_psc->segmentation_classes.size(); ++ i)
    if (m_psc->segmentation_classes[i]->id() == "ground")
      index_ground = i;
    else if (m_psc->segmentation_classes[i]->id() == "roof")
      index_roof = i;
    else if (m_psc->segmentation_classes[i]->id() == "facade")
      index_facade = i;
    else if (m_psc->segmentation_classes[i]->id() == "vegetation")
      index_vege = i;

  for (std::size_t i = 0; i < m_psc->HPS.size(); ++ i)
    {
      if (m_psc->HPS[i].neighbor == (unsigned char)(-1))
        {
          m_psc->HPS[i].confidence = 0.;
          continue;
        }

      std::size_t c0 = m_psc->HPS[i].AE_label;
      std::size_t c1 = m_psc->HPS[i].neighbor;
      
      if ((c0 == index_ground && c1 == index_facade) || (c0 == index_facade && c1 == index_ground) ||
          (c0 == index_ground && c1 == index_vege) || (c0 == index_vege && c1 == index_ground) ||
          (c0 == index_facade && c1 == index_roof) || (c0 == index_roof && c1 == index_facade))
        m_psc->HPS[i].confidence = 1.;
      else
        m_psc->HPS[i].confidence = 0.2;
    }
  
  // for (std::size_t i = 0; i < m_psc->clusters.size(); ++ i)
  //   {
  //     std::size_t c0 = m_psc->HPS[m_psc->clusters[i].indices[0]].AE_label;
  //     std::size_t nb_good = 1;
  //     std::size_t nb_tot = 1;
      
  //     for (std::set<std::size_t>::iterator it = m_psc->clusters[i].neighbors.begin ();
  //          it != m_psc->clusters[i].neighbors.end (); ++ it)
  //       {
  //         std::size_t c1 = m_psc->HPS[m_psc->clusters[*it].indices[0]].AE_label;
  //         std::size_t weight = 1;//std::min (m_psc->clusters[*it].indices.size(), m_psc->clusters[i].indices.size());
  //         nb_tot += weight;
          
  //         if ((c0 == index_ground && c1 == index_facade) || (c0 == index_facade && c1 == index_ground) ||
  //             (c0 == index_ground && c1 == index_vege) || (c0 == index_vege && c1 == index_ground) ||
  //             (c0 == index_facade && c1 == index_roof) || (c0 == index_roof && c1 == index_facade))
  //           nb_good += weight;
  //       }
      
  //     double confidence = nb_good / (double)nb_tot;
  //     for (std::size_t j = 0; j < m_psc->clusters[i].indices.size(); ++ j)
  //       m_psc->HPS[m_psc->clusters[i].indices[j]].confidence = confidence;
  //   }

  invalidateOpenGLBuffers();
}


void Scene_point_set_classification_item::save_2d_image()
{
  std::ofstream f ("out.ppm");
  f << "P3" << std::endl
    << m_psc->grid_HPS.width() << " " << m_psc->grid_HPS.height() << std::endl
    << "255" << std::endl;

  std::size_t cnt = 0;
  for (std::size_t i = 0; i < m_psc->grid_HPS.height(); ++ i)
    for (std::size_t j = 0; j < m_psc->grid_HPS.width(); ++ j)
      {
        std::vector<std::size_t> scores (m_psc->segmentation_classes.size(), 0);

        if (m_psc->grid_HPS(i,j).empty())
          f << "255 255 255 ";
        else
          {
            for (std::size_t k = 0; k < m_psc->grid_HPS(i,j).size(); ++ k)
              scores[m_psc->HPS[m_psc->grid_HPS(i,j)[k]].AE_label] ++;

            std::size_t max = 0;
            std::size_t c = 0;
            for (std::size_t k = 0; k < scores.size(); ++ k)
              if (scores[k] > max)
                {
                  c = k;
                  max = scores[k];
                }
        
            if (m_psc->segmentation_classes[c]->id() == "vegetation")
              f << "0 255 27 ";
            else if (m_psc->segmentation_classes[c]->id() == "ground")
              f << "245 180 0 ";
            else if (m_psc->segmentation_classes[c]->id() == "road")
              f << "114 114 130 ";
            else if (m_psc->segmentation_classes[c]->id() == "roof")
              f << "255 0 170 ";
            else if (m_psc->segmentation_classes[c]->id() == "facade")
              f << "100 0 255 ";
            else if (m_psc->segmentation_classes[c]->id() == "building")
              f << "0 114 255 ";
            else
              f << "255 255 255 ";
          }
        
        ++ cnt;
        if (cnt == 5)
          {
            f << std::endl;
            cnt = 0;
          }
      }

  f.close();
}

void Scene_point_set_classification_item::extract_2d_outline (double radius,
                                                              std::vector<Kernel::Point_3>& outline)
{
    
  // Projection
  double mean_z = 0.;
  std::size_t nb = 0;
  std::vector<Kernel::Point_2> pts;
  
  for (std::size_t i = 0; i < m_psc->HPS.size(); ++ i)
    if (m_psc->HPS[i].AE_label == 1)
      {
        mean_z += m_psc->HPS[i].position.z();
        nb ++;
      }
    else if (m_psc->HPS[i].AE_label > 1)
      pts.push_back (Kernel::Point_2 (m_psc->HPS[i].position.x(),
                                      m_psc->HPS[i].position.y()));

  mean_z /= nb;
  

  // Alpha shapes
  typedef CGAL::Alpha_shape_vertex_base_2<Kernel> Vb;
  typedef CGAL::Alpha_shape_face_base_2<Kernel>  Fb;
  typedef CGAL::Triangulation_data_structure_2<Vb,Fb> Tds;
  typedef CGAL::Delaunay_triangulation_2<Kernel,Tds> Triangulation_2;
  typedef CGAL::Alpha_shape_2<Triangulation_2>  Alpha_shape_2;

  Alpha_shape_2 ashape (pts.begin (), pts.end (), radius);

  // Output
  for (Alpha_shape_2::Finite_edges_iterator it = ashape.finite_edges_begin ();
       it != ashape.finite_edges_end (); ++ it)
    {
      if (ashape.classify(*it) != Alpha_shape_2::REGULAR)
        continue;

      for (std::size_t i = 0; i < 2; ++ i)
        {
          const Kernel::Point_2& p = it->first->vertex ((it->second + 1 + i) % 3)->point();
          outline.push_back (Kernel::Point_3 (p.x(), p.y(), mean_z));
        }

    }

}

void Scene_point_set_classification_item::extract_building_map (double,
                                                                std::vector<Kernel::Triangle_3>& faces)
{
  std::size_t index_ground = (std::size_t)(-1);
  std::size_t index_roof = (std::size_t)(-1);
  std::size_t index_facade = (std::size_t)(-1);
  for (std::size_t i = 0; i < m_psc->segmentation_classes.size(); ++ i)
    if (m_psc->segmentation_classes[i]->id() == "ground")
      index_ground = i;
    else if (m_psc->segmentation_classes[i]->id() == "roof")
      index_roof = i;
    else if (m_psc->segmentation_classes[i]->id() == "facade")
      index_facade = i;

  bool separate_only = (index_roof == (std::size_t)(-1));

  // Estimate ground Z
  std::cerr << "Estimating ground..." << std::endl;
  std::vector<double> z_ground;
  for (std::size_t i = 0; i < m_psc->HPS.size(); ++ i)
    if (m_psc->HPS[i].AE_label == index_ground)
      z_ground.push_back (m_psc->HPS[i].position.z());

  std::sort (z_ground.begin(), z_ground.end());
  double z_median = z_ground[z_ground.size() / 2];
  Kernel::Plane_3 ground (Kernel::Point_3 (0., 0., z_median), Kernel::Vector_3 (0., 0., 1.));
  
  // Project facade planes and get lines
  std::cerr << "Projecting facade planes..." << std::endl;
  std::vector<bool> is_plane_facade (m_psc->groups.size(), false);
  std::vector<std::vector<Kernel::Point_3> > facade_pts (m_psc->groups.size());
  
  for (std::size_t i = 0; i < m_psc->HPS.size(); ++ i)
    if (m_psc->HPS[i].AE_label == index_facade && m_psc->HPS[i].group != (std::size_t)(-1))
      {
        is_plane_facade[m_psc->HPS[i].group] = true;
        facade_pts[m_psc->HPS[i].group].push_back (m_psc->HPS[i].position);
      }


  std::vector<Kernel::Segment_3> borders;
  borders.push_back (Kernel::Segment_3 (Kernel::Point_3 (_bbox.xmin, _bbox.ymin, z_median),
                                        Kernel::Point_3 (_bbox.xmin, _bbox.ymax, z_median)));
  borders.push_back (Kernel::Segment_3 (Kernel::Point_3 (_bbox.xmin, _bbox.ymax, z_median),
                                        Kernel::Point_3 (_bbox.xmax, _bbox.ymax, z_median)));
  borders.push_back (Kernel::Segment_3 (Kernel::Point_3 (_bbox.xmax, _bbox.ymax, z_median),
                                        Kernel::Point_3 (_bbox.xmax, _bbox.ymin, z_median)));
  borders.push_back (Kernel::Segment_3 (Kernel::Point_3 (_bbox.xmax, _bbox.ymin, z_median),
                                        Kernel::Point_3 (_bbox.xmin, _bbox.ymin, z_median)));

  std::vector<Kernel::Segment_3> lines;
  for (std::size_t i = 0; i < is_plane_facade.size(); ++ i)
    if (is_plane_facade[i])
      {
        Kernel::Plane_3 plane = m_psc->groups[i];
        Kernel::Vector_3 normal = plane.orthogonal_vector();
        normal = normal / std::sqrt (normal * normal);
        double horiz = normal * Kernel::Vector_3 (0., 0., 1.);
        // if (horiz > 0.17 || horiz < -0.17)
        //   continue;

        Kernel::Point_3 c = CGAL::centroid (facade_pts[i].begin(), facade_pts[i].end());
        normal = Kernel::Vector_3 (normal.x(), normal.y(), 0.);
        plane = Kernel::Plane_3 (c, normal);
        
        //#define INFINITE_LINES
#if defined(INFINITE_LINES)

        Kernel::Point_3 source, target;
        bool source_found = false;

        for (std::size_t j = 0; j < borders.size(); ++ j)
          {
            CGAL::Object obj = CGAL::intersection (borders[j], plane);
            Kernel::Point_3 inter;
            if (CGAL::assign (inter, obj))
              {
                if (source_found)
                  {
                    target = inter;
                    lines.push_back (Kernel::Segment_3 (source, target));
                    break;
                  }
                else
                  {
                    source = inter;
                    source_found = true;
                  }
              }
          }
#else
        CGAL::Object obj = CGAL::intersection (ground, plane);
        Kernel::Line_3 line;

        // if (facade_pts[i].size() < 100)
        //   continue;
        
        if (CGAL::assign (line, obj))
          {
            Kernel::Point_3 origin = line.projection (CGAL::ORIGIN);
            Kernel::Vector_3 v (line);
            double min = std::numeric_limits<double>::max();
            double max = -std::numeric_limits<double>::max();
            Kernel::Point_3 source = CGAL::ORIGIN;
            Kernel::Point_3 target = CGAL::ORIGIN;
            for (std::size_t j = 0; j < facade_pts[i].size(); ++ j)
              {
                Kernel::Vector_3 dir (origin, facade_pts[i][j]);
                double prod = v * dir;
                if (prod > max)
                  {
                    max = prod;
                    target = facade_pts[i][j];
                  }
                if (prod < min)
                  {
                    min = prod;
                    source = facade_pts[i][j];
                  }
              }
            source = line.projection (source);
            target = line.projection (target);
            Kernel::Point_3 middle (0.5 * (source.x() + target.x()),
                                    0.5 * (source.y() + target.y()),
                                    0.5 * (source.z() + target.z()));

            source = middle + 1.1 * (source - middle);
            target = middle + 1.1 * (target - middle);
            lines.push_back (Kernel::Segment_3 (source, target));
          }
#endif
      }

  std::cerr << " -> Found " << lines.size() << " line(s)" << std::endl;
  
  // Build 2D constrained Delaunay triangulation
  std::cerr << "Building 2D constrained Delaunay triangulation..." << std::endl;
    
  typedef CGAL::Triangulation_vertex_base_2<Kernel>                Vb;
  typedef CGAL::Triangulation_face_base_with_info_2<int,Kernel>    Fbwi;
  typedef CGAL::Constrained_triangulation_face_base_2<Kernel,Fbwi> Fb;
  typedef CGAL::Triangulation_data_structure_2<Vb,Fb>              TDS;
  typedef CGAL::Exact_predicates_tag                               Itag;
  typedef CGAL::Constrained_Delaunay_triangulation_2<Kernel, TDS, Itag> CDT;

  CDT cdt;
  cdt.insert (Kernel::Point_2 (_bbox.xmin, _bbox.ymin));
  cdt.insert (Kernel::Point_2 (_bbox.xmin, _bbox.ymax));
  cdt.insert (Kernel::Point_2 (_bbox.xmax, _bbox.ymin));
  cdt.insert (Kernel::Point_2 (_bbox.xmax, _bbox.ymax));
  for (std::size_t i = 0; i < lines.size(); ++ i)
    cdt.insert_constraint (Kernel::Point_2 (lines[i].source().x(), lines[i].source().y()),
                           Kernel::Point_2 (lines[i].target().x(), lines[i].target().y()));

  std::cerr << " -> " << cdt.number_of_faces () << " face(s) created" << std::endl;
  if (separate_only)
    {
      std::vector<std::vector<double> > heights;
      
      for (CDT::Finite_faces_iterator it = cdt.finite_faces_begin(); it != cdt.finite_faces_end(); ++ it)
        {
          it->info() = heights.size();
          heights.push_back (std::vector<double>());
          
        }

      for (std::size_t i = 0; i < m_psc->HPS.size(); ++ i)
        if (m_psc->HPS[i].AE_label == index_ground)
          {
            CDT::Face_handle f = cdt.locate (Kernel::Point_2 (m_psc->HPS[i].position.x(),
                                                              m_psc->HPS[i].position.y()));

            heights[f->info()].push_back (m_psc->HPS[i].position.z());
          }
      
      std::vector<double> h;
      for (std::size_t i = 0; i < heights.size(); ++ i)
        {
          if (heights[i].empty())
            {
              h.push_back (0.);
              continue;
            }
          std::sort (heights[i].begin(), heights[i].end());
          h.push_back (heights[i][heights[i].size() / 10]);
        }
      
      for (CDT::Finite_faces_iterator it = cdt.finite_faces_begin(); it != cdt.finite_faces_end(); ++ it)
        {
          faces.push_back (Kernel::Triangle_3 (Kernel::Point_3 (it->vertex (0)->point().x(),
                                                                it->vertex (0)->point().y(),
                                                                h[it->info()]),
                                               Kernel::Point_3 (it->vertex (1)->point().x(),
                                                                it->vertex (1)->point().y(),
                                                                h[it->info()]),
                                               Kernel::Point_3 (it->vertex (2)->point().x(),
                                                                it->vertex (2)->point().y(),
                                                                h[it->info()])));
        }

      for (CDT::Finite_edges_iterator it = cdt.finite_edges_begin(); it != cdt.finite_edges_end(); ++ it)
        {
          if (cdt.is_infinite (it->first) ||
              cdt.is_infinite (it->first->neighbor(it->second)))
            continue;
          Kernel::Point_3 a (it->first->vertex ((it->second + 1)%3)->point().x(),
                             it->first->vertex ((it->second + 1)%3)->point().y(),
                             h[it->first->info()]);
          Kernel::Point_3 b (it->first->vertex ((it->second + 1)%3)->point().x(),
                             it->first->vertex ((it->second + 1)%3)->point().y(),
                             h[it->first->neighbor(it->second)->info()]);
          Kernel::Point_3 c (it->first->vertex ((it->second + 2)%3)->point().x(),
                             it->first->vertex ((it->second + 2)%3)->point().y(),
                             h[it->first->info()]);
          Kernel::Point_3 d (it->first->vertex ((it->second + 2)%3)->point().x(),
                             it->first->vertex ((it->second + 2)%3)->point().y(),
                             h[it->first->neighbor(it->second)->info()]);
          faces.push_back (Kernel::Triangle_3 (a, b, c));
          faces.push_back (Kernel::Triangle_3 (b, c, d));
        }
    }
  else
    {
      for (CDT::Finite_faces_iterator it = cdt.finite_faces_begin(); it != cdt.finite_faces_end(); ++ it)
        it->info() = -10;

      // Project ground + roof points on faces
      std::cerr << "Projecting ground and roof points on faces..." << std::endl;
      for (std::size_t i = 0; i < m_psc->HPS.size(); ++ i)
        {
          int iter = 0;
          if (m_psc->HPS[i].AE_label == index_roof)
            iter = 1;
          else if (m_psc->HPS[i].AE_label == index_ground)
            iter = -1;
          else
            continue;

          CDT::Face_handle f = cdt.locate (Kernel::Point_2 (m_psc->HPS[i].position.x(),
                                                            m_psc->HPS[i].position.y()));

          f->info() += iter;
        }

  
      // Label and extract faces
      for (CDT::Finite_faces_iterator it = cdt.finite_faces_begin(); it != cdt.finite_faces_end(); ++ it)
        if (it->info() > 0)
          faces.push_back (Kernel::Triangle_3 (Kernel::Point_3 (it->vertex (0)->point().x(),
                                                                it->vertex (0)->point().y(),
                                                                z_median),
                                               Kernel::Point_3 (it->vertex (1)->point().x(),
                                                                it->vertex (1)->point().y(),
                                                                z_median),
                                               Kernel::Point_3 (it->vertex (2)->point().x(),
                                                                it->vertex (2)->point().y(),
                                                                z_median)));
      std::cerr << " -> Found " << faces.size() << " building face(s)" << std::endl;
    }
}


void Scene_point_set_classification_item::extract_facades (double radius,
                                                           std::vector<Kernel::Triangle_3>& faces)
{
  std::size_t index_facade = 0;
  for (std::size_t i = 0; i < m_psc->segmentation_classes.size(); ++ i)
    if (m_psc->segmentation_classes[i]->id() == "facade")
      {
        index_facade = i;
        break;
      }

  typedef CGAL::Alpha_shape_vertex_base_2<Kernel> Vb;
  typedef CGAL::Alpha_shape_face_base_2<Kernel>  Fb;
  typedef CGAL::Triangulation_data_structure_2<Vb,Fb> Tds;
  typedef CGAL::Delaunay_triangulation_2<Kernel,Tds> Triangulation_2;
  typedef CGAL::Alpha_shape_2<Triangulation_2>  Alpha_shape_2;
  
  std::vector<std::vector<Kernel::Point_3> > facade_pts (m_psc->groups.size());
    
  for (std::size_t i = 0; i < m_psc->HPS.size(); ++ i)
    if (m_psc->HPS[i].AE_label == index_facade && m_psc->HPS[i].group != (std::size_t)(-1))
      {
        facade_pts[m_psc->HPS[i].group].push_back (m_psc->HPS[i].position);
      }

  for (std::size_t i = 0; i < facade_pts.size(); ++ i)
    {
      Kernel::Plane_3 plane = m_psc->groups[i];
      Kernel::Vector_3 normal = plane.orthogonal_vector();
      Kernel::Point_3 c = CGAL::centroid (facade_pts[i].begin(), facade_pts[i].end());
      normal = Kernel::Vector_3 (normal.x(), normal.y(), 0.);
      plane = Kernel::Plane_3 (c, normal);

      std::vector<Kernel::Point_2> projections;
      projections.reserve (facade_pts[i].size ());

      for (std::size_t j = 0; j < facade_pts[i].size (); ++ j)
        projections.push_back (plane.to_2d (facade_pts[i][j]));

      Alpha_shape_2 ashape (projections.begin (), projections.end (), radius);
  
      for (Alpha_shape_2::Finite_faces_iterator it = ashape.finite_faces_begin ();
           it != ashape.finite_faces_end (); ++ it)
        {
          if (ashape.classify (it) != Alpha_shape_2::INTERIOR)
            continue;
          faces.push_back (Kernel::Triangle_3 (plane.to_3d (it->vertex(0)->point()),
                                               plane.to_3d (it->vertex(1)->point()),
                                               plane.to_3d (it->vertex(2)->point())));
        }
    }
}
