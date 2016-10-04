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
  m_index_color = 1;

  m_points->point_set()->unselect_all();
  m_points->point_set()->collect_garbage();
  m_psc = new PSC(m_points->point_set()->begin(), m_points->point_set()->end(), m_points->point_set()->point_map());
  
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
  m_index_color = 1;
  
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
      if (m_predefined_types.empty())
        {
          for (Point_set::const_iterator it = m_points->point_set()->begin();
               it != m_points->point_set()->first_selected(); ++ it)
            {
              QColor color (0, 0, 0);
              Type_handle c = m_psc->classification_type_of(*it);
          
              if (c != Type_handle())
                {        
                  if (c->id() == "vegetation")
                    color = QColor(0, 255, 27);
                  else if (c->id() == "ground")
                    color = QColor(245, 180, 0);
                  else if (c->id() == "road")
                    color = QColor(114, 114, 130);
                  else if (c->id() == "roof")
                    color = QColor(255, 0, 170);
                  else if (c->id() == "facade")
                    color = QColor(100, 0, 255);
                  else if (c->id() == "building")
                    color = QColor(0, 114, 225);
                }

              colors_points.push_back ((double)(color.red()) / 255.);
              colors_points.push_back ((double)(color.green()) / 255.);
              colors_points.push_back ((double)(color.blue()) / 255.);
            }
        }
      else
        {
          std::map<Type_handle, QColor> map_colors;
          for (std::size_t i = 0; i < m_predefined_types.size(); ++ i)
            map_colors.insert (m_predefined_types[i]);
          
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
    }
  else if (index_color == 2) // confidence
    {
      for (Point_set::const_iterator it = m_points->point_set()->begin();
           it != m_points->point_set()->first_selected(); ++ it)
        {
          QColor color (0, 0, 0);
          Type_handle c = m_psc->classification_type_of(*it);
          
          if (c != Type_handle())
            {        
              if (c->id() == "vegetation")
                color = QColor(0, 255, 27);
              else if (c->id() == "ground")
                color = QColor(245, 180, 0);
              else if (c->id() == "road")
                color = QColor(114, 114, 130);
              else if (c->id() == "roof")
                color = QColor(255, 0, 170);
              else if (c->id() == "facade")
                color = QColor(100, 0, 255);
              else if (c->id() == "building")
                color = QColor(0, 114, 225);
            }

          double confidence = m_psc->confidence_of(*it);
          colors_points.push_back (1. - confidence * (1. - (double)(color.red()) / 255.));
          colors_points.push_back (1. - confidence * (1. - (double)(color.green()) / 255.));
          colors_points.push_back (1. - confidence * (1. - (double)(color.blue()) / 255.));
        }
    }
  else if (index_color == 3) // scatter
    {
      double weight = m_disp->weight;
      m_disp->weight = m_disp->max;
      for (Point_set::const_iterator it = m_points->point_set()->begin();
           it != m_points->point_set()->first_selected(); ++ it)
        {
          colors_points.push_back (ramp.r(m_disp->value(*it)));
          colors_points.push_back (ramp.g(m_disp->value(*it)));
          colors_points.push_back (ramp.b(m_disp->value(*it)));
        }
      m_disp->weight = weight;
    }
  else if (index_color == 4) // planarity
    {
      double weight = m_d2p->weight;
      m_d2p->weight = m_d2p->max;
      for (Point_set::const_iterator it = m_points->point_set()->begin();
           it != m_points->point_set()->first_selected(); ++ it)
        {
          colors_points.push_back (ramp.r(m_d2p->value(*it)));
          colors_points.push_back (ramp.g(m_d2p->value(*it)));
          colors_points.push_back (ramp.b(m_d2p->value(*it)));
        }
      m_d2p->weight = weight;
    }
  else if (index_color == 5) // horizontality
    {
      double weight = m_verti->weight;
      m_verti->weight = m_verti->max;
      for (Point_set::const_iterator it = m_points->point_set()->begin();
           it != m_points->point_set()->first_selected(); ++ it)
        {
          colors_points.push_back (ramp.r(m_verti->value(*it)));
          colors_points.push_back (ramp.g(m_verti->value(*it)));
          colors_points.push_back (ramp.b(m_verti->value(*it)));
        }
      m_verti->weight = weight;
    }
  else if (index_color == 6) // elevation
    {
      double weight = m_elev->weight;
      m_elev->weight = m_elev->max;
      for (Point_set::const_iterator it = m_points->point_set()->begin();
           it != m_points->point_set()->first_selected(); ++ it)
        {
          double v = std::max(0., m_elev->value(*it));
          colors_points.push_back (ramp.r(v));
          colors_points.push_back (ramp.g(v));
          colors_points.push_back (ramp.b(v));
        }
      m_elev->weight = weight;
    }
  else if (index_color == 7) // color
    {
      double weight = m_col_att->weight;
      m_col_att->weight = m_col_att->max;
      for (Point_set::const_iterator it = m_points->point_set()->begin();
           it != m_points->point_set()->first_selected(); ++ it)
        {
          colors_points.push_back (ramp.r(m_col_att->value(*it)));
          colors_points.push_back (ramp.g(m_col_att->value(*it)));
          colors_points.push_back (ramp.b(m_col_att->value(*it)));
        }
      m_col_att->weight = weight;
    }
  else if (index_color == 8) // RANSAC
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
  else if (index_color == 9) // Clusters
    {
      // for (std::size_t i = 0; i < m_psc->clusters().size(); ++ i)
      //   {
      //     QColor color (0, 0, 0);
      //     CGAL::Data_classification::Type *c = m_psc->classification_type_of(m_psc->clusters()[i].indices[0]);
          
      //     if (c != NULL)
      //       {        
      //         if (c->id() == "vegetation")
      //           color = QColor(0, 255, 27);
      //         else if (c->id() == "ground")
      //           color = QColor(245, 180, 0);
      //         else if (c->id() == "road")
      //           color = QColor(114, 114, 130);
      //         else if (c->id() == "roof")
      //           color = QColor(255, 0, 170);
      //         else if (c->id() == "facade")
      //           color = QColor(100, 0, 255);
      //         else if (c->id() == "building")
      //           color = QColor(0, 114, 225);
      //       }


      //     colors_points.push_back ((double)(color.red()) / 255.);
      //     colors_points.push_back ((double)(color.green()) / 255.);
      //     colors_points.push_back ((double)(color.blue()) / 255.);
      //   }
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
bool Scene_point_set_classification_item::write_ply_point_set(std::ostream& stream) const
{
  if (m_psc->classification_type_of(0) == Type_handle())
    {
      std::cerr << "Error: classification was not performed." << std::endl;
      return false;
    }
    
  stream.precision (std::numeric_limits<double>::digits10 + 2);
  
  stream << "ply" << std::endl
         << "format ascii 1.0" << std::endl
         << "element vertex " << m_points->point_set()->size() << std::endl
         << "property float x" << std::endl
         << "property float y" << std::endl
         << "property float z" << std::endl
         << "property uchar red" << std::endl
         << "property uchar green" << std::endl
         << "property uchar blue" << std::endl
         << "end_header" << std::endl;

  for (Point_set::const_iterator it = m_points->point_set()->begin();
       it != m_points->point_set()->end(); ++ it)
    {
      QColor color (0, 0, 0);
      Type_handle c = m_psc->classification_type_of(m_psc->clusters()[*it].indices[0]);
          
      if (c != Type_handle())
        {        
          if (c->id() == "vegetation")
            color = QColor(0, 255, 27);
          else if (c->id() == "ground")
            color = QColor(245, 180, 0);
          else if (c->id() == "road")
            color = QColor(114, 114, 130);
          else if (c->id() == "roof")
            color = QColor(255, 0, 170);
          else if (c->id() == "facade")
            color = QColor(100, 0, 255);
          else if (c->id() == "building")
                      color = QColor(0, 114, 225);
        }
      
      stream << m_points->point_set()->point(*it) << " "
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
  
  if (out == 3 && m_disp == Attribute_handle())
    out = 0;
  if (out == 4 && m_d2p == Attribute_handle())
    out = 0;
  if (out == 5 && m_verti == Attribute_handle())
    out = 0;
  if (out == 6 && m_elev == Attribute_handle())
    out = 0;
  if (out == 7 && m_col_att == Attribute_handle())
    out = 0;

  if (out == 0 && !(m_points->point_set()->has_colors()))
    out = -1;
  return out;
}

void Scene_point_set_classification_item::estimate_parameters (double& grid_resolution,
                                                               double& radius_neighbors,
                                                               double& radius_dtm)
{
  double average_spacing = CGAL::compute_average_spacing<CGAL::Sequential_tag> (m_points->point_set()->begin(),
                                                                                m_points->point_set()->end(),
                                                                                m_points->point_set()->point_map(), 6);

  grid_resolution = average_spacing;
  radius_neighbors = 5 * grid_resolution;
  radius_dtm = 5 * radius_neighbors;
}

void Scene_point_set_classification_item::compute_features (const double& grid_resolution,
                                                            const double& radius_neighbors,
                                                            const double& radius_dtm,
                                                            const QColor& c)
{
  Q_ASSERT (!(m_points->point_set()->empty()));
  Q_ASSERT (m_psc != NULL);
  m_psc->clear_attributes();

  compute_bbox();
  if (m_helper != NULL) delete m_helper;
  m_helper = new Helper (m_points->point_set()->begin(),
                         m_points->point_set()->end(),
                         m_points->point_set()->point_map(),
                         grid_resolution, radius_neighbors, radius_dtm);

  
  m_disp = Attribute_handle (new Dispersion (m_points->point_set()->begin(),
                           m_points->point_set()->end(),
                           m_points->point_set()->point_map(),
                           m_helper->grid(),
                           grid_resolution,
                           radius_neighbors,
                           1.));
  m_psc->add_attribute (m_disp);
  
  m_d2p = Attribute_handle (new Distance_to_plane (m_points->point_set()->begin(),
                                 m_points->point_set()->end(),
                                 m_points->point_set()->point_map(),
                                 m_helper->eigen(),
                                 1.));
  m_psc->add_attribute (m_d2p);
  
  if (m_points->point_set()->has_normal_map())
    m_verti = Attribute_handle (new Verticality (m_points->point_set()->begin(),
                               m_points->point_set()->end(),
                               m_points->point_set()->normal_map(),
                               1.));
  else
    m_verti = Attribute_handle (new Verticality (m_points->point_set()->begin(),
                               m_points->point_set()->end(),
                               m_helper->eigen(),
                               1.));

  m_psc->add_attribute (m_verti);

  m_elev = Attribute_handle (new Elevation (m_points->point_set()->begin(),
                          m_points->point_set()->end(),
                          m_points->point_set()->point_map(),
                          _bbox,
                          m_helper->grid(),
                          grid_resolution,
                          radius_neighbors,
                          radius_dtm,
                          1.));
  m_psc->add_attribute (m_elev);

  if (m_points->point_set()->has_colors())
    m_col_att = Attribute_handle (new Color_att (m_points->point_set()->begin(),
                               m_points->point_set()->end(),
                               Point_set_color_map<Point_set>(m_points->point_set()),
                               1.,
                               c.hue(), c.saturation(), c.value()));
  else
    m_col_att = Attribute_handle (new Empty_color ());

  m_psc->add_attribute (m_col_att);

}


void Scene_point_set_classification_item::compute_ransac (const double& radius_neighbors)
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
  
  op.epsilon = radius_neighbors;
  op.cluster_epsilon = radius_neighbors;
  op.normal_threshold = 0.9;

  std::cerr << "Computing RANSAC..." << std::endl;
  shape_detection.detect(op);
  
  CGAL::regularize_planes (shape_detection, true, true, true, false,
                           10., radius_neighbors);
  
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


void Scene_point_set_classification_item::compute_clusters (const double& radius_neighbors)
{
  Q_ASSERT (m_psc != NULL);
  if (m_helper == NULL)
    {
      std::cerr << "Error: no neighborhood" << std::endl;
      return;
    }

  m_psc->cluster_points (m_helper->neighborhood(), radius_neighbors);

  invalidateOpenGLBuffers();
}



void Scene_point_set_classification_item::extract_2d_outline (double radius,
                                                              std::vector<Kernel::Point_3>& outline)
{
    
  // Projection
  double mean_z = 0.;
  std::size_t nb = 0;
  std::vector<Kernel::Point_2> pts;
  
  for (Point_set::const_iterator it = m_points->point_set()->begin();
       it != m_points->point_set()->end(); ++ it)
    if (m_psc->classification_type_of(*it)->id() == "ground")
      {
        mean_z += m_points->point_set()->point(*it).z();
        nb ++;
      }
    else if (m_psc->classification_type_of(*it)->id() != "vegetation")
      pts.push_back (Kernel::Point_2 (m_points->point_set()->point(*it).x(),
                                      m_points->point_set()->point(*it).y()));

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
                                                                std::vector<Kernel::Triangle_3>& /*faces*/)
{
//   // Estimate ground Z
//   std::cerr << "Estimating ground..." << std::endl;
//   std::vector<double> z_ground;
//   for (std::size_t i = 0; i < m_points->point_set().size(); ++ i)
//     if (m_psc->classification_type_of(i)->id() == "ground")
//       z_ground.push_back (m_points->point_set()[i].z());

//   std::sort (z_ground.begin(), z_ground.end());
//   double z_median = z_ground[z_ground.size() / 2];
//   Kernel::Plane_3 ground (Kernel::Point_3 (0., 0., z_median), Kernel::Vector_3 (0., 0., 1.));
  
//   // Project facade planes and get lines
//   std::cerr << "Projecting facade planes..." << std::endl;
//   std::vector<bool> is_plane_facade (m_psc->groups.size(), false);
//   std::vector<std::vector<Kernel::Point_3> > facade_pts (m_psc->groups.size());
  
//   for (std::size_t i = 0; i < m_psc->HPS.size(); ++ i)
//     if (m_psc->HPS[i].AE_label == index_facade && m_psc->HPS[i].group != (std::size_t)(-1))
//       {
//         is_plane_facade[m_psc->HPS[i].group] = true;
//         facade_pts[m_psc->HPS[i].group].push_back (m_psc->HPS[i].position);
//       }


//   std::vector<Kernel::Segment_3> borders;
//   borders.push_back (Kernel::Segment_3 (Kernel::Point_3 (_bbox.xmin(), _bbox.ymin(), z_median),
//                                         Kernel::Point_3 (_bbox.xmin(), _bbox.ymax(), z_median)));
//   borders.push_back (Kernel::Segment_3 (Kernel::Point_3 (_bbox.xmin(), _bbox.ymax(), z_median),
//                                         Kernel::Point_3 (_bbox.xmax(), _bbox.ymax(), z_median)));
//   borders.push_back (Kernel::Segment_3 (Kernel::Point_3 (_bbox.xmax(), _bbox.ymax(), z_median),
//                                         Kernel::Point_3 (_bbox.xmax(), _bbox.ymin(), z_median)));
//   borders.push_back (Kernel::Segment_3 (Kernel::Point_3 (_bbox.xmax(), _bbox.ymin(), z_median),
//                                         Kernel::Point_3 (_bbox.xmin(), _bbox.ymin(), z_median)));

//   std::vector<Kernel::Segment_3> lines;
//   for (std::size_t i = 0; i < is_plane_facade.size(); ++ i)
//     if (is_plane_facade[i])
//       {
//         Kernel::Plane_3 plane = m_psc->groups[i];
//         Kernel::Vector_3 normal = plane.orthogonal_vector();
//         normal = normal / std::sqrt (normal * normal);
//         // double horiz = normal * Kernel::Vector_3 (0., 0., 1.);
//         // if (horiz > 0.17 || horiz < -0.17)
//         //   continue;

//         Kernel::Point_3 c = CGAL::centroid (facade_pts[i].begin(), facade_pts[i].end());
//         normal = Kernel::Vector_3 (normal.x(), normal.y(), 0.);
//         plane = Kernel::Plane_3 (c, normal);
        
//         //#define INFINITE_LINES
// #if defined(INFINITE_LINES)

//         Kernel::Point_3 source, target;
//         bool source_found = false;

//         for (std::size_t j = 0; j < borders.size(); ++ j)
//           {
//             CGAL::Object obj = CGAL::intersection (borders[j], plane);
//             Kernel::Point_3 inter;
//             if (CGAL::assign (inter, obj))
//               {
//                 if (source_found)
//                   {
//                     target = inter;
//                     lines.push_back (Kernel::Segment_3 (source, target));
//                     break;
//                   }
//                 else
//                   {
//                     source = inter;
//                     source_found = true;
//                   }
//               }
//           }
// #else
//         CGAL::Object obj = CGAL::intersection (ground, plane);
//         Kernel::Line_3 line;

//         // if (facade_pts[i].size() < 100)
//         //   continue;
        
//         if (CGAL::assign (line, obj))
//           {
//             Kernel::Point_3 origin = line.projection (CGAL::ORIGIN);
//             Kernel::Vector_3 v (line);
//             double min = std::numeric_limits<double>::max();
//             double max = -std::numeric_limits<double>::max();
//             Kernel::Point_3 source = CGAL::ORIGIN;
//             Kernel::Point_3 target = CGAL::ORIGIN;
//             for (std::size_t j = 0; j < facade_pts[i].size(); ++ j)
//               {
//                 Kernel::Vector_3 dir (origin, facade_pts[i][j]);
//                 double prod = v * dir;
//                 if (prod > max)
//                   {
//                     max = prod;
//                     target = facade_pts[i][j];
//                   }
//                 if (prod < min)
//                   {
//                     min = prod;
//                     source = facade_pts[i][j];
//                   }
//               }
//             source = line.projection (source);
//             target = line.projection (target);
//             Kernel::Point_3 middle (0.5 * (source.x() + target.x()),
//                                     0.5 * (source.y() + target.y()),
//                                     0.5 * (source.z() + target.z()));

//             source = middle + 1.1 * (source - middle);
//             target = middle + 1.1 * (target - middle);
//             lines.push_back (Kernel::Segment_3 (source, target));
//           }
// #endif
//       }

//   std::cerr << " -> Found " << lines.size() << " line(s)" << std::endl;
  
//   // Build 2D constrained Delaunay triangulation
//   std::cerr << "Building 2D constrained Delaunay triangulation..." << std::endl;
    
//   typedef CGAL::Triangulation_vertex_base_2<Kernel>                Vb;
//   typedef CGAL::Triangulation_face_base_with_info_2<int,Kernel>    Fbwi;
//   typedef CGAL::Constrained_triangulation_face_base_2<Kernel,Fbwi> Fb;
//   typedef CGAL::Triangulation_data_structure_2<Vb,Fb>              TDS;
//   typedef CGAL::Exact_predicates_tag                               Itag;
//   typedef CGAL::Constrained_Delaunay_triangulation_2<Kernel, TDS, Itag> CDT;

//   CDT cdt;
//   cdt.insert (Kernel::Point_2 (_bbox.xmin(), _bbox.ymin()));
//   cdt.insert (Kernel::Point_2 (_bbox.xmin(), _bbox.ymax()));
//   cdt.insert (Kernel::Point_2 (_bbox.xmax(), _bbox.ymin()));
//   cdt.insert (Kernel::Point_2 (_bbox.xmax(), _bbox.ymax()));
//   for (std::size_t i = 0; i < lines.size(); ++ i)
//     cdt.insert_constraint (Kernel::Point_2 (lines[i].source().x(), lines[i].source().y()),
//                            Kernel::Point_2 (lines[i].target().x(), lines[i].target().y()));

//   std::cerr << " -> " << cdt.number_of_faces () << " face(s) created" << std::endl;
//   for (CDT::Finite_faces_iterator it = cdt.finite_faces_begin(); it != cdt.finite_faces_end(); ++ it)
//     it->info() = -10;

//   // Project ground + roof points on faces
//   std::cerr << "Projecting ground and roof points on faces..." << std::endl;
//   for (std::size_t i = 0; i < m_psc->HPS.size(); ++ i)
//     {
//       int iter = 0;
//       if (m_psc->HPS[i].AE_label == index_roof)
//         iter = 1;
//       else if (m_psc->HPS[i].AE_label == index_ground)
//         iter = -1;
//       else
//         continue;

//       CDT::Face_handle f = cdt.locate (Kernel::Point_2 (m_psc->HPS[i].position.x(),
//                                                         m_psc->HPS[i].position.y()));

//       f->info() += iter;
//     }

  
//   // Label and extract faces
//   for (CDT::Finite_faces_iterator it = cdt.finite_faces_begin(); it != cdt.finite_faces_end(); ++ it)
//     if (it->info() > 0)
//       faces.push_back (Kernel::Triangle_3 (Kernel::Point_3 (it->vertex (0)->point().x(),
//                                                             it->vertex (0)->point().y(),
//                                                             z_median),
//                                            Kernel::Point_3 (it->vertex (1)->point().x(),
//                                                             it->vertex (1)->point().y(),
//                                                             z_median),
//                                            Kernel::Point_3 (it->vertex (2)->point().x(),
//                                                             it->vertex (2)->point().y(),
//                                                             z_median)));
//   std::cerr << " -> Found " << faces.size() << " building face(s)" << std::endl;
}


void Scene_point_set_classification_item::extract_facades (double /*radius*/,
                                                           std::vector<Kernel::Triangle_3>& /*faces*/)
{
  // std::size_t index_facade = 0;
  // for (std::size_t i = 0; i < m_psc->segmentation_classes.size(); ++ i)
  //   if (m_psc->get_classification_type(i)->id() == "facade")
  //     {
  //       index_facade = i;
  //       break;
  //     }

  // typedef CGAL::Alpha_shape_vertex_base_2<Kernel> Vb;
  // typedef CGAL::Alpha_shape_face_base_2<Kernel>  Fb;
  // typedef CGAL::Triangulation_data_structure_2<Vb,Fb> Tds;
  // typedef CGAL::Delaunay_triangulation_2<Kernel,Tds> Triangulation_2;
  // typedef CGAL::Alpha_shape_2<Triangulation_2>  Alpha_shape_2;
  
  // std::vector<std::vector<Kernel::Point_3> > facade_pts (m_psc->groups.size());
    
  // for (std::size_t i = 0; i < m_psc->HPS.size(); ++ i)
  //   if (m_psc->HPS[i].AE_label == index_facade && m_psc->HPS[i].group != (std::size_t)(-1))
  //     {
  //       facade_pts[m_psc->HPS[i].group].push_back (m_psc->HPS[i].position);
  //     }

  // for (std::size_t i = 0; i < facade_pts.size(); ++ i)
  //   {
  //     Kernel::Plane_3 plane = m_psc->groups[i];
  //     Kernel::Vector_3 normal = plane.orthogonal_vector();
  //     Kernel::Point_3 c = CGAL::centroid (facade_pts[i].begin(), facade_pts[i].end());
  //     normal = Kernel::Vector_3 (normal.x(), normal.y(), 0.);
  //     plane = Kernel::Plane_3 (c, normal);

  //     std::vector<Kernel::Point_2> projections;
  //     projections.reserve (facade_pts[i].size ());

  //     for (std::size_t j = 0; j < facade_pts[i].size (); ++ j)
  //       projections.push_back (plane.to_2d (facade_pts[i][j]));

  //     Alpha_shape_2 ashape (projections.begin (), projections.end (), radius);
  
  //     for (Alpha_shape_2::Finite_faces_iterator it = ashape.finite_faces_begin ();
  //          it != ashape.finite_faces_end (); ++ it)
  //       {
  //         if (ashape.classify (it) != Alpha_shape_2::INTERIOR)
  //           continue;
  //         faces.push_back (Kernel::Triangle_3 (plane.to_3d (it->vertex(0)->point()),
  //                                              plane.to_3d (it->vertex(1)->point()),
  //                                              plane.to_3d (it->vertex(2)->point())));
  //       }
  //   }
}

void Scene_point_set_classification_item::train(std::vector<std::string>& classes,
                                                std::vector<QColor>& colors,
                                                std::size_t nb_trials)
{

  if (m_helper == NULL)
    {
      std::cerr << "Error: features not computed" << std::endl;
      return;
    }

  m_psc->clear_attributes();

  m_helper->generate_attributes (*m_psc,
                                 m_points->point_set()->begin(),
                                 m_points->point_set()->end(),
                                 m_points->point_set()->point_map());
  
  m_psc->clear_classification_types();
  std::vector<std::pair<Type_handle, QColor> > predef;
  for (std::size_t i = 0; i < classes.size(); ++ i)
    {
      predef.push_back (std::make_pair (get_or_add_classification_type(classes[i].c_str(), colors[i]),
                                        colors[i]));
      if (m_psc->number_of_classification_types() != i + 1)
        m_psc->add_classification_type (predef.back().first);

    }
  m_predefined_types.swap (predef);

  m_psc->training(nb_trials);  
  m_psc->run();
  
  invalidateOpenGLBuffers();
  Q_EMIT itemChanged();
}

bool Scene_point_set_classification_item::run_auto (int method, double weight,
                                                    double radius_neighbors)
{
  std::cerr << "Run auto" << std::endl;
  m_points->point_set()->collect_garbage();
  if (method == 0)
    m_psc->run();
  else if (method == 1)
    m_psc->run_with_graphcut (m_helper->neighborhood(), weight);
  else if (method == 2)
    m_psc->run_with_groups (m_helper->neighborhood(), radius_neighbors);

  return true;
}
