#ifndef POINT_SET_CLASSIFICATION_ITEM_H
#define POINT_SET_CLASSIFICATION_ITEM_H

#include <CGAL/Three/Scene_item.h>
#include <CGAL/IO/read_ply_points.h>
#include <CGAL/Point_set_classification.h>
#include <CGAL/Data_classification/Segmentation_attribute_scatter.h>
#include <CGAL/Data_classification/Segmentation_attribute_elevation.h>
#include <CGAL/Data_classification/Segmentation_attribute_horizontality.h>
#include <CGAL/Data_classification/Segmentation_attribute_nonplanarity.h>
#include <CGAL/Data_classification/Segmentation_attribute_color.h>


#include "Scene_point_set_classification_item_config.h"
#include "Scene_points_with_normal_item.h"
#include "Polyhedron_type_fwd.h"
#include "Kernel_type.h"
#include "Point_set_3.h"

#include <iostream>

typedef CGAL::Point_set_classification<Kernel> PSC;

typedef CGAL::Segmentation_attribute_scatter<Kernel> Scatter;
typedef CGAL::Segmentation_attribute_elevation<Kernel> Elevation;
typedef CGAL::Segmentation_attribute_horizontality<Kernel> Horizontality;
typedef CGAL::Segmentation_attribute_nonplanarity<Kernel> NonPlanarity;
typedef CGAL::Segmentation_attribute_color<Kernel> ColorSeg;

typedef CGAL::Data_classification::RGB_Color Color;

class QMenu;
class QAction;

// This class represents a point set in the OpenGL scene
class SCENE_POINT_SET_CLASSIFICATION_ITEM_EXPORT Scene_point_set_classification_item
  : public CGAL::Three::Scene_item
{
  Q_OBJECT

public:
  
  Scene_point_set_classification_item(PSC* psc = NULL);
  Scene_point_set_classification_item(const Scene_points_with_normal_item* points,
                                      double grid_resolution);
  Scene_point_set_classification_item(const Scene_point_set_classification_item& toCopy);
  ~Scene_point_set_classification_item();
  
  Scene_point_set_classification_item* clone() const;

  struct Region
  {
    std::vector<std::size_t> indices;
    double x_min, x_max, y_min, y_max;
    Region () : indices(),
                x_min(std::numeric_limits<double>::max()),
                x_max(-std::numeric_limits<double>::max()),
                y_min(std::numeric_limits<double>::max()),
                y_max(-std::numeric_limits<double>::max()) { }
  };

  template <bool use_y>
  struct Sort_by_coordinate
  {
    std::vector<PSC::HPoint>& m_hps;
    Sort_by_coordinate (std::vector<PSC::HPoint>& hps) : m_hps (hps) { }
    bool operator() (const std::size_t& a, const std::size_t& b)
    {
      if (use_y)
        return m_hps[a].position.y() < m_hps[b].position.y();
      else
        return m_hps[a].position.x() < m_hps[b].position.x();
    }
  };
  
  template <typename OutputIterator>
  bool segment_point_set (std::size_t nb_max_pt, OutputIterator output)
  {
    if (m_psc->HPS.size() < nb_max_pt)
      return false;

    std::list<Region> queue;
    queue.push_front (Region());
    
    for (std::size_t i = 0; i < m_psc->HPS.size(); ++ i)
      {
        queue.front().indices.push_back (i);
        if (m_psc->HPS[i].position.x() < queue.front().x_min)
          queue.front().x_min = m_psc->HPS[i].position.x();
        if (m_psc->HPS[i].position.x() > queue.front().x_max)
          queue.front().x_max = m_psc->HPS[i].position.x();
        if (m_psc->HPS[i].position.y() < queue.front().y_min)
          queue.front().y_min = m_psc->HPS[i].position.y();
        if (m_psc->HPS[i].position.y() > queue.front().y_max)
          queue.front().y_max = m_psc->HPS[i].position.y();
      }

    while (!(queue.empty()))
      {
        Region& current = queue.front();

        if (current.indices.size() < nb_max_pt)
          {
            std::vector<Kernel::Point_3> dummy;

            PSC* psc = new PSC (dummy.begin(), dummy.end(), m_psc->m_grid_resolution);
            for (std::size_t i = 0; i < current.indices.size(); ++ i)
              {
                psc->HPS.push_back (PSC::HPoint());
                psc->HPS.back().position = m_psc->HPS[current.indices[i]].position;
                psc->HPS.back().echo = m_psc->HPS[current.indices[i]].echo;
                psc->HPS.back().ind_x = m_psc->HPS[current.indices[i]].ind_x;
                psc->HPS.back().ind_y = m_psc->HPS[current.indices[i]].ind_y;
                psc->HPS.back().group = m_psc->HPS[current.indices[i]].group;
                psc->HPS.back().AE_label = m_psc->HPS[current.indices[i]].AE_label;
                psc->HPS.back().neighbor = m_psc->HPS[current.indices[i]].neighbor;
                psc->HPS.back().confidence = m_psc->HPS[current.indices[i]].confidence;
                psc->HPS.back().color = m_psc->HPS[current.indices[i]].color;
              }
            *(output ++) = new Scene_point_set_classification_item (psc);
          }
        else
          {
            if (current.x_max - current.x_min > current.y_max - current.y_min)
              {
                std::sort (current.indices.begin(), current.indices.end(),
                           Sort_by_coordinate<false>(m_psc->HPS));
                
                queue.push_back(Region());
                queue.push_back(Region());
                std::list<Region>::iterator it = queue.end();
                Region& positive = *(-- it);
                Region& negative = *(-- it);
                
                negative.y_min = current.y_min;
                positive.y_min = current.y_min;
                negative.y_max = current.y_max;
                positive.y_max = current.y_max;
                
                negative.x_min = current.x_min;
                positive.x_max = current.x_max;
                std::size_t i = 0;
                
                for (; i < current.indices.size() / 2; ++ i)
                  negative.indices.push_back (current.indices[i]);
                double med_x = 0.5 * (m_psc->HPS[current.indices[i]].position.x()
                                      + m_psc->HPS[current.indices[i+1]].position.x());
                negative.x_max = med_x;
                positive.x_min = med_x;
                
                for (; i < current.indices.size(); ++ i)
                  positive.indices.push_back (current.indices[i]);
              }
            else
              {
                std::sort (current.indices.begin(), current.indices.end(),
                           Sort_by_coordinate<true>(m_psc->HPS));

                queue.push_back(Region());
                queue.push_back(Region());
                std::list<Region>::iterator it = queue.end();
                Region& positive = *(-- it);
                Region& negative = *(-- it);

                negative.x_min = current.x_min;
                positive.x_min = current.x_min;
                negative.x_max = current.x_max;
                positive.x_max = current.x_max;

                negative.y_min = current.y_min;
                positive.y_max = current.y_max;
                std::size_t i = 0;
                
                for (; i < current.indices.size() / 2; ++ i)
                  negative.indices.push_back (current.indices[i]);
                double med_y = 0.5 * (m_psc->HPS[current.indices[i]].position.y()
                                      + m_psc->HPS[current.indices[i+1]].position.y());
                negative.y_max = med_y;
                positive.y_min = med_y;
                
                for (; i < current.indices.size(); ++ i)
                  positive.indices.push_back (current.indices[i]);
              }
          }
        
        queue.pop_front ();
      }
    
    return true;
  }
  
  // Function to override the context menu
  QMenu* contextMenu();

  // IO
  bool read_ply_point_set(std::istream& in, double grid_resolution);
  bool write_ply_point_set(std::ostream& out) const;
  // Function for displaying meta-data of the item
  virtual QString toolTip() const;

  virtual void invalidateOpenGLBuffers();

  // Indicate if rendering mode is supported
  virtual bool supportsRenderingMode(RenderingMode m) const;

  virtual void draw_edges(CGAL::Three::Viewer_interface* viewer) const;
  virtual void draw_points(CGAL::Three::Viewer_interface*) const;
  virtual void draw_splats(CGAL::Three::Viewer_interface*) const;
  

  // Gets dimensions
  virtual bool isFinite() const { return true; }
  virtual bool isEmpty() const;
  virtual void compute_bbox() const;

  virtual void setRenderingMode(RenderingMode m);

  void change_color (int index);
  int real_index_color() const;
  
  void estimate_parameters (double& grid_resolution,
                            double& radius_neighbors,
                            double& radius_dtm);
  void compute_features (const double& grid_resolution,
                         const double& radius_neighbors,
                         const double& radius_dtm,
                         const QColor& c);
  void compute_ransac (const double& radius_neighbors);
  void compute_clusters (const double& radius_neighbors);

  void save_2d_image ();
  
  template <typename Classes>
  bool run (double weight_scat, double weight_plan,
            double weight_hori, double weight_elev, double weight_colo,
            bool multiply, Classes& classes, int method, double weight)
  {
  if (m_scat == NULL || m_plan == NULL || m_hori == NULL || m_elev == NULL || m_colo == NULL)
    return false;

  m_psc->m_grid_resolution = weight;

  m_scat->weight = std::tan ((1. - weight_scat) * (CGAL_PI/2));
  m_plan->weight = std::tan ((1. - weight_plan) * (CGAL_PI/2));
  m_hori->weight = std::tan ((1. - weight_hori) * (CGAL_PI/2));
  m_elev->weight = std::tan ((1. - weight_elev) * (CGAL_PI/2));
  m_colo->weight = std::tan ((1. - weight_colo) * (CGAL_PI/2));
  
  // m_scat->weight = 2 * (1. - weight_scat) * m_scat->max;
  // m_plan->weight = 2 * (1. - weight_plan) * m_plan->max;
  // m_hori->weight = 2 * (1. - weight_hori) * m_hori->max;
  // m_elev->weight = 2 * (1. - weight_elev) * m_elev->mean;
  // m_colo->weight = 2 * (1. - weight_colo) * m_colo->max;

  
  if (!(m_psc->segmentation_classes.empty()))
    {
      for (std::size_t i = 0; i < m_psc->segmentation_classes.size(); ++ i)
        delete m_psc->segmentation_classes[i];
      m_psc->segmentation_classes.clear();
    }

  std::cerr << "Running with:" << std::endl;
  for (std::size_t i = 0; i < classes.size(); ++ i)
    {
      if (!(classes[i].checkbox->isChecked()))
        continue;
      
      m_psc->segmentation_classes.push_back (new CGAL::Classification_type
                                             (classes[i].label->text().toLower().toStdString().c_str()));
      m_psc->segmentation_classes.back()->change_attribute_effect (m_scat, (CGAL::Classification_type::Attribute_side)(classes[i].combo[0]->currentIndex()));
      m_psc->segmentation_classes.back()->change_attribute_effect (m_plan, (CGAL::Classification_type::Attribute_side)(classes[i].combo[1]->currentIndex()));
      m_psc->segmentation_classes.back()->change_attribute_effect (m_hori, (CGAL::Classification_type::Attribute_side)(classes[i].combo[2]->currentIndex()));
      m_psc->segmentation_classes.back()->change_attribute_effect (m_elev, (CGAL::Classification_type::Attribute_side)(classes[i].combo[3]->currentIndex()));
      m_psc->segmentation_classes.back()->change_attribute_effect (m_colo, (CGAL::Classification_type::Attribute_side)(classes[i].combo[4]->currentIndex()));

      std::cerr << " * ";
      m_psc->segmentation_classes.back()->info();
    }
  

  std::cerr << "Weights: " << m_scat->weight << " " << m_plan->weight << " " << m_hori->weight
            << " " << m_elev->weight << " " << m_colo->weight << std::endl;

  m_psc->m_multiplicative = multiply;
  
  m_psc->point_cloud_classification(method);

  // if (method != 0)
  //   save_2d_image();
  
  invalidateOpenGLBuffers();
  return false;
}

  template <typename ItemContainer>
  void generate_point_sets (ItemContainer& items)
  {
    for (std::size_t i = 0; i < m_psc->HPS.size(); ++ i)
      {
        int c = m_psc->HPS[i].AE_label;
        if (c == -1)
          continue;
        
        if (m_psc->segmentation_classes[c]->id() == "vegetation")
          items[0]->point_set()->push_back (m_psc->HPS[i].position);
        else if (m_psc->segmentation_classes[c]->id() == "ground")
          items[1]->point_set()->push_back (m_psc->HPS[i].position);
        else if (m_psc->segmentation_classes[c]->id() == "road")
          items[2]->point_set()->push_back (m_psc->HPS[i].position);
        else if (m_psc->segmentation_classes[c]->id() == "roof")
          items[3]->point_set()->push_back (m_psc->HPS[i].position);
        else if (m_psc->segmentation_classes[c]->id() == "facade")
          items[4]->point_set()->push_back (m_psc->HPS[i].position);
        else if (m_psc->segmentation_classes[c]->id() == "building")
          items[5]->point_set()->push_back (m_psc->HPS[i].position);
      }

  }

  void extract_2d_outline (double radius,
                           std::vector<Kernel::Point_3>& outline);
  void extract_building_map (double radius,
                             std::vector<Kernel::Triangle_3>& faces);
  void extract_facades (double radius,
                        std::vector<Kernel::Triangle_3>& faces);

                                   
public Q_SLOTS:

// Data
private:
  PSC* m_psc;
  Scatter* m_scat;
  Elevation* m_elev;
  Horizontality* m_hori;
  NonPlanarity* m_plan;
  ColorSeg* m_colo;
  std::vector<Color> m_real_colors;

  int m_index_color;
  
  enum VAOs {
      Edges=0,
      ThePoints,
      NbOfVaos = ThePoints+1
  };
  enum VBOs {
    Edges_vertices = 0,
    Points_vertices,
    Points_colors,
    NbOfVbos = Points_colors+1
  };

  mutable std::vector<double> positions_lines;
  mutable std::vector<double> positions_points;
  mutable std::vector<double> colors_points;
  mutable std::size_t nb_points;
  mutable std::size_t nb_lines;
   

  mutable QOpenGLShaderProgram *program;

  using CGAL::Three::Scene_item::initialize_buffers;
  void initialize_buffers(CGAL::Three::Viewer_interface *viewer) const;

  void compute_normals_and_vertices() const;


}; // end class Scene_point_set_classification_item


template <typename Floating>
class My_ply_interpreter
{
  std::vector<Kernel::Point_3>& points;
  std::vector<Color>& colors;
  std::vector<unsigned char>& echo;
  bool echo_float;
  
public:
  My_ply_interpreter (std::vector<Kernel::Point_3>& points,
                      std::vector<Color>& colors,
                      std::vector<unsigned char>& echo)
    : points (points), colors (colors), echo (echo), echo_float (false)
  { }

  // Init and test if input file contains the right properties
  bool is_applicable (CGAL::Ply_reader& reader)
  {
    if (reader.does_tag_exist<float> ("scalar_Return_Number"))
      echo_float = true;
    
    return reader.does_tag_exist<Floating> ("x")
      && reader.does_tag_exist<Floating> ("y")
      && reader.does_tag_exist<Floating> ("z");
  }

  // Describes how to process one line (= one point object)
  void process_line (CGAL::Ply_reader& reader)
  {
    Floating x = (Floating)0., y = (Floating)0., z = (Floating)0.;
    Color c = {{ 0, 0, 0 }};
    unsigned char e = 0;

    reader.assign (x, "x");
    reader.assign (y, "y");
    reader.assign (z, "z");
    reader.assign (c[0], "red");
    reader.assign (c[1], "green");
    reader.assign (c[2], "blue");
    if (echo_float)
      {
        float f = 0.;
        reader.assign (f, "scalar_Return_Number");
        e = (unsigned char)f;
      }
    else
      reader.assign (e, "echo");

    points.push_back (Kernel::Point_3 (x, y, z));
    colors.push_back (c);
    echo.push_back (e);
  }

};

template <typename HPS>
struct HPS_property_map{
  const std::vector<HPS>* points;

  typedef Kernel::Point_3 value_type;
  typedef const value_type& reference;
  typedef std::size_t key_type;
  typedef boost::lvalue_property_map_tag category;
  
  HPS_property_map () : points(NULL)
  {
  }
  HPS_property_map (std::vector<HPS>* pts) : points(pts)
  {
  }
  ~HPS_property_map()
  {
  }

  HPS_property_map (const HPS_property_map& pmap) : points(pmap.points)
  {
  }
  
  HPS_property_map& operator= (const HPS_property_map& pmap)
  {
    this->points = pmap.points;
    return *this;
  }
  
  reference operator[](key_type k) const
  {
    return (*points)[k].position;
  }
  
  friend reference get(const HPS_property_map& ppmap,key_type i)
  {return ppmap[i];}
};



#endif // POINT_SET_CLASSIFICATION_ITEM_H
