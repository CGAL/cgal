#ifndef POINT_SET_CLASSIFICATION_ITEM_H
#define POINT_SET_CLASSIFICATION_ITEM_H

#include <CGAL/Three/Scene_item.h>
#include <CGAL/Point_set_classification.h>
#include <CGAL/IO/read_ply_points.h>

#include "Scene_point_set_classification_item_config.h"
#include "Scene_points_with_normal_item.h"
#include "Polyhedron_type_fwd.h"
#include "Kernel_type.h"
#include "Point_set_3.h"

#include <iostream>

typedef CGAL::Point_set_classification<Kernel> PSC;

typedef CGAL::Scatter_segmentation_attribute<Kernel> Scatter;
typedef CGAL::Elevation_segmentation_attribute<Kernel> Elevation;
typedef CGAL::Horizontality_segmentation_attribute<Kernel> Horizontality;
typedef CGAL::Distance_to_plane_segmentation_attribute<Kernel> Planarity;
typedef CGAL::Color_segmentation_attribute<Kernel> ColorSeg;

typedef CGAL::Data_classification::RGB_Color Color;

class QMenu;
class QAction;

// This class represents a point set in the OpenGL scene
class SCENE_POINT_SET_CLASSIFICATION_ITEM_EXPORT Scene_point_set_classification_item
  : public CGAL::Three::Scene_item
{
  Q_OBJECT

public:
  
  Scene_point_set_classification_item();
  Scene_point_set_classification_item(const Scene_points_with_normal_item* points,
                                      double grid_resolution);
  Scene_point_set_classification_item(const Scene_point_set_classification_item& toCopy);
  ~Scene_point_set_classification_item();
  Scene_point_set_classification_item* clone() const;

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
                         const double& radius_dtm);
  void compute_ransac (const double& radius_neighbors);

  template <typename Classes>
  bool run (double weight_scat, double weight_plan,
            double weight_hori, double weight_elev, double weight_colo,
            bool multiply, Classes& classes, int method)
  {
  if (m_scat == NULL || m_plan == NULL || m_hori == NULL || m_elev == NULL || m_colo == NULL)
    return false;

  m_scat->weight = 2 * (1. - weight_scat) * m_scat->max;
  m_plan->weight = 2 * (1. - weight_plan) * m_plan->max;
  m_hori->weight = 2 * (1. - weight_hori) * m_hori->max;
  m_elev->weight = 2 * (1. - weight_elev) * m_elev->mean;
  m_colo->weight = 2 * (1. - weight_colo) * m_colo->max;
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
  Planarity* m_plan;
  ColorSeg* m_colo;
  std::vector<Color> m_real_colors;

  int m_index_color;
  
  enum VAOs {
      ThePoints = 0,
      NbOfVaos = ThePoints+1
  };
  enum VBOs {
    Points_vertices = 0,
    Points_colors,
    NbOfVbos = Points_colors+1
  };

  mutable std::vector<double> positions_points;
  mutable std::vector<double> colors_points;
  mutable std::size_t nb_points;
 

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
    
public:
  My_ply_interpreter (std::vector<Kernel::Point_3>& points,
                      std::vector<Color>& colors)
    : points (points), colors (colors)
  { }

  // Init and test if input file contains the right properties
  bool is_applicable (CGAL::Ply_reader& reader)
  {
    return reader.does_tag_exist<Floating> ("x")
      && reader.does_tag_exist<Floating> ("y")
      && reader.does_tag_exist<Floating> ("z");
  }

  // Describes how to process one line (= one point object)
  void process_line (CGAL::Ply_reader& reader)
  {
    Floating x = (Floating)0., y = (Floating)0., z = (Floating)0.;
    Color c = {{ 0, 0, 0 }};

    reader.assign (x, "x");
    reader.assign (y, "y");
    reader.assign (z, "z");
    reader.assign (c[0], "red");
    reader.assign (c[1], "green");
    reader.assign (c[2], "blue");

    points.push_back (Kernel::Point_3 (x, y, z));
    colors.push_back (c);
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
