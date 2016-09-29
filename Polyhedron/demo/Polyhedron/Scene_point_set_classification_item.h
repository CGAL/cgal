#ifndef POINT_SET_CLASSIFICATION_ITEM_H
#define POINT_SET_CLASSIFICATION_ITEM_H

#include <CGAL/Three/Scene_item.h>
#include <CGAL/Point_set_classification.h>
#include <CGAL/Data_classification/Segmentation_attribute_vertical_dispersion.h>
#include <CGAL/Data_classification/Segmentation_attribute_elevation.h>
#include <CGAL/Data_classification/Segmentation_attribute_verticality.h>
#include <CGAL/Data_classification/Segmentation_attribute_distance_to_plane.h>
#include <CGAL/Data_classification/Segmentation_attribute_color.h>
#include <CGAL/Shape_detection_3/Shape_base.h>

#include "Scene_point_set_classification_item_config.h"
#include "Scene_points_with_normal_item.h"
#include "Polyhedron_type_fwd.h"
#include "Kernel_type.h"
#include "Point_set_3.h"

#include <iostream>


template <typename Kernel, typename Iterator, typename Color_map>
class Segmentation_attribute_empty_color : public CGAL::Segmentation_attribute_color<Kernel, Iterator, Color_map>
{
  typedef CGAL::Segmentation_attribute_color<Kernel, Iterator, Color_map> Base;
public:
  Segmentation_attribute_empty_color ()
    : Base (Iterator(), Iterator(),
            Color_map(), 1.)
  { }
  virtual double value (std::size_t) { return 1.; }
};


class QMenu;
class QAction;

// This class represents a point set in the OpenGL scene
class SCENE_POINT_SET_CLASSIFICATION_ITEM_EXPORT Scene_point_set_classification_item
  : public CGAL::Three::Scene_item
{
  Q_OBJECT

  
  template <typename PointSet>
    class Point_set_color_map
  {
  public:
    typedef CGAL::Data_classification::RGB_Color Color;
    typedef typename PointSet::Index key_type;
    typedef Color value_type;
    typedef value_type& reference;
    typedef boost::lvalue_property_map_tag category;
    
    PointSet* ps;

    Point_set_color_map(PointSet* ps = NULL)
      : ps(ps) { }

    friend value_type get (const Point_set_color_map& pm, const key_type& i)
    {
      Color out = {{ pm.ps->red(i),
                     pm.ps->green(i),
                     pm.ps->blue(i) }};
                     
      return out;
    }
  };


  typedef typename Kernel::Point_3 Point_3;
  typedef typename Kernel::Vector_3 Vector_3;
  typedef CGAL::Data_classification::RGB_Color Color;
  
  typedef Point_set::iterator Iterator;
  typedef Point_set::Point_map Point_map;
  typedef Point_set::Vector_map Vector_map;
  typedef Point_set_color_map<Point_set> Color_map;

  typedef CGAL::Point_set_classification<Kernel, Iterator, Point_map>                   PSC;
  typedef CGAL::Data_classification::Planimetric_grid<Kernel, Iterator, Point_map>      Planimetric_grid;
  typedef CGAL::Data_classification::Neighborhood<Kernel, Iterator, Point_map>          Neighborhood;
  typedef CGAL::Data_classification::Local_eigen_analysis<Kernel, Iterator, Point_map>  Local_eigen_analysis;
  typedef CGAL::Segmentation_attribute_vertical_dispersion<Kernel, Iterator, Point_map> Dispersion;
  typedef CGAL::Segmentation_attribute_elevation<Kernel, Iterator, Point_map>           Elevation;
  typedef CGAL::Segmentation_attribute_verticality<Kernel, Iterator, Point_map>         Verticality;
  typedef CGAL::Segmentation_attribute_distance_to_plane<Kernel, Iterator, Point_map>   Distance_to_plane;
  typedef CGAL::Segmentation_attribute_color<Kernel, Iterator, Color_map>               Color_att;
  typedef Segmentation_attribute_empty_color<Kernel, Iterator, Color_map>               Empty_color;


public:
  
  Scene_point_set_classification_item(PSC* psc = NULL);
  Scene_point_set_classification_item(Scene_points_with_normal_item* points);
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
    std::vector<Point_3>& m_hps;
    Sort_by_coordinate (std::vector<Point_3>& hps) : m_hps (hps) { }
    bool operator() (const std::size_t& a, const std::size_t& b)
    {
      if (use_y)
        return m_hps[a].y() < m_hps[b].y();
      else
        return m_hps[a].x() < m_hps[b].x();
    }
  };
  
  template <typename OutputIterator>
  bool segment_point_set (std::size_t nb_max_pt, OutputIterator output)
  {
    return false;
    // if (m_points->point_set()->size() < nb_max_pt)
    //   return false;

    // std::list<Region> queue;
    // queue.push_front (Region());
    
    // for (Points_set::const_iterator it = m_points->point_set()->begin();
    //      it != m_points->point_set()->end(); ++ it)
    //   {
    //     queue.front().indices.push_back (*it);
    //     if (m_points->point_set()[i].x() < queue.front().x_min)
    //       queue.front().x_min = m_points->point_set()->point(it).x();
    //     if (m_points->point_set()->point(it).x() > queue.front().x_max)
    //       queue.front().x_max = m_points->point_set()->point(it).x();
    //     if (m_points->point_set()->point(it).y() < queue.front().y_min)
    //       queue.front().y_min = m_points->point_set()->point(it).y();
    //     if (m_points->point_set()->point(it).y() > queue.front().y_max)
    //       queue.front().y_max = m_points->point_set()->point(it).y();
    //   }

    // while (!(queue.empty()))
    //   {
    //     Region& current = queue.front();

    //     if (current.indices.size() < nb_max_pt)
    //       {
    //         std::vector<Kernel::Point_3> dummy;

    //         Scene_point_set_classification_item* new_item
    //           = new Scene_point_set_classification_item ();
            
    //         for (std::size_t i = 0; i < current.indices.size(); ++ i)
    //           {
    //             new_item->m_points->point_set().push_back(m_points->point_set()[current.indices[i]]);
    //             new_item->m_normals.push_back(m_normals[current.indices[i]]);
    //             new_item->m_colors.push_back(m_colors[current.indices[i]]);
    //           }

    //         *(output ++) = new_item;
    //       }
    //     else
    //       {
    //         if (current.x_max - current.x_min > current.y_max - current.y_min)
    //           {
    //             std::sort (current.indices.begin(), current.indices.end(),
    //                        Sort_by_coordinate<false>(m_points->point_set()));
                
    //             queue.push_back(Region());
    //             queue.push_back(Region());
    //             std::list<Region>::iterator it = queue.end();
    //             Region& positive = *(-- it);
    //             Region& negative = *(-- it);
                
    //             negative.y_min = current.y_min;
    //             positive.y_min = current.y_min;
    //             negative.y_max = current.y_max;
    //             positive.y_max = current.y_max;
                
    //             negative.x_min = current.x_min;
    //             positive.x_max = current.x_max;
    //             std::size_t i = 0;
                
    //             for (; i < current.indices.size() / 2; ++ i)
    //               negative.indices.push_back (current.indices[i]);
    //             double med_x = 0.5 * (m_points->point_set()[current.indices[i]].x()
    //                                   + m_points->point_set()[current.indices[i+1]].x());
    //             negative.x_max = med_x;
    //             positive.x_min = med_x;
                
    //             for (; i < current.indices.size(); ++ i)
    //               positive.indices.push_back (current.indices[i]);
    //           }
    //         else
    //           {
    //             std::sort (current.indices.begin(), current.indices.end(),
    //                        Sort_by_coordinate<true>(m_points->point_set()));

    //             queue.push_back(Region());
    //             queue.push_back(Region());
    //             std::list<Region>::iterator it = queue.end();
    //             Region& positive = *(-- it);
    //             Region& negative = *(-- it);

    //             negative.x_min = current.x_min;
    //             positive.x_min = current.x_min;
    //             negative.x_max = current.x_max;
    //             positive.x_max = current.x_max;

    //             negative.y_min = current.y_min;
    //             positive.y_max = current.y_max;
    //             std::size_t i = 0;
                
    //             for (; i < current.indices.size() / 2; ++ i)
    //               negative.indices.push_back (current.indices[i]);
    //             double med_y = 0.5 * (m_points->point_set()[current.indices[i]].y()
    //                                   + m_points->point_set()[current.indices[i+1]].y());
    //             negative.y_max = med_y;
    //             positive.y_min = med_y;
                
    //             for (; i < current.indices.size(); ++ i)
    //               positive.indices.push_back (current.indices[i]);
    //           }
    //       }
        
    //     queue.pop_front ();
    //   }
    
    // return true;
  }
  
  // Function to override the context menu
  QMenu* contextMenu();

  // IO
  bool write_ply_point_set(std::ostream& out) const;
  // Function for displaying meta-data of the item
  virtual QString toolTip() const;

  virtual void invalidateOpenGLBuffers();

  // Indicate if rendering mode is supported
  virtual bool supportsRenderingMode(RenderingMode m) const;

  virtual void drawEdges(CGAL::Three::Viewer_interface* viewer) const;
  virtual void drawPoints(CGAL::Three::Viewer_interface*) const;
  virtual void drawSplats(CGAL::Three::Viewer_interface*) const;
  

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

  template <typename Classes>
  bool run (double weight_scat, double weight_plan,
            double weight_hori, double weight_elev, double weight_colo,
            Classes& classes, int method, double weight,
            double radius_neighbors)
  {
  if (m_disp == NULL || m_d2p == NULL || m_verti == NULL || m_elev == NULL || m_col_att == NULL)
    return false;

  m_disp->weight = std::tan ((1. - weight_scat) * (CGAL_PI/2));
  m_d2p->weight = std::tan ((1. - weight_plan) * (CGAL_PI/2));
  m_verti->weight = std::tan ((1. - weight_hori) * (CGAL_PI/2));
  m_elev->weight = std::tan ((1. - weight_elev) * (CGAL_PI/2));
  m_col_att->weight = std::tan ((1. - weight_colo) * (CGAL_PI/2));

  m_psc->clear_and_delete_classification_types();

  std::cerr << "Running with:" << std::endl;
  for (std::size_t i = 0; i < classes.size(); ++ i)
    {
      if (!(classes[i].checkbox->isChecked()))
        continue;

      CGAL::Classification_type* ct = new CGAL::Classification_type
        (classes[i].label->text().toLower().toStdString().c_str());
      ct->set_attribute_effect
        (m_disp, (CGAL::Classification_type::Attribute_effect)(classes[i].combo[0]->currentIndex()));
      ct->set_attribute_effect
        (m_d2p, (CGAL::Classification_type::Attribute_effect)(classes[i].combo[1]->currentIndex()));
      ct->set_attribute_effect
        (m_verti, (CGAL::Classification_type::Attribute_effect)(classes[i].combo[2]->currentIndex()));
      ct->set_attribute_effect
        (m_elev, (CGAL::Classification_type::Attribute_effect)(classes[i].combo[3]->currentIndex()));
      ct->set_attribute_effect
        (m_col_att, (CGAL::Classification_type::Attribute_effect)(classes[i].combo[4]->currentIndex()));

      std::cerr << " * ";
      ct->info();

      m_psc->add_classification_type (ct);
    }
  

  std::cerr << "Weights: " << m_disp->weight << " " << m_d2p->weight << " " << m_verti->weight
            << " " << m_elev->weight << " " << m_col_att->weight << std::endl;
  
  if (method == 0)
    m_psc->run();
  else if (method == 1)
    m_psc->run_with_graphcut (*m_neighborhood, weight);
  else if (method == 2)
    m_psc->run_with_groups (*m_neighborhood, radius_neighbors);
  
  invalidateOpenGLBuffers();
  return false;
}

  template <typename ItemContainer>
  void generate_point_sets (ItemContainer& items)
  {
    for (Point_set::const_iterator it = m_points->point_set()->begin();
         it != m_points->point_set()->end(); ++ it)
      {
        CGAL::Classification_type* c = m_psc->classification_type_of (*it);
        if (c == NULL)
          continue;
        
        if (c->id() == "vegetation")
          items[0]->point_set()->insert (m_points->point_set()->point(*it));
        else if (c->id() == "ground")
          items[1]->point_set()->insert (m_points->point_set()->point(*it));
        else if (c->id() == "road")
          items[2]->point_set()->insert (m_points->point_set()->point(*it));
        else if (c->id() == "roof")
          items[3]->point_set()->insert (m_points->point_set()->point(*it));
        else if (c->id() == "facade")
          items[4]->point_set()->insert (m_points->point_set()->point(*it));
        else if (c->id() == "building")
          items[5]->point_set()->insert (m_points->point_set()->point(*it));
      }
  }

  void extract_2d_outline (double radius,
                           std::vector<Kernel::Point_3>& outline);
  void extract_building_map (double radius,
                             std::vector<Kernel::Triangle_3>& faces);
  void extract_facades (double radius,
                        std::vector<Kernel::Triangle_3>& faces);

  Scene_points_with_normal_item* points_item() { return m_points; }

  CGAL::Classification_type* get_or_add_classification_type
    (const char* name, const QColor& color)
  {
    for (std::size_t i = 0; i < m_predefined_types.size(); ++ i)
      if (m_predefined_types[i].first->id() == name)
        {
          m_predefined_types[i].second = color;
          return m_predefined_types[i].first;
        }
    m_predefined_types.push_back
      (std::make_pair (new CGAL::Classification_type (name), color));
    m_psc->add_classification_type (m_predefined_types.back().first);
    return m_predefined_types.back().first;
  }

  void reset_training_sets()
  {
    for (std::size_t i = 0; i < m_predefined_types.size(); ++ i)
      m_predefined_types[i].first->training_set().clear();

    for (Point_set::const_iterator it = m_points->point_set()->begin();
         it != m_points->point_set()->end(); ++ it)
      m_psc->set_classification_type_of(*it, NULL);
    invalidateOpenGLBuffers();
    Q_EMIT itemChanged();
  }

  void train(std::vector<std::string>& classes,
             std::vector<QColor>& colors,
             std::size_t nb_trials);
  
  void add_selection_to_training_set (const char* name, const QColor& color)
  {
    CGAL::Classification_type* class_type = get_or_add_classification_type (name, color);
    if (!(m_psc->classification_prepared()))
      m_psc->prepare_classification();
    
    for (Point_set::const_iterator it = m_points->point_set()->first_selected();
         it != m_points->point_set()->end(); ++ it)
      {
        class_type->training_set().push_back (*it);
        m_psc->set_classification_type_of(*it, class_type);
      }
    m_points->resetSelection();
  }
                                   
public Q_SLOTS:

// Data
private:
  
  Scene_points_with_normal_item* m_points;

  PSC* m_psc;
  std::vector<std::pair<CGAL::Classification_type*, QColor> > m_predefined_types;

  Planimetric_grid* m_grid;
  Neighborhood* m_neighborhood;
  Local_eigen_analysis* m_eigen;
  
  Dispersion* m_disp;
  Elevation* m_elev;
  Verticality* m_verti;
  Distance_to_plane* m_d2p;
  Color_att* m_col_att;

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

  using CGAL::Three::Scene_item::initializeBuffers;
  void initializeBuffers(CGAL::Three::Viewer_interface *viewer) const;

  void compute_normals_and_vertices() const;


}; // end class Scene_point_set_classification_item


template <class Traits>
class My_plane : public CGAL::Shape_detection_3::Shape_base<Traits> {
public:
  /// \cond SKIP_IN_MANUAL
  typedef typename Traits::Point_map Point_map;
  ///< property map to access the location of an input point.
  typedef typename Traits::Normal_map Normal_map;
  ///< property map to access the unoriented normal of an input point.
  typedef typename Traits::FT FT; ///< number type.
  typedef typename Traits::Point_3 Point_3; ///< point type.
  typedef typename Traits::Vector_3 Vector_3;
  /// \endcond

  typedef typename Traits::Plane_3 Plane_3;///< plane type for conversion operator.

public:
  My_plane() : CGAL::Shape_detection_3::Shape_base<Traits>() {}

  /*!
    Conversion operator to `Plane_3` type.
  */
  operator Plane_3() const {
    return Plane_3(this->get_x(m_normal), this->get_y(m_normal), this->get_z(m_normal), m_d);
  }
            
  /*!
    Normal vector of the plane.
  */
  Vector_3 plane_normal() const {
    return m_normal;
  }
            
  /*!
    Signed distance from the origin.
  */
  FT d() const {
    return m_d;
  }
    
  /// \cond SKIP_IN_MANUAL
  /*!
    Computes squared Euclidean distance from query point to the shape.
  */
  FT squared_distance(const Point_3 &p) const {
    FT d = this->scalar_pdct(
                             this->constr_vec(p, m_point_on_primitive), m_normal);
    return d * d;
  }

    
  /*!
    Helper function to write the plane equation and
    number of assigned points into a string.
  */
  std::string info() const {
    std::stringstream sstr;
    sstr << "Type: plane (" << this->get_x(m_normal) << ", " << this->get_y(m_normal) 
         << ", " << this->get_z(m_normal) << ")x - " << m_d << "= 0"
         << " #Pts: " << this->m_indices.size();

    return sstr.str();
  }
  /// \endcond

protected:
  /// \cond SKIP_IN_MANUAL
  virtual void create_shape(const std::vector<std::size_t> &indices) {
    Point_3 p1 = this->point(indices[0]);
    Point_3 p2 = this->point(indices[1]);
    Point_3 p3 = this->point(indices[2]);

    m_normal = this->cross_pdct(
                                this->constr_vec(p2, p1), this->constr_vec(p3, p1));

    if (std::fabs (m_normal * Vector_3 (0., 0., 1.)) > 0.5)
      m_normal = Vector_3 (0., 0., 1.);
    else
      m_normal = Vector_3 (m_normal.x(), m_normal.y(), 0.);
    
    FT length = CGAL::sqrt(this->sqlen(m_normal));

    // Are the points almost singular?
    if (length < (FT)0.0001) {
      return;
    }

    m_normal = this->scale(m_normal, (FT)1.0 / length);
    m_d = -(this->get_x(p1) * this->get_x(m_normal) 
            + this->get_y(p1) * this->get_y(m_normal) 
            + this->get_z(p1) * this->get_z(m_normal));

    //check deviation of the 3 normal
    Vector_3 l_v = this->constr_vec(); 
    for (std::size_t i = 0;i<3;i++) {
      l_v = this->normal(indices[i]);

      if (CGAL::abs(this->scalar_pdct(l_v, m_normal))
          < this->m_normal_threshold * CGAL::sqrt(this->sqlen(l_v))) {
        this->m_is_valid = false;
        return;
      }

      m_point_on_primitive = p1;
      m_base1 = this->cross_pdct(this->constr_vec(p2, p1), m_normal);
      m_base1 = this->scale(m_base1, ((FT)1.0 / CGAL::sqrt(this->sqlen(m_base1))));

      m_base2 = this->cross_pdct(m_base1, m_normal);
      m_base2 = this->scale(m_base2, ((FT)1.0 / CGAL::sqrt(this->sqlen(m_base2))));
    }

    this->m_is_valid = true;
  }

  virtual void parameters(const std::vector<std::size_t> &indices,
                          std::vector<std::pair<FT, FT> > &parameterSpace,
                          FT &,                    
                          FT min[2],
                          FT max[2]) const {
    // Transform first point before to initialize min/max
    Vector_3 p = this->constr_vec(
                                  m_point_on_primitive, this->point(indices[0]));
    FT u = this->scalar_pdct(p, m_base1);
    FT v = this->scalar_pdct(p, m_base2);
    parameterSpace[0] = std::pair<FT, FT>(u, v);
    min[0] = max[0] = u;
    min[1] = max[1] = v;

    for (std::size_t i = 1;i<indices.size();i++) {
      p = this->constr_vec(m_point_on_primitive, this->point(indices[i]));
      u = this->scalar_pdct(p, m_base1);
      v = this->scalar_pdct(p, m_base2);
      min[0] = (std::min<FT>)(min[0], u);
      max[0] = (std::max<FT>)(max[0], u);
      min[1] = (std::min<FT>)(min[1], v);
      max[1] = (std::max<FT>)(max[1], v);
      parameterSpace[i] = std::pair<FT, FT>(u, v);
    }
  }
    
  virtual void squared_distance(const std::vector<std::size_t> &indices,
                                std::vector<FT> &dists) const {
    for (std::size_t i = 0;i<indices.size();i++) {
      const FT d = this->scalar_pdct(
                                     this->constr_vec(m_point_on_primitive, this->point(indices[i])), 
                                     m_normal);
      dists[i] = d * d;
    }
  }

  virtual void cos_to_normal(const std::vector<std::size_t> &indices, 
                             std::vector<FT> &angles) const {
    for (std::size_t i = 0;i<indices.size();i++) {
      angles[i] = CGAL::abs(
                            this->scalar_pdct(this->normal(indices[i]), m_normal));
    }
  }

  FT cos_to_normal(const Point_3 &, const Vector_3 &n) const{
    return CGAL::abs(this->scalar_pdct(n, m_normal));
  } 
    
  virtual std::size_t minimum_sample_size() const {
    return 3;
  }

  virtual bool supports_connected_component() const {
    return true;
  }

private:
  Point_3 m_point_on_primitive;
  Vector_3 m_base1, m_base2, m_normal;
  FT m_d;
  /// \endcond
};



#endif // POINT_SET_CLASSIFICATION_ITEM_H
