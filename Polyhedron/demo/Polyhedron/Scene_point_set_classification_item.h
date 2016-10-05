#ifndef POINT_SET_CLASSIFICATION_ITEM_H
#define POINT_SET_CLASSIFICATION_ITEM_H

#include <CGAL/Three/Scene_item.h>
#include <CGAL/Point_set_classification.h>
#include <CGAL/Data_classification/Attribute.h>
#include <CGAL/Data_classification/Attribute_vertical_dispersion.h>
#include <CGAL/Data_classification/Attribute_elevation.h>
#include <CGAL/Data_classification/Attribute_verticality.h>
#include <CGAL/Data_classification/Attribute_distance_to_plane.h>
#include <CGAL/Data_classification/Attribute_color.h>
#include <CGAL/Data_classification/Attribute_echo_scatter.h>
#include <CGAL/Data_classification/Attributes_eigen.h>
#include <CGAL/Data_classification/Helper.h>
#include <CGAL/Shape_detection_3/Shape_base.h>

#include "Scene_point_set_classification_item_config.h"
#include "Scene_points_with_normal_item.h"
#include "Polyhedron_type_fwd.h"
#include "Kernel_type.h"
#include "Point_set_3.h"

#include <iostream>


template <typename Kernel, typename Iterator, typename Color_map>
class Attribute_empty_color : public CGAL::Data_classification::Attribute_color<Kernel, Iterator, Color_map>
{
  typedef CGAL::Data_classification::Attribute_color<Kernel, Iterator, Color_map> Base;
public:
  Attribute_empty_color ()
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
  typedef CGAL::Data_classification::Helper<Kernel, Iterator, Point_map>                Helper;
  typedef CGAL::Data_classification::Type_handle                                                Type_handle;
  typedef CGAL::Data_classification::Attribute_handle                                           Attribute_handle;
  typedef CGAL::Data_classification::Attribute_vertical_dispersion<Kernel, Iterator, Point_map> Dispersion;
  typedef CGAL::Data_classification::Attribute_elevation<Kernel, Iterator, Point_map>           Elevation;
  typedef CGAL::Data_classification::Attribute_verticality<Kernel, Iterator, Point_map>         Verticality;
  typedef CGAL::Data_classification::Attribute_distance_to_plane<Kernel, Iterator, Point_map>   Distance_to_plane;
  typedef CGAL::Data_classification::Attribute_color<Kernel, Iterator, Color_map>               Color_att;
  typedef Attribute_empty_color<Kernel, Iterator, Color_map>                                    Empty_color;


 public:
  
  Scene_point_set_classification_item(PSC* psc = NULL);
  Scene_point_set_classification_item(Scene_points_with_normal_item* points);
  Scene_point_set_classification_item(const Scene_point_set_classification_item& toCopy);
  ~Scene_point_set_classification_item();
  
  Scene_point_set_classification_item* clone() const;
  
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
  
  void estimate_parameters ();
  void compute_features ();
  void compute_ransac ();
  void compute_clusters ();
  bool run (int method);

  template <typename Item>
  void generate_point_set_items(std::vector<Item>& )
  {
    // TODO
  }

  Scene_points_with_normal_item* points_item() { return m_points; }

  void reset_training_sets()
  {
    m_psc->prepare_classification();
    m_psc->reset_training_sets();
    invalidateOpenGLBuffers();
    Q_EMIT itemChanged();
  }

  void train();


  Type_handle get_classification_type (const char* name)
  {
    for (std::size_t i = 0; i < m_types.size(); ++ i)
      if (m_types[i].first->id() == name)
        return m_types[i].first;
    return Type_handle();
  }
  
  void add_selection_to_training_set (const char* name)
  {
    Type_handle class_type = get_classification_type (name);
    
    if (!(m_psc->classification_prepared()))
      m_psc->prepare_classification();

    for (Point_set::const_iterator it = m_points->point_set()->first_selected();
         it != m_points->point_set()->end(); ++ it)
      m_psc->set_classification_type_of(*it, class_type);

    m_psc->add_training_set(class_type,
                            m_points->point_set()->first_selected(),
                            m_points->point_set()->end());
    m_points->resetSelection();
  }

  double& grid_resolution() { return m_grid_resolution; }
  double& radius_neighbors() { return m_radius_neighbors; }
  double& radius_dtm() { return m_radius_dtm; }
  std::size_t& number_of_trials() { return m_nb_trials; }
  bool features_computed() const { return (m_helper != NULL); }

  const std::vector<std::pair<Type_handle, QColor> >& types() { return m_types; }
  void add_new_type (const char* name, const QColor& color)
  {
    m_types.push_back (std::make_pair (m_psc->add_classification_type(name),
                                       color));
  }

  void remove_type (const char* name)
  {
    for (std::size_t i = 0; i < m_types.size(); ++ i)
      if (m_types[i].first->id() == name)
        {
          m_types.erase (m_types.begin() + i);
          m_psc->remove_classification_type (m_types[i].first);
          break;
        }
  }

  void change_type_color (const char* name, const QColor& color)
  {
    for (std::size_t i = 0; i < m_types.size(); ++ i)
      if (m_types[i].first->id() == name)
        {
          m_types[i].second = color;
          break;
        }
  }

  template <typename ComboBox>
  void fill_display_combo_box (ComboBox* cb)
  {
    for (std::size_t i = 0; i < m_psc->number_of_attributes(); ++ i)
      {
        std::ostringstream oss;
        oss << "Attribute " << m_psc->get_attribute(i)->id();
        cb->addItem (oss.str().c_str());
      }
  }

 public Q_SLOTS:

  // Data
 private:
  
  Scene_points_with_normal_item* m_points;
  PSC* m_psc;
  Helper* m_helper;

  double m_grid_resolution;
  double m_radius_neighbors;
  double m_radius_dtm;

  std::vector<std::pair<Type_handle, QColor> > m_types;
  std::size_t m_nb_trials;
  double m_smoothing;
  
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
