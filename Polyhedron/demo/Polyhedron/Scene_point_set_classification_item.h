#ifndef POINT_SET_CLASSIFICATION_ITEM_H
#define POINT_SET_CLASSIFICATION_ITEM_H

#include <CGAL/Three/Scene_item.h>
#include <CGAL/Classifier.h>
#include <CGAL/Classification/Attribute.h>
#include <CGAL/Classification/Attribute_vertical_dispersion.h>
#include <CGAL/Classification/Attribute_elevation.h>
#include <CGAL/Classification/Attribute_verticality.h>
#include <CGAL/Classification/Attribute_distance_to_plane.h>
#include <CGAL/Classification/Attribute_color.h>
#include <CGAL/Classification/Attribute_echo_scatter.h>
#include <CGAL/Classification/Attributes_eigen.h>
#include <CGAL/Classification/Helper.h>

#include "Scene_point_set_classification_item_config.h"
#include "Scene_points_with_normal_item.h"
#include "Polyhedron_type_fwd.h"
#include "Kernel_type.h"
#include "Point_set_3.h"

#include <iostream>


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
    typedef CGAL::Classification::RGB_Color Color;
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

 public:
  typedef Kernel::Point_3 Point_3;
  typedef Kernel::Vector_3 Vector_3;
  typedef CGAL::Classification::RGB_Color Color;
  
  typedef Point_set::iterator Iterator;
  typedef Point_set::Point_map Point_map;
  typedef Point_set::Vector_map Vector_map;
  typedef Point_set_color_map<Point_set> Color_map;

  typedef CGAL::Classifier<Iterator, Point_map>                                         PSC;
  typedef CGAL::Classification::Helper<Kernel, Iterator, Point_map>                Helper;
  typedef CGAL::Classification::Type_handle                                                Type_handle;
  typedef CGAL::Classification::Attribute_handle                                           Attribute_handle;
  typedef CGAL::Classification::Attribute_vertical_dispersion<Kernel, Iterator, Point_map> Dispersion;
  typedef CGAL::Classification::Attribute_elevation<Kernel, Iterator, Point_map>           Elevation;
  typedef CGAL::Classification::Attribute_verticality<Kernel, Iterator, Point_map>         Verticality;
  typedef CGAL::Classification::Attribute_distance_to_plane<Kernel, Iterator, Point_map>   Distance_to_plane;


 public:
  
  Scene_point_set_classification_item(PSC* psc = NULL);
  Scene_point_set_classification_item(Scene_points_with_normal_item* points);
  Scene_point_set_classification_item(const Scene_point_set_classification_item& toCopy);
  ~Scene_point_set_classification_item();
  
  Scene_point_set_classification_item* clone() const;
  
  // Function to override the context menu
  QMenu* contextMenu();

  // IO
  bool write_ply_point_set(std::ostream& out);
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

  void reset_indices();
  void compute_features ();
  bool run (int method);

  template <typename Item>
  void generate_point_set_items(std::vector<Item*>& items,
                                const char* name)
  {
    std::map<Type_handle, std::size_t> map_types;
    for (std::size_t i = 0; i < m_types.size(); ++ i)
      {
        items.push_back (new Item);
        items.back()->setName (QString("%1 (%2)").arg(name).arg(m_types[i].first->id().c_str()));
        items.back()->setColor (m_types[i].second);
        map_types[m_types[i].first] = i;
      }

    for (Point_set::const_iterator it = m_points->point_set()->begin();
         it != m_points->point_set()->end(); ++ it)
      {
        Type_handle c = m_psc->classification_type_of(*it);
        if (c != Type_handle())
          items[map_types[c]]->point_set()->insert (m_points->point_set()->point(*it));
      }
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

  void callback()
  {
    invalidateOpenGLBuffers();
    Q_EMIT itemChanged();
  }
  

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

  void validate_selection ()
  {
    if (!(m_psc->classification_prepared()))
      m_psc->prepare_classification();

    for (Point_set::const_iterator it = m_points->point_set()->first_selected();
         it != m_points->point_set()->end(); ++ it)
      {
        Type_handle t = m_psc->classification_type_of(*it);
        m_psc->add_training_index (t, *it);
      }

    m_points->resetSelection();
  }

  double& smoothing() { return m_smoothing; }
  std::size_t& nb_scales() { return m_nb_scales; }
  std::size_t& number_of_trials() { return m_nb_trials; }
  bool features_computed() const { return (m_helper != NULL); }

  std::vector<std::pair<Type_handle, QColor> >& types() { return m_types; }
  std::size_t number_of_attributes() const { return m_psc->number_of_attributes(); }
  Attribute_handle attribute(std::size_t i) { return m_psc->get_attribute(i); }
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
          m_psc->remove_classification_type (m_types[i].first);
          m_types.erase (m_types.begin() + i);
          break;
        }
    invalidateOpenGLBuffers();
    Q_EMIT itemChanged();
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
 void fill_display_combo_box (ComboBox* cb, ComboBox* cb1)
  {
    for (std::size_t i = 0; i < m_psc->number_of_attributes(); ++ i)
      {
        std::size_t scale = m_helper->scale_of_attribute(m_psc->get_attribute(i));
        std::ostringstream oss;
        oss << "Attribute " << m_psc->get_attribute(i)->id() << "_" << scale;
        cb->addItem (oss.str().c_str());
        cb1->addItem (oss.str().c_str());
      }
  }

  void save_config(const char* filename)
  {
    if (m_helper == NULL)
      {
        std::cerr << "Error: features not computed" << std::endl;
        return;
      }

    m_helper->save (filename, *m_psc);
  }

  void load_config(const char* filename)
  {
    if (m_helper != NULL)
      delete m_helper;

    m_psc->clear_classification_types();
    m_psc->clear_attributes();
    reset_indices();
    
    bool normals = m_points->point_set()->has_normal_map();
    bool colors = m_points->point_set()->has_colors();
    Point_set::Property_map<boost::uint8_t> echo_map;
    bool echo;
    boost::tie (echo_map, echo) = m_points->point_set()->template property_map<boost::uint8_t>("echo");

    if (!normals && !colors && !echo)
      m_helper = new Helper (filename, *m_psc,
                             m_points->point_set()->begin(),
                             m_points->point_set()->end(),
                             m_points->point_set()->point_map());
    else if (!normals && !colors && echo)
      m_helper = new Helper (filename, *m_psc,
                             m_points->point_set()->begin(),
                             m_points->point_set()->end(),
                             m_points->point_set()->point_map(),
                             CGAL::Empty_property_map<Iterator, Vector_3>(),
                             CGAL::Empty_property_map<Iterator, Color>(),
                             echo_map);
    else if (!normals && colors && !echo)
      m_helper = new Helper (filename, *m_psc,
                             m_points->point_set()->begin(),
                             m_points->point_set()->end(),
                             m_points->point_set()->point_map(),
                             CGAL::Empty_property_map<Iterator, Vector_3>(),
                             Color_map(m_points->point_set()));
    else if (!normals && colors && echo)
      m_helper = new Helper (filename, *m_psc,
                             m_points->point_set()->begin(),
                             m_points->point_set()->end(),
                             m_points->point_set()->point_map(),
                             CGAL::Empty_property_map<Iterator, Vector_3>(),
                             Color_map(m_points->point_set()),
                             echo_map);
    else if (normals && !colors && !echo)
      m_helper = new Helper (filename, *m_psc,
                             m_points->point_set()->begin(),
                             m_points->point_set()->end(),
                             m_points->point_set()->point_map(),
                             m_points->point_set()->normal_map());
    else if (normals && !colors && echo)
      m_helper = new Helper (filename, *m_psc,
                             m_points->point_set()->begin(),
                             m_points->point_set()->end(),
                             m_points->point_set()->point_map(),
                             m_points->point_set()->normal_map(),
                             CGAL::Empty_property_map<Iterator, Color>(),
                             echo_map);
    else if (normals && colors && !echo)
      m_helper = new Helper (filename, *m_psc,
                             m_points->point_set()->begin(),
                             m_points->point_set()->end(),
                             m_points->point_set()->point_map(),
                             m_points->point_set()->normal_map(),
                             CGAL::Empty_property_map<Iterator, Color>());
    else
      m_helper = new Helper (filename, *m_psc,
                             m_points->point_set()->begin(),
                             m_points->point_set()->end(),
                             m_points->point_set()->point_map(),
                             m_points->point_set()->normal_map(),
                             Color_map(m_points->point_set()),
                             echo_map);
    
    std::vector<std::pair<Type_handle, QColor> > new_types;
    for (std::size_t i = 0; i < m_psc->number_of_classification_types(); ++ i)
      {
        Type_handle t = m_psc->get_classification_type(i);
        QColor color (192 + rand() % 60,
                      192 + rand() % 60,
                      192 + rand() % 60);

        for (std::size_t j = 0; j < m_types.size(); ++ j)
          if (t->id() == m_types[j].first->id())
            {
              color = m_types[j].second;
              break;
            }

        new_types.push_back (std::make_pair (t, color));
      }
    m_types.swap (new_types);
  }

 public Q_SLOTS:

  // Data
 private:
  
  Scene_points_with_normal_item* m_points;
  PSC* m_psc;
  Helper* m_helper;

  std::size_t m_nb_scales;

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




#endif // POINT_SET_CLASSIFICATION_ITEM_H
