#ifndef POINT_SET_CLASSIFICATION_ITEM_H
#define POINT_SET_CLASSIFICATION_ITEM_H

//#define CGAL_DO_NOT_USE_BOYKOV_KOLMOGOROV_MAXFLOW_SOFTWARE
#define CGAL_CLASSIFICATION_VERBOSE

#include <CGAL/Three/Scene_item.h>
#include <CGAL/Point_set_classifier.h>
#include <CGAL/Classification/Trainer.h>
#include <CGAL/Classification/Feature_base.h>
#include <CGAL/Classification/Feature/Vertical_dispersion.h>
#include <CGAL/Classification/Feature/Elevation.h>
#include <CGAL/Classification/Feature/Verticality.h>
#include <CGAL/Classification/Feature/Distance_to_plane.h>
#include <CGAL/Classification/Feature/Hsv.h>
#include <CGAL/Classification/Feature/Echo_scatter.h>
#include <CGAL/Classification/Feature/Eigen.h>


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
      Color out = {{ (unsigned char)(255 * pm.ps->red(i)),
                     (unsigned char)(255 * pm.ps->green(i)),
                     (unsigned char)(255 * pm.ps->blue(i)) }};
                     
      return out;
    }
  };

 public:
  typedef Kernel::Point_3 Point_3;
  typedef Kernel::Vector_3 Vector_3;
  typedef CGAL::Classification::RGB_Color Color;
  
  typedef Point_set::Point_map Point_map;
  typedef Point_set::Vector_map Vector_map;
  typedef Point_set_color_map<Point_set> Color_map;

  typedef CGAL::Point_set_classifier<Kernel, Point_set, Point_map>                     PSC;
  typedef CGAL::Classification::Trainer<Point_set, Point_map>        Trainer;
  typedef CGAL::Classification::Label_handle                                               Label_handle;
  typedef CGAL::Classification::Feature_handle                                             Feature_handle;
  typedef CGAL::Classification::Feature::Vertical_dispersion<Kernel, Point_set, Point_map> Dispersion;
  typedef CGAL::Classification::Feature::Elevation<Kernel, Point_set, Point_map>           Elevation;
  typedef CGAL::Classification::Feature::Verticality<Kernel, Point_set, Point_map>         Verticality;
  typedef CGAL::Classification::Feature::Distance_to_plane<Kernel, Point_set, Point_map>   Distance_to_plane;


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
    std::map<Label_handle, std::size_t> map_labels;
    for (std::size_t i = 0; i < m_labels.size(); ++ i)
      {
        items.push_back (new Item);
        items.back()->setName (QString("%1 (%2)").arg(name).arg(m_labels[i].first->name().c_str()));
        items.back()->setColor (m_labels[i].second);
        map_labels[m_labels[i].first] = i;
      }

    for (Point_set::const_iterator it = m_points->point_set()->begin();
         it != m_points->point_set()->end(); ++ it)
      {
        Label_handle c = m_psc->label_of(*it);
        if (c != Label_handle())
          items[map_labels[c]]->point_set()->insert (m_points->point_set()->point(*it));
      }
  }

  Scene_points_with_normal_item* points_item() { return m_points; }

  void reset_training_sets()
  {
    m_trainer->reset_inlier_sets();
    invalidateOpenGLBuffers();
    Q_EMIT itemChanged();
  }

  void train();

  void callback()
  {
    invalidateOpenGLBuffers();
    Q_EMIT itemChanged();
  }
  

  Label_handle get_label (const char* name)
  {
    for (std::size_t i = 0; i < m_labels.size(); ++ i)
      if (m_labels[i].first->name() == name)
        return m_labels[i].first;
    return Label_handle();
  }
  
  void add_selection_to_training_set (const char* name)
  {
    Label_handle label = get_label (name);

    for (Point_set::const_iterator it = m_points->point_set()->first_selected();
         it != m_points->point_set()->end(); ++ it)
      m_psc->set_label_of(*it, label);

    m_trainer->set_inliers(label,
                           boost::make_iterator_range(m_points->point_set()->first_selected(),
                                                      m_points->point_set()->end()));

    m_points->resetSelection();
  }

  void validate_selection ()
  {
    for (Point_set::const_iterator it = m_points->point_set()->first_selected();
         it != m_points->point_set()->end(); ++ it)
      {
        Label_handle t = m_psc->label_of(*it);
        m_trainer->set_inlier (t, *it);
      }

    m_points->resetSelection();
  }

  double& smoothing() { return m_smoothing; }
  std::size_t& nb_scales() { return m_nb_scales; }
  std::size_t& number_of_trials() { return m_nb_trials; }
  bool features_computed() const { return (m_psc->number_of_features() != 0); }

  std::vector<std::pair<Label_handle, QColor> >& labels() { return m_labels; }
  std::size_t number_of_features() const { return m_psc->number_of_features(); }
  Feature_handle feature(std::size_t i) { return m_psc->feature(i); }
  void add_new_label (const char* name, const QColor& color)
  {
    m_labels.push_back (std::make_pair (m_psc->add_label(name),
                                       color));
  }

  void remove_label (const char* name)
  {
    for (std::size_t i = 0; i < m_labels.size(); ++ i)
      if (m_labels[i].first->name() == name)
        {
          m_psc->remove_label (m_labels[i].first);
          m_labels.erase (m_labels.begin() + i);
          break;
        }
    invalidateOpenGLBuffers();
    Q_EMIT itemChanged();
  }

  void change_label_color (const char* name, const QColor& color)
  {
    for (std::size_t i = 0; i < m_labels.size(); ++ i)
      if (m_labels[i].first->name() == name)
        {
          m_labels[i].second = color;
          break;
        }
  }

  template <typename ComboBox>
 void fill_display_combo_box (ComboBox* cb, ComboBox* cb1)
  {
    for (std::size_t i = 0; i < m_psc->number_of_features(); ++ i)
      {
        std::size_t scale = m_psc->scale_of_feature(m_psc->feature(i));
        std::ostringstream oss;
        oss << "Feature " << m_psc->feature(i)->name() << "_" << scale;
        cb->addItem (oss.str().c_str());
        cb1->addItem (oss.str().c_str());
      }
  }

  void save_config(const char* filename)
  {
    if (m_psc->number_of_features() == 0)
      {
        std::cerr << "Error: features not computed" << std::endl;
        return;
      }

    std::ofstream f (filename);
    m_psc->save_configuration (f);
  }

  void load_config(const char* filename)
  {
    if (m_psc->number_of_features() != 0)
      m_psc->clear();

    reset_indices();
    
    bool normals = m_points->point_set()->has_normal_map();
    bool colors = m_points->point_set()->has_colors();
    Point_set::Property_map<boost::uint8_t> echo_map;
    bool echo;
    boost::tie (echo_map, echo) = m_points->point_set()->template property_map<boost::uint8_t>("echo");

    std::ifstream f (filename);
    if (!normals && !colors && !echo)
      m_psc->load_configuration (f);
    else if (!normals && !colors && echo)
      m_psc->load_configuration (f, CGAL::Default(), CGAL::Default(), echo_map);
    else if (!normals && colors && !echo)
      m_psc->load_configuration (f, CGAL::Default(), Color_map(m_points->point_set()));
    else if (!normals && colors && echo)
      m_psc->load_configuration (f, CGAL::Default(), Color_map(m_points->point_set()), echo_map);
    else if (normals && !colors && !echo)
      m_psc->load_configuration (f, m_points->point_set()->normal_map());
    else if (normals && !colors && echo)
      m_psc->load_configuration (f, m_points->point_set()->normal_map(), CGAL::Default(), echo_map);
    else if (normals && colors && !echo)
      m_psc->load_configuration (f, m_points->point_set()->normal_map(), Color_map(m_points->point_set()));
    else
      m_psc->load_configuration (f, m_points->point_set()->normal_map(), Color_map(m_points->point_set()), echo_map);
    
    std::vector<std::pair<Label_handle, QColor> > new_labels;
    for (std::size_t i = 0; i < m_psc->number_of_labels(); ++ i)
      {
        Label_handle t = m_psc->label(i);
        QColor color (192 + rand() % 60,
                      192 + rand() % 60,
                      192 + rand() % 60);

        for (std::size_t j = 0; j < m_labels.size(); ++ j)
          if (t->name() == m_labels[j].first->name())
            {
              color = m_labels[j].second;
              break;
            }

        new_labels.push_back (std::make_pair (t, color));
      }
    m_labels.swap (new_labels);
  }

 public Q_SLOTS:

  // Data
 private:
  
  Scene_points_with_normal_item* m_points;
  PSC* m_psc;
  Trainer* m_trainer;

  std::size_t m_nb_scales;

  std::vector<std::pair<Label_handle, QColor> > m_labels;
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
