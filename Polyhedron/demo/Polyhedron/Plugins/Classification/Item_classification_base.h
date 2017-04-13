#ifndef ITEM_CLASSIFICATION_BASE_H
#define ITEM_CLASSIFICATION_BASE_H

#include <CGAL/Three/Scene_item.h>

#include <QComboBox>

#include <CGAL/Classification/Feature_set.h>
#include <CGAL/Classification/Label_set.h>
#include <CGAL/Classification/Sum_of_weighted_features_predicate.h>

#ifdef CGAL_LINKED_WITH_OPENCV
#include <CGAL/Classification/Random_forest_predicate.h>
#endif

class Item_classification_base
{
public:
  typedef CGAL::Classification::Label_handle   Label_handle;
  typedef CGAL::Classification::Feature_handle Feature_handle;
  typedef CGAL::Classification::Label_set   Label_set;
  typedef CGAL::Classification::Feature_set Feature_set;
  typedef CGAL::Classification::Sum_of_weighted_features_predicate Sum_of_weighted_features;

#ifdef CGAL_LINKED_WITH_OPENCV
  typedef CGAL::Classification::Random_forest_predicate Random_forest;
#endif
  
public:
  
  Item_classification_base() { }
  virtual ~Item_classification_base() { }

  virtual CGAL::Three::Scene_item* item() = 0;
  virtual void erase_item() = 0;

  virtual void compute_features () = 0;

  virtual void add_selection_to_training_set (const char* name) = 0;
  virtual void reset_training_sets() = 0;

  virtual void validate_selection () = 0;
  virtual void train(int predicate) = 0;
  virtual bool run (int method, int predicate) = 0;
  
  virtual void change_color (int index) = 0;
  virtual void generate_one_item_per_label(std::vector<CGAL::Three::Scene_item*>& items,
                                           const char* name) const = 0;

  virtual bool write_output(std::ostream& out) = 0;

  bool features_computed() const { return (m_features.size() != 0); }
  std::size_t number_of_features() const { return m_features.size(); }  
  Feature_handle feature(std::size_t i) { return m_features[i]; }
  float weight (Feature_handle f) const { return m_sowf->weight(f); }
  void set_weight (Feature_handle f, float w) const { m_sowf->set_weight(f,w); }
  Sum_of_weighted_features::Effect effect (Label_handle l, Feature_handle f) const { return m_sowf->effect(l,f); }
  void set_effect (Label_handle l, Feature_handle f, Sum_of_weighted_features::Effect e)
  { m_sowf->set_effect (l, f, e); }
  
  void add_new_label (const char* name, const QColor& color)
  {
    m_labels.add(name);
    m_label_colors.push_back (color);
    delete m_sowf;
    m_sowf = new Sum_of_weighted_features (m_labels, m_features);
    
#ifdef CGAL_LINKED_WITH_OPENCV
    delete m_random_forest;
    m_random_forest = new Random_forest (m_labels, m_features);
#endif
  }
  void remove_label (const char* name)
  {
    for (std::size_t i = 0; i < m_labels.size(); ++ i)
      if (m_labels[i]->name() == name)
        {
          m_labels.remove(m_labels[i]);
          m_label_colors.erase (m_label_colors.begin() + i);
          break;
        }
    delete m_sowf;
    m_sowf = new Sum_of_weighted_features (m_labels, m_features);

#ifdef CGAL_LINKED_WITH_OPENCV
    delete m_random_forest;
    m_random_forest = new Random_forest (m_labels, m_features);
#endif
  }
  std::size_t number_of_labels() const { return m_labels.size(); }
  Label_handle label(std::size_t i) { return m_labels[i]; }

  void fill_display_combo_box (QComboBox* cb, QComboBox* cb1) const
  {
    for (std::size_t i = 0; i < m_features.size(); ++ i)
      {
        std::ostringstream oss;
        oss << "Feature " << m_features[i]->name();
        cb->addItem (oss.str().c_str());
        cb1->addItem (oss.str().c_str());
      }
  }

  void save_config(const char* filename, int predicate)
  {
    if (m_features.size() == 0)
      {
        std::cerr << "Error: features not computed" << std::endl;
        return;
      }

    if (predicate == 0)
    {
      std::ofstream f (filename);
      m_sowf->save_configuration (f);
    }
    else
    {
#ifdef CGAL_LINKED_WITH_OPENCV
      m_random_forest->save_configuration (filename);
#endif
    }
  }
  void load_config(const char* filename, int predicate)
  {
    if (m_features.size() == 0)
      {
        std::cerr << "Error: features not computed" << std::endl;
        return;
      }

    if (predicate == 0)
    {
      std::ifstream f (filename);
      m_sowf->load_configuration (f, true);
    }
    else
    {
#ifdef CGAL_LINKED_WITH_OPENCV
      m_random_forest->load_configuration (filename);
#endif
    }
  }

  std::size_t& nb_scales() { return m_nb_scales; }
  std::size_t& number_of_trials() { return m_nb_trials; }
  float& smoothing() { return m_smoothing; }
  std::size_t& subdivisions() { return m_subdivisions; }
  const QColor& label_color(std::size_t i) const { return m_label_colors[i]; }
  std::size_t get_label (const char* name)
  {
    for (std::size_t i = 0; i < m_labels.size(); ++ i)
      if (m_labels[i]->name() == name)
        return i;
    return std::size_t(-1);
  }
  void change_label_color (const char* name, const QColor& color)
  {
    m_label_colors[get_label(name)] = color;
  }

protected:

  Label_set m_labels;
  Feature_set m_features;
  std::vector<QColor> m_label_colors;
  Sum_of_weighted_features* m_sowf;
#ifdef CGAL_LINKED_WITH_OPENCV
  Random_forest* m_random_forest;
#endif

  std::size_t m_nb_scales;
  std::size_t m_nb_trials;
  float m_smoothing;
  std::size_t m_subdivisions;
  
};




#endif // ITEM_CLASSIFICATION_BASE_H
