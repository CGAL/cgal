#ifndef ITEM_CLASSIFICATION_BASE_H
#define ITEM_CLASSIFICATION_BASE_H

#include <CGAL/Three/Scene_item.h>

#include <QComboBox>
#include <QLineEdit>
#include <QSpinBox>
#include <QMultipleInputDialog.h>

#include <CGAL/Classification/Feature_set.h>
#include <CGAL/Classification/Label_set.h>
#include <CGAL/Classification/Sum_of_weighted_features_classifier.h>
#include <CGAL/Classification/ETHZ/Random_forest_classifier.h>

#ifdef CGAL_LINKED_WITH_OPENCV
#include <CGAL/Classification/OpenCV/Random_forest_classifier.h>
#endif
#ifdef CGAL_LINKED_WITH_TENSORFLOW
#include <CGAL/Classification/TensorFlow/Neural_network_classifier.h>
#endif

class Item_classification_base
{
public:
  typedef CGAL::Classification::Label_handle   Label_handle;
  typedef CGAL::Classification::Feature_handle Feature_handle;
  typedef CGAL::Classification::Label_set   Label_set;
  typedef CGAL::Classification::Feature_set Feature_set;
  typedef CGAL::Classification::Sum_of_weighted_features_classifier Sum_of_weighted_features;
  typedef CGAL::Classification::ETHZ::Random_forest_classifier ETHZ_random_forest;

#ifdef CGAL_LINKED_WITH_OPENCV
  typedef CGAL::Classification::OpenCV::Random_forest_classifier Random_forest;
#endif
#ifdef CGAL_LINKED_WITH_TENSORFLOW
  typedef CGAL::Classification::TensorFlow::Neural_network_classifier<> Neural_network;
#endif
  
public:
  
  Item_classification_base() { }
  virtual ~Item_classification_base() { }

  virtual CGAL::Three::Scene_item* item() = 0;
  virtual void erase_item() = 0;

  virtual CGAL::Bbox_3 bbox() { return item()->bbox(); }

  virtual void compute_features (std::size_t nb_scales, float voxel_size) = 0;

  virtual std::string feature_statistics () const { return std::string(); }
  
  virtual void add_selection_to_training_set (std::size_t label) = 0;
  virtual void reset_training_set (std::size_t label) = 0;
  virtual void reset_training_set_of_selection() = 0;
  virtual void reset_training_sets() = 0;

  virtual void select_random_region() = 0;
  virtual void validate_selection () = 0;
  virtual void train(int classifier, const QMultipleInputDialog&) = 0;
  virtual bool run (int method, int classifier, std::size_t subdivisions, double smoothing) = 0;
  
  virtual void update_color () = 0;
  virtual void change_color (int index, float* vmin = NULL, float* vmax = NULL) = 0;
  virtual CGAL::Three::Scene_item* generate_one_item (const char* name,
                                                      int label) const = 0;
  virtual void generate_one_item_per_label(std::vector<CGAL::Three::Scene_item*>& items,
                                           const char* name) const = 0;

  bool features_computed() const { return (m_features.size() != 0); }
  std::size_t number_of_features() const { return m_features.size(); }  
  Feature_handle feature(std::size_t i) { return m_features[i]; }
  float weight (Feature_handle f) const { return m_sowf->weight(f); }
  void set_weight (Feature_handle f, float w) const { m_sowf->set_weight(f,w); }
  Sum_of_weighted_features::Effect effect (Label_handle l, Feature_handle f) const { return m_sowf->effect(l,f); }
  void set_effect (Label_handle l, Feature_handle f, Sum_of_weighted_features::Effect e)
  { m_sowf->set_effect (l, f, e); }
  
  virtual QColor add_new_label (const char* name)
  {
    m_labels.add(name);
    m_label_colors.push_back (get_new_label_color (name));
    
    delete m_sowf;
    m_sowf = new Sum_of_weighted_features (m_labels, m_features);

    delete m_ethz;
    m_ethz = new ETHZ_random_forest (m_labels, m_features);

#ifdef CGAL_LINKED_WITH_OPENCV
    delete m_random_forest;
    m_random_forest = new Random_forest (m_labels, m_features);
#endif

#ifdef CGAL_LINKED_WITH_TENSORFLOW
    delete m_neural_network;
    m_neural_network = new Neural_network (m_labels, m_features);
#endif
    
    return m_label_colors.back();
  }
  virtual void remove_label (std::size_t position)
  {
    m_labels.remove(m_labels[position]);
    m_label_colors.erase (m_label_colors.begin() + position);

    delete m_sowf;
    m_sowf = new Sum_of_weighted_features (m_labels, m_features);

    delete m_ethz;
    m_ethz = new ETHZ_random_forest (m_labels, m_features);

#ifdef CGAL_LINKED_WITH_OPENCV
    delete m_random_forest;
    m_random_forest = new Random_forest (m_labels, m_features);
#endif

#ifdef CGAL_LINKED_WITH_TENSORFLOW
    delete m_neural_network;
    m_neural_network = new Neural_network (m_labels, m_features);
#endif
  }
  
  virtual void clear_labels ()
  {
    m_labels.clear();
    m_label_colors.clear();

    delete m_sowf;
    m_sowf = new Sum_of_weighted_features (m_labels, m_features);

    delete m_ethz;
    m_ethz = new ETHZ_random_forest (m_labels, m_features);

#ifdef CGAL_LINKED_WITH_OPENCV
    delete m_random_forest;
    m_random_forest = new Random_forest (m_labels, m_features);
#endif

#ifdef CGAL_LINKED_WITH_TENSORFLOW
    delete m_neural_network;
    m_neural_network = new Neural_network (m_labels, m_features);
#endif
  }
  std::size_t number_of_labels() const { return m_labels.size(); }
  Label_handle label(std::size_t i) { return m_labels[i]; }

  virtual void fill_display_combo_box (QComboBox* cb, QComboBox* cb1) const
  {
    for (std::size_t i = 0; i < m_labels.size(); ++ i)
      {
        std::ostringstream oss;
        oss << "Label " << m_labels[i]->name();
        cb->addItem (oss.str().c_str());
        cb1->addItem (oss.str().c_str());
      }
    for (std::size_t i = 0; i < m_features.size(); ++ i)
      {
        std::ostringstream oss;
        oss << "Feature " << m_features[i]->name();
        cb->addItem (oss.str().c_str());
        cb1->addItem (oss.str().c_str());
      }
  }

  void save_config(const char* filename, int classifier)
  {
    if (m_features.size() == 0)
      {
        std::cerr << "Error: features not computed" << std::endl;
        return;
      }

    if (classifier == 0)
    {
      std::ofstream f (filename);
      m_sowf->save_configuration (f);
    }
    else if (classifier == 1)
    {
      std::ofstream f (filename, std::ios_base::out | std::ios_base::binary);
      m_ethz->save_configuration (f);
    }
    else if (classifier == 2)
    {
#ifdef CGAL_LINKED_WITH_OPENCV
      m_random_forest->save_configuration (filename);
#endif
    }
    else if (classifier == 3)
    {
#ifdef CGAL_LINKED_WITH_TENSORFLOW
      std::ofstream f (filename);
      m_neural_network->save_configuration (f);
#endif
    }
  }
  void load_config(const char* filename, int classifier)
  {
    if (m_features.size() == 0)
      {
        std::cerr << "Error: features not computed" << std::endl;
        return;
      }

    if (classifier == 0)
    {
      std::ifstream f (filename);
      m_sowf->load_configuration (f, true);
    }
    else if (classifier == 1)
    {
      if (m_ethz == NULL)
        m_ethz = new ETHZ_random_forest (m_labels, m_features);
      std::ifstream f (filename, std::ios_base::in | std::ios_base::binary);
      m_ethz->load_configuration (f);
    }
    else if (classifier == 2)
    {
#ifdef CGAL_LINKED_WITH_OPENCV
      m_random_forest->load_configuration (filename);
#endif
    }
    else if (classifier == 3)
    {
#ifdef CGAL_LINKED_WITH_TENSORFLOW
      if (m_neural_network == NULL)
        m_neural_network = new Neural_network (m_labels, m_features);
      std::ifstream f (filename);
      m_neural_network->load_configuration (f, true);
#endif
    }
  }

  const QColor& label_color(std::size_t i) const { return m_label_colors[i]; }
  void change_label_color (std::size_t position, const QColor& color)
  {
    m_label_colors[position] = color;
  }
  void change_label_name (std::size_t position, const std::string& name)
  {
    m_labels[position]->set_name (name);
  }

  QColor get_new_label_color (const std::string& name)
  {
    QColor color (64 + rand() % 192,
                  64 + rand() % 192,
                  64 + rand() % 192);
      
    if (name == "ground")
      color = QColor (186, 189, 182);
    else if (name == "low_veget")
      color = QColor (78, 154, 6);
    else if (name == "med_veget"
             || name == "vegetation")
      color = QColor (138, 226, 52);
    else if (name == "high_veget")
      color = QColor (204, 255, 201);
    else if (name == "building"
             || name == "roof")
      color = QColor (245, 121, 0);
    else if (name == "noise")
      color = QColor (0, 0, 0);
    else if (name == "reserved")
      color = QColor (233, 185, 110);
    else if (name == "water")
      color = QColor (114, 159, 207);
    else if (name == "rail")
      color = QColor (136, 46, 25);
    else if (name == "road_surface")
      color = QColor (56, 56, 56);
    else if (name == "reserved_2")
      color = QColor (193, 138, 51);
    else if (name == "wire_guard")
      color = QColor (37, 61, 136);
    else if (name == "wire_conduct")
      color = QColor (173, 127, 168);
    else if (name == "trans_tower")
      color = QColor (136, 138, 133);
    else if (name == "wire_connect")
      color = QColor (145, 64, 236);
    else if (name == "bridge_deck")
      color = QColor (213, 93, 93);
    else if (name == "high_noise")
      color = QColor (255, 0, 0);
    else if (name == "facade")
      color = QColor (77, 131, 186);
    
    return color;
  }

protected:

  Label_set m_labels;
  Feature_set m_features;
  std::vector<QColor> m_label_colors;
  Sum_of_weighted_features* m_sowf;
  ETHZ_random_forest* m_ethz;
#ifdef CGAL_LINKED_WITH_OPENCV
  Random_forest* m_random_forest;
#endif
#ifdef CGAL_LINKED_WITH_TENSORFLOW
  Neural_network* m_neural_network;
#endif
  
};




#endif // ITEM_CLASSIFICATION_BASE_H
