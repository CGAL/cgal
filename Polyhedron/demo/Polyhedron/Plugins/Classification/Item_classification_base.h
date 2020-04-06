#ifndef ITEM_CLASSIFICATION_BASE_H
#define ITEM_CLASSIFICATION_BASE_H

#include <CGAL/Three/Scene_item.h>

#include <QComboBox>
#include <QLineEdit>
#include <QSpinBox>
#include <QMultipleInputDialog.h>

#include <CGAL/Classification/Feature_set.h>
#include <CGAL/Classification/Label_set.h>


#include <CGAL/Classification/ETHZ/Random_forest_classifier.h>
#include <CGAL/Classification/Sum_of_weighted_features_classifier.h>

#ifdef CGAL_LINKED_WITH_TENSORFLOW
#  include <CGAL/Classification/TensorFlow/Neural_network_classifier.h>
#endif

#ifdef CGAL_LINKED_WITH_OPENCV
#  include <CGAL/Classification/OpenCV/Random_forest_classifier.h>
#endif

#define CGAL_CLASSIFICATION_ETHZ_ID "Random Forest (ETHZ)"
#define CGAL_CLASSIFICATION_ETHZ_NUMBER 0

#define CGAL_CLASSIFICATION_TENSORFLOW_ID "Neural Network (TensorFlow)"
#define CGAL_CLASSIFICATION_TENSORFLOW_NUMBER 1

#define CGAL_CLASSIFICATION_OPENCV_ID "Random Forest (OpenCV)"
#define CGAL_CLASSIFICATION_OPENCV_NUMBER 2

#define CGAL_CLASSIFICATION_SOWF_ID "Sum of Weighted Features"
#define CGAL_CLASSIFICATION_SOWF_NUMBER 3


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

  QColor label_qcolor (Label_handle l) const
  {
    return QColor (l->color().red(), l->color().green(), l->color().blue());
  }

  virtual QColor add_new_label (const char* name)
  {
    m_labels.add(name);

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

    return label_qcolor (m_labels[m_labels.size() - 1]);
  }
  virtual void remove_label (std::size_t position)
  {
    m_labels.remove(m_labels[position]);

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

    if (classifier == CGAL_CLASSIFICATION_SOWF_NUMBER)
    {
      std::ofstream f (filename);
      m_sowf->save_configuration (f);
    }
    else if (classifier == CGAL_CLASSIFICATION_ETHZ_NUMBER)
    {
      std::cerr << "D ";
      std::ofstream f (filename, std::ios_base::out | std::ios_base::binary);
      m_ethz->save_configuration (f);
      std::cerr << "E ";
    }
    else if (classifier == CGAL_CLASSIFICATION_OPENCV_NUMBER)
    {
#ifdef CGAL_LINKED_WITH_OPENCV
      m_random_forest->save_configuration (filename);
#endif
    }
    else if (classifier == CGAL_CLASSIFICATION_TENSORFLOW_NUMBER)
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

    if (classifier == CGAL_CLASSIFICATION_SOWF_NUMBER)
    {
      std::ifstream f (filename);
      m_sowf->load_configuration (f, true);
    }
    else if (classifier == CGAL_CLASSIFICATION_ETHZ_NUMBER)
    {
      if (m_ethz == NULL)
        m_ethz = new ETHZ_random_forest (m_labels, m_features);
      std::ifstream f (filename, std::ios_base::in | std::ios_base::binary);

      // Handle deprecated files
      if (std::string(filename).find(".gz") != std::string::npos)
        m_ethz->load_deprecated_configuration(f);
      else
        m_ethz->load_configuration (f);
    }
    else if (classifier == CGAL_CLASSIFICATION_OPENCV_NUMBER)
    {
#ifdef CGAL_LINKED_WITH_OPENCV
      m_random_forest->load_configuration (filename);
#endif
    }
    else if (classifier == CGAL_CLASSIFICATION_TENSORFLOW_NUMBER)
    {
#ifdef CGAL_LINKED_WITH_TENSORFLOW
      if (m_neural_network == NULL)
        m_neural_network = new Neural_network (m_labels, m_features);
      std::ifstream f (filename);
      m_neural_network->load_configuration (f, true);
#endif
    }
  }

  QColor label_color(std::size_t i) const
  {
    return label_qcolor (m_labels[i]);
  }
  void change_label_color (std::size_t position, const QColor& color)
  {
    m_labels[position]->set_color
      (CGAL::Color (color.red(), color.green(), color.blue()));
  }
  void change_label_name (std::size_t position, const std::string& name)
  {
    m_labels[position]->set_name (name);
  }

protected:

  Label_set m_labels;
  Feature_set m_features;
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
