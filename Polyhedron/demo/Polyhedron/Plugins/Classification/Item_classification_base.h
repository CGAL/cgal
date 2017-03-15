#ifndef ITEM_CLASSIFICATION_BASE_H
#define ITEM_CLASSIFICATION_BASE_H

#include <CGAL/Three/Scene_item.h>

#include <QComboBox>

#include <CGAL/Classification/Trainer.h>
#include <CGAL/Classification/Feature_base.h>

class Item_classification_base
{
public:
  typedef CGAL::Classification::Label_handle   Label_handle;
  typedef CGAL::Classification::Feature_handle Feature_handle;

public:
  
  Item_classification_base() { }
  virtual ~Item_classification_base() { }

  virtual CGAL::Three::Scene_item* item() = 0;
  virtual void erase_item() = 0;

  virtual void compute_features () = 0;
  virtual bool features_computed() const = 0;
  virtual std::size_t number_of_features() const = 0;
  virtual Feature_handle feature(std::size_t i) = 0;

  virtual void add_new_label (const char* name, const QColor& color) = 0;
  virtual void remove_label (const char* name) = 0;
  
  virtual void add_selection_to_training_set (const char* name) = 0;
  virtual void reset_training_sets() = 0;
  virtual void validate_selection () = 0;
  virtual void train() = 0;
  virtual bool run (int method) = 0;
  
  virtual void change_color (int index) = 0;
  virtual void fill_display_combo_box (QComboBox* cb, QComboBox* cb1) const = 0;
  virtual void generate_one_item_per_label(std::vector<CGAL::Three::Scene_item*>& items,
                                           const char* name) const = 0;

  virtual bool write_output(std::ostream& out) = 0;
  virtual void save_config(const char* filename) = 0;
  virtual void load_config(const char* filename) = 0;

  std::size_t& nb_scales() { return m_nb_scales; }
  std::size_t& number_of_trials() { return m_nb_trials; }
  double& smoothing() { return m_smoothing; }
  std::size_t& subdivisions() { return m_subdivisions; }
  std::vector<std::pair<Label_handle, QColor> >& labels() { return m_labels; }
  void change_label_color (const char* name, const QColor& color)
  {
    for (std::size_t i = 0; i < m_labels.size(); ++ i)
      if (m_labels[i].first->name() == name)
        {
          m_labels[i].second = color;
          break;
        }
  }

protected:

  std::size_t m_nb_scales;
  std::vector<std::pair<Label_handle, QColor> > m_labels;
  std::size_t m_nb_trials;
  double m_smoothing;
  std::size_t m_subdivisions;
  
};




#endif // ITEM_CLASSIFICATION_BASE_H
