#ifndef POINT_SET_ITEM_CLASSIFICATION_H
#define POINT_SET_ITEM_CLASSIFICATION_H

//#define CGAL_DO_NOT_USE_BOYKOV_KOLMOGOROV_MAXFLOW_SOFTWARE
#define CGAL_CLASSIFICATION_VERBOSE

#include <CGAL/Three/Scene_item.h>

#include "Scene_points_with_normal_item.h"
#include "Item_classification_base.h"
#include "Polyhedron_type_fwd.h"
#include "Kernel_type.h"
#include "Point_set_3.h"

#include <CGAL/Classification.h>

#include <iostream>

#ifdef CGAL_LINKED_WITH_TBB
typedef CGAL::Parallel_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif



// This class represents a point set in the OpenGL scene
class Point_set_item_classification : public Item_classification_base
{
 public:
  typedef Kernel::Point_3 Point_3;
  typedef Kernel::Vector_3 Vector_3;
  typedef CGAL::Classification::RGB_Color Color;
  
  typedef Point_set::Point_map Point_map;
  typedef Point_set::Vector_map Vector_map;

  typedef CGAL::Classification::Point_set_feature_generator<Kernel, Point_set, Point_map>               Generator;
  
 public:
  
  Point_set_item_classification(Scene_points_with_normal_item* points);
  ~Point_set_item_classification();

  CGAL::Three::Scene_item* item() { return m_points; }
  void erase_item() { m_points = NULL; }

  void compute_features ();
  
  void add_selection_to_training_set (const char* name)
  {
    std::size_t label = get_label (name);

    for (Point_set::const_iterator it = m_points->point_set()->first_selected();
         it != m_points->point_set()->end(); ++ it)
      {
        m_training[*it] = label;
        m_classif[*it] = label;
      }

    m_points->resetSelection();
    if (m_index_color == 1 || m_index_color == 2)
      change_color (m_index_color);
  }
  void reset_training_sets()
  {
    for (Point_set::const_iterator it = m_points->point_set()->begin();
         it != m_points->point_set()->end(); ++ it)
      m_training[*it] = std::size_t(-1);
  }
  void validate_selection ()
  {
    for (Point_set::const_iterator it = m_points->point_set()->first_selected();
         it != m_points->point_set()->end(); ++ it)
      m_training[*it] = m_classif[*it];

    m_points->resetSelection();
    if (m_index_color == 1 || m_index_color == 2)
      change_color (m_index_color);
  }
  void train(int predicate);
  bool run (int method, int predicate);

  void change_color (int index);
  void generate_one_item_per_label(std::vector<CGAL::Three::Scene_item*>& items,
                                   const char* name) const
  {
    std::vector<Scene_points_with_normal_item*> points_item
      (m_labels.size(), NULL);
    for (std::size_t i = 0; i < m_labels.size(); ++ i)
      {
        points_item[i] = new Scene_points_with_normal_item;
        points_item[i]->setName (QString("%1 (%2)").arg(name).arg(m_labels[i]->name().c_str()));
        points_item[i]->setColor (m_label_colors[i]);
        items.push_back (points_item[i]);
      }

    for (Point_set::const_iterator it = m_points->point_set()->begin();
         it != m_points->point_set()->end(); ++ it)
      {
        std::size_t c = m_classif[*it];
        if (c != std::size_t(-1))
          points_item[c]->point_set()->insert (m_points->point_set()->point(*it));
      }
  }
  
  bool write_output(std::ostream& out);

  int real_index_color() const;
  void reset_indices();
  void backup_existing_colors_and_add_new();
  void reset_colors();

 private:
  
  template <typename Predicate>
  bool run (int method, const Predicate& predicate)
  {
    std::vector<std::size_t> indices;

    if (method == 0)
      CGAL::Classification::classify<Concurrency_tag> (*(m_points->point_set()),
                                                       m_labels, predicate,
                                                       indices);
    else if (method == 1)
      CGAL::Classification::classify_with_local_smoothing<Concurrency_tag>
        (*(m_points->point_set()), m_points->point_set()->point_map(), m_labels, predicate,
         m_generator->neighborhood().range_neighbor_query(m_generator->radius_neighbors()),
         indices);
    else if (method == 2)
      CGAL::Classification::classify_with_graphcut<Concurrency_tag>
        (*(m_points->point_set()), m_points->point_set()->point_map(),
         m_labels, predicate,
         m_generator->neighborhood().k_neighbor_query(12),
         m_smoothing, m_subdivisions, indices);

    std::vector<std::size_t> ground_truth(m_points->point_set()->size(), std::size_t(-1));
    for (Point_set::const_iterator it = m_points->point_set()->begin();
         it != m_points->point_set()->first_selected(); ++ it)
      {
        m_classif[*it] = indices[*it];
        ground_truth[*it] = m_training[*it];
      }
  
    if (m_index_color == 1 || m_index_color == 2)
      change_color (m_index_color);

    std::cerr << "Precision, recall, F1 scores and IoU:" << std::endl;
    
    CGAL::Classification::Evaluation eval (m_labels, ground_truth, indices);
  
    for (std::size_t i = 0; i < m_labels.size(); ++ i)
      {
        std::cerr << " * " << m_labels[i]->name() << ": "
                  << eval.precision(m_labels[i]) << " ; "
                  << eval.recall(m_labels[i]) << " ; "
                  << eval.f1_score(m_labels[i]) << " ; "
                  << eval.intersection_over_union(m_labels[i]) << std::endl;
      }

    std::cerr << "Accuracy = " << eval.accuracy() << std::endl
              << "Mean F1 score = " << eval.mean_f1_score() << std::endl
              << "Mean IoU = " << eval.mean_intersection_over_union() << std::endl;


    return true;
  }

  Scene_points_with_normal_item* m_points;

  Point_set::Property_map<unsigned char> m_red;
  Point_set::Property_map<unsigned char> m_green;
  Point_set::Property_map<unsigned char> m_blue;
  Point_set::Property_map<Color> m_color;
  Point_set::Property_map<std::size_t> m_training;
  Point_set::Property_map<std::size_t> m_classif;

  Generator* m_generator;
  
  int m_index_color;
  
}; // end class Point_set_item_classification




#endif // POINT_SET_ITEM_CLASSIFICATION_H
