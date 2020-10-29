#ifndef CLUSTER_CLASSIFICATION_H
#define CLUSTER_CLASSIFICATION_H

//#define CGAL_DO_NOT_USE_BOYKOV_KOLMOGOROV_MAXFLOW_SOFTWARE
#define CGAL_CLASSIFICATION_VERBOSE

#include <CGAL/Three/Scene_item.h>

#include "Scene_points_with_normal_item.h"
#include "Item_classification_base.h"
#include "Kernel_type.h"
#include "Point_set_3.h"

#include <CGAL/Classification.h>

#include <iostream>

typedef CGAL::Parallel_if_available_tag Concurrency_tag;

class Cluster_classification : public Item_classification_base
{
 public:
  typedef Kernel::Point_3 Point_3;
  typedef Kernel::Vector_3 Vector_3;

  typedef Point_set::Point_map Point_map;
  typedef Point_set::Vector_map Vector_map;

  typedef CGAL::Classification::Point_set_feature_generator<Kernel, Point_set, Point_map>               Generator;
  typedef CGAL::Classification::Local_eigen_analysis Local_eigen_analysis;
  typedef CGAL::Classification::Feature::Cluster_mean_of_feature Mean_of_feature;
  typedef CGAL::Classification::Feature::Cluster_variance_of_feature Variance_of_feature;

  typedef CGAL::Classification::Cluster<Point_set, Point_map> Cluster;

  struct Point_set_with_cluster_info
    : public CGAL::cpp98::unary_function<const Point_set::Index&,
                                         std::pair<Kernel::Point_3, int> >
  {
    Point_set* point_set;
    Point_set::Property_map<int>* cluster_id;

    Point_set_with_cluster_info (Point_set* point_set,
                                 Point_set::Property_map<int>& cluster_id)
      : point_set (point_set)
      , cluster_id (&cluster_id)
    { }

    std::pair<Kernel::Point_3, int> operator() (const Point_set::Index& idx) const
    {
      return std::make_pair (point_set->point(idx), (*cluster_id)[idx]);
    }

  };

 public:

  Cluster_classification(Scene_points_with_normal_item* points);
  ~Cluster_classification();

  CGAL::Three::Scene_item* item() { return m_points; }
  void erase_item() { m_points = NULL; }

  CGAL::Bbox_3 bbox()
  {
    if (m_points->point_set()->nb_selected_points() == 0)
      return m_points->bbox();

    CGAL::Bbox_3 bb = CGAL::bbox_3 (boost::make_transform_iterator
                                    (m_points->point_set()->first_selected(),
                                     CGAL::Property_map_to_unary_function<Point_set::Point_map>
                                     (m_points->point_set()->point_map())),
                                    boost::make_transform_iterator
                                    (m_points->point_set()->end(),
                                     CGAL::Property_map_to_unary_function<Point_set::Point_map>
                                     (m_points->point_set()->point_map())));

    double xcenter = (bb.xmax() + bb.xmin()) / 2.;
    double ycenter = (bb.ymax() + bb.ymin()) / 2.;
    double zcenter = (bb.zmax() + bb.zmin()) / 2.;

    double dx = bb.xmax() - bb.xmin();
    double dy = bb.ymax() - bb.ymin();
    double dz = bb.zmax() - bb.zmin();

    dx *= 10.;
    dy *= 10.;
    dz *= 10.;

    return CGAL::Bbox_3 (xcenter - dx, ycenter - dy, zcenter - dz,
                         xcenter + dx, ycenter + dy, zcenter + dz);
  }

  void compute_features (std::size_t nb_scales, float voxel_size);
  void add_remaining_point_set_properties_as_features(Feature_set& feature_set);

  void select_random_region();

  void add_cluster_features ()
  {
    m_eigen = boost::make_shared<Local_eigen_analysis>
      (Local_eigen_analysis::create_from_point_clusters(m_clusters,
                                                        Concurrency_tag()));

    m_features.template add<CGAL::Classification::Feature::Cluster_size> (m_clusters);
    m_features.template add<CGAL::Classification::Feature::Cluster_vertical_extent> (m_clusters);

    m_features.template add<CGAL::Classification::Feature::Eigenvalue> (m_clusters, *m_eigen, 0);
    m_features.template add<CGAL::Classification::Feature::Eigenvalue> (m_clusters, *m_eigen, 1);
    m_features.template add<CGAL::Classification::Feature::Eigenvalue> (m_clusters, *m_eigen, 2);

    // m_features.template add<CGAL::Classification::Feature::Linearity> (m_clusters, *m_eigen);
    // m_features.template add<CGAL::Classification::Feature::Planarity> (m_clusters, *m_eigen);
    // m_features.template add<CGAL::Classification::Feature::Sphericity> (m_clusters, *m_eigen);
    // m_features.template add<CGAL::Classification::Feature::Omnivariance> (m_clusters, *m_eigen);
    // m_features.template add<CGAL::Classification::Feature::Anisotropy> (m_clusters, *m_eigen);
    // m_features.template add<CGAL::Classification::Feature::Eigentropy> (m_clusters, *m_eigen);
    // m_features.template add<CGAL::Classification::Feature::Sum_eigenvalues> (m_clusters, *m_eigen);
    // m_features.template add<CGAL::Classification::Feature::Surface_variation> (m_clusters, *m_eigen);
    // m_features.template add<CGAL::Classification::Feature::Verticality<Kernel> > (m_clusters, *m_eigen);
  }

  template <typename Type>
  bool try_adding_simple_feature (Feature_set& feature_set, const std::string& name)
  {
    typedef typename Point_set::template Property_map<Type> Pmap;
    bool okay = false;
    Pmap pmap;
    boost::tie (pmap, okay) = m_points->point_set()->template property_map<Type>(name.c_str());
    if (okay)
      feature_set.template add<CGAL::Classification::Feature::Simple_feature <Point_set, Pmap> >
        (*(m_points->point_set()), pmap, name.c_str());

    return okay;
  }

  void add_selection_to_training_set (std::size_t label)
  {
    for (Point_set::const_iterator it = m_points->point_set()->first_selected();
         it != m_points->point_set()->end(); ++ it)
    {
      int cid = m_cluster_id[*it];
      if (cid != -1)
      {
        m_clusters[cid].training() = int(label);
        m_clusters[cid].label() = int(label);
      }
    }

    m_points->resetSelection();
    if (m_index_color == 1 || m_index_color == 2)
      change_color (m_index_color);
  }
  void reset_training_set(std::size_t label)
  {
    for (std::size_t i = 0; i < m_clusters.size(); ++ i)
      if (m_clusters[i].training() == int(label))
        m_clusters[i].training() = -1;

    if (m_index_color == 1 || m_index_color == 2)
      change_color (m_index_color);
  }
  void reset_training_set_of_selection()
  {
    for (Point_set::const_iterator it = m_points->point_set()->first_selected();
         it != m_points->point_set()->end(); ++ it)
    {
      int cid = m_cluster_id[*it];
      if (cid != -1)
      {
        m_clusters[cid].training() = -1;
        m_clusters[cid].label() = -1;
      }
    }
    if (m_index_color == 1 || m_index_color == 2)
      change_color (m_index_color);
  }
  void reset_training_sets()
  {
    for (std::size_t i = 0; i < m_clusters.size(); ++ i)
      m_clusters[i].training() = -1;
    if (m_index_color == 1 || m_index_color == 2)
      change_color (m_index_color);
  }
  void validate_selection ()
  {
    for (Point_set::const_iterator it = m_points->point_set()->first_selected();
         it != m_points->point_set()->end(); ++ it)
    {
      int cid = m_cluster_id[*it];
      if (cid != -1)
        m_clusters[cid].training() = m_clusters[cid].label();
    }

    m_points->resetSelection();
    if (m_index_color == 1 || m_index_color == 2)
      change_color (m_index_color);
  }
  void train(int classifier, const QMultipleInputDialog& dialog);
  bool run (int method, int classifier, std::size_t subdivisions, double smoothing);

  void update_color () { change_color (m_index_color); }
  void change_color (int index, float* vmin = NULL, float* vmax = NULL);
  CGAL::Three::Scene_item* generate_one_item (const char* name,
                                              int label) const
  {
    Scene_points_with_normal_item* points_item
      = new Scene_points_with_normal_item;

    points_item->setName (QString("%1 (%2)").arg(name).arg(m_labels[label]->name().c_str()));
    points_item->setColor (label_qcolor (m_labels[label]));
    for (Point_set::const_iterator it = m_points->point_set()->begin();
         it != m_points->point_set()->end(); ++ it)
    {
      int cid = m_cluster_id[*it];
      if (cid != -1)
      {
        int c = m_clusters[cid].label();
        if (c == label)
          points_item->point_set()->insert (m_points->point_set()->point(*it));
      }
    }
    return points_item;
  }
  void generate_one_item_per_label(std::vector<CGAL::Three::Scene_item*>& items,
                                   const char* name) const
  {
    std::vector<Scene_points_with_normal_item*> points_item
      (m_labels.size(), NULL);
    for (std::size_t i = 0; i < m_labels.size(); ++ i)
    {
      points_item[i] = new Scene_points_with_normal_item;
      points_item[i]->setName (QString("%1 (%2)").arg(name).arg(m_labels[i]->name().c_str()));
      points_item[i]->setColor (label_qcolor (m_labels[i]));
      items.push_back (points_item[i]);
    }

    for (Point_set::const_iterator it = m_points->point_set()->begin();
         it != m_points->point_set()->end(); ++ it)
    {
      int cid = m_cluster_id[*it];
      if (cid != -1)
      {
        int c = m_clusters[cid].label();
        points_item[c]->point_set()->insert (m_points->point_set()->point(*it));
      }
    }

  }

  QColor add_new_label (const char* name)
  {
    QColor out = Item_classification_base::add_new_label (name);
    update_comments_of_point_set_item();
    return out;
  }

  void remove_label (std::size_t position)
  {
    Item_classification_base::remove_label (position);

    for (std::size_t i = 0; i < m_clusters.size(); ++ i)
    {
      if (m_clusters[i].training() == int(position))
        m_clusters[i].training() = -1;
      else if (m_clusters[i].training() > int(position))
        m_clusters[i].training() --;

      if (m_clusters[i].label() == int(position))
        m_clusters[i].label() = -1;
      else if (m_clusters[i].label() > int(position))
        m_clusters[i].label() --;
    }

    update_comments_of_point_set_item();
  }

  void fill_display_combo_box (QComboBox* cb, QComboBox* cb1) const
  {
    cb->addItem ("Clusters");
    Item_classification_base::fill_display_combo_box(cb, cb1);
  }

  int real_index_color() const;
  void reset_indices();
  void backup_existing_colors_and_add_new();
  void reset_colors();

 private:

  void update_comments_of_point_set_item()
  {
    std::string& comments = m_points->comments();

    // Remove previously registered labels from comments
    std::string new_comment;

    std::istringstream stream (comments);
    std::string line;
    while (getline(stream, line))
    {
      std::stringstream iss (line);
      std::string tag;
      if (iss >> tag && tag == "label")
        continue;
      new_comment += line + "\n";
    }
    comments = new_comment;

    comments += "label -1 unclassified\n";
    for (std::size_t i = 0; i < m_labels.size(); ++ i)
    {
      std::ostringstream oss;
      oss << "label " << i << " " << m_labels[i]->name() << std::endl;
      comments += oss.str();
    }
  }

  template <typename Classifier>
  bool run (int method, const Classifier& classifier,
            std::size_t subdivisions, double smoothing)
  {
    std::vector<int> indices (m_clusters.size(), -1);

    if (method == 0)
      CGAL::Classification::classify<Concurrency_tag> (m_clusters,
                                                       m_labels, classifier,
                                                       indices);
    else if (method == 1)
      CGAL::Classification::classify_with_local_smoothing<Concurrency_tag>
        (m_clusters,
         CGAL::Identity_property_map<Cluster>(),
         m_labels, classifier,
         Cluster::Neighbor_query(),
         indices);
    else if (method == 2)
      CGAL::Classification::classify_with_graphcut<Concurrency_tag>
        (m_clusters,
         CGAL::Identity_property_map<Cluster>(),
         m_labels, classifier,
         Cluster::Neighbor_query(),
         smoothing, subdivisions, indices);

    std::vector<int> ground_truth(m_clusters.size(), -1);
    for (std::size_t i = 0; i < m_clusters.size(); ++ i)
    {
      m_clusters[i].label() = indices[i];
      ground_truth[i] = m_clusters[i].training();
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

  std::vector<Cluster> m_clusters;

  Point_set::Property_map<unsigned char> m_red;
  Point_set::Property_map<unsigned char> m_green;
  Point_set::Property_map<unsigned char> m_blue;
  Point_set::Property_map<CGAL::Color> m_color;
  Point_set::Property_map<int> m_cluster_id;
  Point_set::Property_map<int> m_training;
  Point_set::Property_map<int> m_classif;

  std::vector<std::vector<float> > m_label_probabilities;

  int m_index_color;

  boost::shared_ptr<Local_eigen_analysis> m_eigen;

  bool m_input_is_las;

}; // end class Cluster_classification




#endif // CLUSTER_CLASSIFICATION_H
