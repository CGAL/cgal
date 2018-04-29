#ifndef CLUSTER_CLASSIFICATION_H
#define CLUSTER_CLASSIFICATION_H

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


class Cluster_classification : public Item_classification_base
{
 public:
  typedef Kernel::Point_3 Point_3;
  typedef Kernel::Vector_3 Vector_3;
  typedef CGAL::Classification::RGB_Color Color;
  
  typedef Point_set::Point_map Point_map;
  typedef Point_set::Vector_map Vector_map;

  typedef CGAL::Classification::Point_set_feature_generator<Kernel, Point_set, Point_map>               Generator;
  typedef CGAL::Classification::Local_eigen_analysis Local_eigen_analysis;

  struct Cluster
  {
    std::vector<Point_set::Index> inliers;
    std::vector<std::size_t> neighbors;
    CGAL::Bbox_3 bounding_box;
    int training;
    int label;
    Cluster() : training (-1), label (-1) { }

    std::size_t size() const { return inliers.size(); }
    const Point_set::Index& operator[] (std::size_t i) const { return inliers[i]; }

    CGAL::Bbox_3 bbox() const
    {
      return bounding_box;
    }
  };

  struct Point_set_with_cluster_info
    : public CGAL::unary_function<const Point_set::Index&,
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

  struct Neighbor_query
  {
    template <typename OutputIterator>
    OutputIterator operator() (const Cluster& cluster, OutputIterator output) const
    {
      return std::copy (cluster.neighbors.begin(), cluster.neighbors.end(), output);
    }
  };
  
  class Mean_feature : public CGAL::Classification::Feature_base
  {
    std::vector<float> m_values;
    
  public:

    template <typename PointRange,
              typename ClusterRange>
    Mean_feature (const PointRange&,
                  ClusterRange& clusters,
                  Feature_handle pointwise_feature)
    {
      std::ostringstream oss;
      oss << "mean_" << pointwise_feature->name();
      this->set_name (oss.str());

      m_values.reserve (clusters.size());
      for (std::size_t i = 0; i < clusters.size(); ++ i)
      {
        double mean = 0.;

        for (std::size_t j = 0; j < clusters[i].inliers.size(); ++ j)
          mean += double(pointwise_feature->value (clusters[i].inliers[j]));
        mean /= clusters[i].inliers.size();
        m_values.push_back (float(mean));
      }
    }

    virtual float value (std::size_t cluster_index)
    {
      return m_values[cluster_index];
    }
    
  };

  class Standard_deviation_feature : public CGAL::Classification::Feature_base
  {
    std::vector<float> m_values;
    
  public:

    template <typename PointRange,
              typename ClusterRange>
    Standard_deviation_feature (const PointRange&,
                                ClusterRange& clusters,
                                Feature_handle pointwise_feature,
                                Feature_handle mean_feature)
    {
      std::ostringstream oss;
      oss << "standard_dev_" << pointwise_feature->name();
      this->set_name (oss.str());

      m_values.reserve (clusters.size());
      for (std::size_t i = 0; i < clusters.size(); ++ i)
      {
        double mean = double (mean_feature->value(i));
        double variance = 0.;

        for (std::size_t j = 0; j < clusters[i].inliers.size(); ++ j)
        {
          double v = double (pointwise_feature->value (clusters[i].inliers[j]));
          variance += (v - mean) * (v - mean);
        }
        variance /= clusters[i].inliers.size();
        m_values.push_back (float(std::sqrt (variance)));
      }
    }

    virtual float value (std::size_t cluster_index)
    {
      return m_values[cluster_index];
    }
    
  };

  class Cluster_size_feature : public CGAL::Classification::Feature_base
  {
    std::vector<float> m_values;
  public:

    template <typename ClusterRange>
    Cluster_size_feature (const ClusterRange& input)
    {
      this->set_name ("cluster_size");
      m_values.reserve (input.size());

      for (std::size_t i = 0; i < input.size(); ++ i)
        m_values.push_back (float(input[i].size()));

      std::vector<float> c = m_values;
      std::sort (c.begin(), c.end());
      std::cerr << " * Min cluster size = " << c.front() << std::endl
                << " * 10% cluster size = " << c[c.size() / 10] << std::endl
                << " * 25% cluster size = " << c[c.size() / 4] << std::endl
                << " * Median cluster size = " << c[c.size() / 2] << std::endl
                << " * 75% cluster size = " << c[3 * c.size() / 4] << std::endl
                << " * 90% cluster size = " << c[9 * c.size() / 10] << std::endl
                << " * Max cluster size = " << c.back() << std::endl;
    }

    virtual float value (std::size_t cluster_index) { return m_values[cluster_index]; }
  };

  class Vertical_extent_feature : public CGAL::Classification::Feature_base
  {
    std::vector<float> m_values;
  public:

    template <typename ClusterRange>
    Vertical_extent_feature (const ClusterRange& input, Point_set* points)
    {
      this->set_name ("vertical_extent");

      m_values.reserve (input.size());
      for (std::size_t i = 0; i < input.size(); ++ i)
      {
        float min_z = std::numeric_limits<float>::max();
        float max_z = -std::numeric_limits<float>::min();
        
        for (std::size_t j = 0; j < input[i].size(); ++ j)
        {
          const Kernel::Point_3& p = points->point(input[i][j]);
          min_z = (std::min) (float(p.z()), min_z);
          max_z = (std::max) (float(p.z()), max_z);
        }
        m_values.push_back ((max_z - min_z));
      }
    }

    virtual float value (std::size_t cluster_index) { return m_values[cluster_index]; }
  };

  class Density_feature : public CGAL::Classification::Feature_base
  {
    std::vector<float> m_values;
  public:

    template <typename ClusterRange>
    Density_feature (const ClusterRange& input, Point_set* points)
    {
      this->set_name ("density");

      m_values.reserve (input.size());
      for (std::size_t i = 0; i < input.size(); ++ i)
      {
        std::vector<Kernel::Point_3> pts;
        for (std::size_t j = 0; j < input[i].size(); ++ j)
          pts.push_back (points->point (input[i][j]));

        Kernel::Point_3 centroid = CGAL::centroid (pts.begin(), pts.end());

        double max_dist = 0.;
        for (std::size_t j = 0; j < pts.size(); ++ j)
        {
          double dist = CGAL::squared_distance (pts[j], centroid);
          if (dist > max_dist)
            max_dist = dist;
        }
        m_values.push_back (input[i].size() / (CGAL_PI * max_dist));
      }
    }

    virtual float value (std::size_t cluster_index) { return m_values[cluster_index]; }
  };

  struct Cluster_index_to_point_map
  {
    typedef Point_set::Index key_type;
    typedef Kernel::Point_3 value_type;
    typedef value_type& reference;
    typedef boost::readable_property_map_tag category;
    
    Point_set* points;
    
    Cluster_index_to_point_map (Point_set* point_set) : points (point_set) { }

    inline friend reference get (const Cluster_index_to_point_map& pmap, key_type idx)
    {
      return pmap.points->point(idx);
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

  void compute_features (std::size_t nb_scales);
  void add_remaining_point_set_properties_as_features(Feature_set& feature_set);
  
  void select_random_region();

  void add_cluster_features ()
  {
    m_eigen = boost::make_shared<Local_eigen_analysis> (Local_eigen_analysis::Input_is_clusters(),
                                                        m_clusters,
                                                        Cluster_index_to_point_map (m_points->point_set()),
                                                        Concurrency_tag());

    m_features.template add<Cluster_size_feature> (m_clusters);
    m_features.template add<Vertical_extent_feature> (m_clusters, m_points->point_set());
    m_features.template add<Density_feature> (m_clusters, m_points->point_set());
    m_features.template add<CGAL::Classification::Feature::Linearity> (m_clusters, *m_eigen);
    m_features.template add<CGAL::Classification::Feature::Planarity> (m_clusters, *m_eigen);
    m_features.template add<CGAL::Classification::Feature::Sphericity> (m_clusters, *m_eigen);
    m_features.template add<CGAL::Classification::Feature::Omnivariance> (m_clusters, *m_eigen);
    m_features.template add<CGAL::Classification::Feature::Anisotropy> (m_clusters, *m_eigen);
    m_features.template add<CGAL::Classification::Feature::Eigentropy> (m_clusters, *m_eigen);
    m_features.template add<CGAL::Classification::Feature::Sum_eigenvalues> (m_clusters, *m_eigen);
    m_features.template add<CGAL::Classification::Feature::Surface_variation> (m_clusters, *m_eigen);
    m_features.template add<CGAL::Classification::Feature::Verticality<Kernel> > (m_clusters, *m_eigen);
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
        m_clusters[cid].training = int(label);
        m_clusters[cid].label = int(label);
      }
    }

    m_points->resetSelection();
    if (m_index_color == 1 || m_index_color == 2)
      change_color (m_index_color);
  }
  void reset_training_set(std::size_t label)
  {
    for (std::size_t i = 0; i < m_clusters.size(); ++ i)
      if (m_clusters[i].training == int(label))
        m_clusters[i].training = -1;

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
        m_clusters[cid].training = -1;
        m_clusters[cid].label = -1;
      }
    }
    if (m_index_color == 1 || m_index_color == 2)
      change_color (m_index_color);
  }
  void reset_training_sets()
  {
    for (std::size_t i = 0; i < m_clusters.size(); ++ i)
        m_clusters[i].training = -1;
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
        m_clusters[cid].training = m_clusters[cid].label;
    }

    m_points->resetSelection();
    if (m_index_color == 1 || m_index_color == 2)
      change_color (m_index_color);
  }
  void train(int classifier, unsigned int nb_trials,
             std::size_t num_trees, std::size_t max_depth);
  bool run (int method, int classifier, std::size_t subdivisions, double smoothing);

  void update_color () { change_color (m_index_color); }
  void change_color (int index);
  CGAL::Three::Scene_item* generate_one_item (const char* name,
                                              int label) const
  {
    Scene_points_with_normal_item* points_item
      = new Scene_points_with_normal_item;
    
    points_item->setName (QString("%1 (%2)").arg(name).arg(m_labels[label]->name().c_str()));
    points_item->setColor (m_label_colors[label]);
    for (Point_set::const_iterator it = m_points->point_set()->begin();
         it != m_points->point_set()->end(); ++ it)
    {
      int cid = m_cluster_id[*it];
      if (cid != -1)
      {
        int c = m_clusters[cid].label;
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
      points_item[i]->setColor (m_label_colors[i]);
      items.push_back (points_item[i]);
    }

    for (Point_set::const_iterator it = m_points->point_set()->begin();
         it != m_points->point_set()->end(); ++ it)
    {
      int cid = m_cluster_id[*it];
      if (cid != -1)
      {
        int c = m_clusters[cid].label;
        points_item[c]->point_set()->insert (m_points->point_set()->point(*it));
      }
    }

  }
  
  bool write_output(std::ostream& out);

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
      if (m_clusters[i].training == int(position))
        m_clusters[i].training = -1;
      else if (m_clusters[i].training > int(position))
        m_clusters[i].training --;
          
      if (m_clusters[i].label == int(position))
        m_clusters[i].label = -1;
      else if (m_clusters[i].label > int(position))
        m_clusters[i].label --;
    }

    update_comments_of_point_set_item();
  }
  
  void fill_display_combo_box (QComboBox* cb, QComboBox* cb1) const
  {
    cb->addItem ("Clusters");
    for (std::size_t i = 0; i < m_features.size(); ++ i)
      {
        std::ostringstream oss;
        oss << "Feature " << m_features[i]->name();
        cb->addItem (oss.str().c_str());
        cb1->addItem (oss.str().c_str());
      }
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
         Neighbor_query(),
         indices);
    else if (method == 2)
      CGAL::Classification::classify_with_graphcut<Concurrency_tag>
        (m_clusters,
         CGAL::Identity_property_map<Cluster>(),
         m_labels, classifier,
         Neighbor_query(),
         smoothing, subdivisions, indices);

    std::vector<int> ground_truth(m_clusters.size(), -1);
    for (std::size_t i = 0; i < m_clusters.size(); ++ i)
    {
      m_clusters[i].label = indices[i];
      ground_truth[i] = m_clusters[i].training;
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
  Point_set::Property_map<Color> m_color;
  Point_set::Property_map<int> m_cluster_id;
  Point_set::Property_map<int> m_training;
  Point_set::Property_map<int> m_classif;
  
  int m_index_color;

  boost::shared_ptr<Local_eigen_analysis> m_eigen;
  
}; // end class Cluster_classification




#endif // CLUSTER_CLASSIFICATION_H
