#if defined (_MSC_VER) && !defined (_WIN64)
#pragma warning(disable:4244) // boost::number_distance::distance()
                              // converts 64 to 32 bits integers
#endif

#include "Cluster_classification.h"
#include "Color_ramp.h"

#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Timer.h>
#include <CGAL/Memory_sizer.h>

#include <QLineEdit>

#include <CGAL/Three/Viewer_interface.h>

#include <set>
#include <stack>
#include <algorithm>
#include <boost/array.hpp>

Cluster_classification::Cluster_classification(Scene_points_with_normal_item* points)
  : m_points (points)
  , m_input_is_las (false)
{
  m_index_color = 1;

  reset_indices();

  backup_existing_colors_and_add_new();

  bool cluster_found = false;
  boost::tie (m_cluster_id, cluster_found) = m_points->point_set()->property_map<int>("shape");
  if (!cluster_found)
  {
    std::cerr << "Error! Cluster not found!" << std::endl;
    abort();
  }

  CGAL::Classification::create_clusters_from_indices (*(m_points->point_set()),
                                                      m_points->point_set()->point_map(),
                                                      m_cluster_id,
                                                      m_clusters);

  std::cerr << m_clusters.size() << " cluster(s) found" << std::endl;

  bool training_found = false;
  boost::tie (m_training, training_found) = m_points->point_set()->add_property_map<int>("training", -1);
  bool classif_found = false;
  boost::tie (m_classif, classif_found) = m_points->point_set()->add_property_map<int>("label", -1);

  training_found = !training_found; // add_property_map returns false if
  classif_found = !classif_found;   // property was already there

  bool las_found = false;

  if (!classif_found)
  {
    Point_set::Property_map<unsigned char> las_classif;
    boost::tie (las_classif, las_found) = m_points->point_set()->property_map<unsigned char>("classification");
    if (las_found)
    {
      m_input_is_las = true;
      for (Point_set::const_iterator it = m_points->point_set()->begin();
           it != m_points->point_set()->first_selected(); ++ it)
      {
        unsigned char uc = las_classif[*it];
        m_classif[*it] = int(uc);
        if (!training_found)
          m_training[*it] = int(uc);
      }
      m_points->point_set()->remove_property_map (las_classif);
      classif_found = true;
      training_found = true;
    }
  }

  if (training_found || classif_found)
  {
    std::vector<int> used_indices;

    for (Point_set::const_iterator it = m_points->point_set()->begin();
         it != m_points->point_set()->first_selected(); ++ it)
    {
      if (training_found)
      {
        int l = m_training[*it];
        if (l >= 0)
        {
          if (std::size_t(l) >= used_indices.size())
            used_indices.resize(std::size_t(l + 1), -1);
          used_indices[std::size_t(l)] = 0;
        }
      }
      if (classif_found)
      {
        int l = m_classif[*it];
        if (l >= 0)
        {
          if (std::size_t(l) >= used_indices.size())
            used_indices.resize(std::size_t(l + 1), -1);
          used_indices[std::size_t(l)] = 0;
        }
      }
    }

    // map indices to filtered indices
    int current_idx = 0;
    for (std::size_t i = 0; i < used_indices.size(); ++ i)
    {
      if (las_found && (i < 2))
      {
        used_indices[i] = -1;
        continue;
      }
      if (used_indices[i] == -1)
        continue;

      used_indices[i] = current_idx;
      ++ current_idx;
    }


    for (Point_set::const_iterator it = m_points->point_set()->begin();
         it != m_points->point_set()->first_selected(); ++ it)
    {
      int c = m_cluster_id[*it];

      if (training_found)
      {
        if (std::size_t(current_idx) != used_indices.size()) // Empty indices -> reorder indices in point set
        {
          if (las_found && (m_training[*it] == 0 || m_training[*it] == 1)) // Unclassified class in LAS
            m_training[*it] = -1;
          else if (m_training[*it] != -1)
            m_training[*it] = used_indices[std::size_t(m_training[*it])];
        }
        if (c != -1 && m_training[*it] != -1)
          m_clusters[c].training() = m_training[*it];
      }
      if (classif_found)
      {
        if (std::size_t(current_idx) != used_indices.size()) // Empty indices -> reorder indices in point set
        {
          if (las_found && (m_classif[*it] == 0 || m_classif[*it] == 1)) // Unclassified class in LAS
            m_classif[*it] = -1;
          else if (m_classif[*it] != -1)
            m_classif[*it] = used_indices[std::size_t(m_classif[*it])];
        }
        if (c != -1 && m_classif[*it] != -1)
          m_clusters[c].label() = m_classif[*it];
      }
    }

    std::map<int, std::string> label_names;
    if (las_found) // Use LAS standard
    {
      // label_names.insert (std::make_pair (0, std::string("never_clfied")));
      // label_names.insert (std::make_pair (1, std::string("unclassified")));
      label_names.insert (std::make_pair (2, std::string("ground")));
      label_names.insert (std::make_pair (3, std::string("low_veget")));
      label_names.insert (std::make_pair (4, std::string("med_veget")));
      label_names.insert (std::make_pair (5, std::string("high_veget")));
      label_names.insert (std::make_pair (6, std::string("building")));
      label_names.insert (std::make_pair (7, std::string("noise")));
      label_names.insert (std::make_pair (8, std::string("reserved")));
      label_names.insert (std::make_pair (9, std::string("water")));
      label_names.insert (std::make_pair (10, std::string("rail")));
      label_names.insert (std::make_pair (11, std::string("road_surface")));
      label_names.insert (std::make_pair (12, std::string("reserved_2")));
      label_names.insert (std::make_pair (13, std::string("wire_guard")));
      label_names.insert (std::make_pair (14, std::string("wire_conduct")));
      label_names.insert (std::make_pair (15, std::string("trans_tower")));
      label_names.insert (std::make_pair (16, std::string("wire_connect")));
      label_names.insert (std::make_pair (17, std::string("bridge_deck")));
      label_names.insert (std::make_pair (18, std::string("high_noise")));
    }
    else // Try to recover label names from PLY comments
    {

      const std::string& comments = m_points->comments();
      std::istringstream stream (comments);
      std::string line;
      while (getline(stream, line))
      {
        std::stringstream iss (line);
        std::string tag;
        if (iss >> tag && tag == "label")
        {
          int idx;
          std::string name;
          if (iss >> idx >> name)
            label_names.insert (std::make_pair (idx, name));
        }
      }
    }

    for (std::size_t i = 0; i < used_indices.size(); ++ i)
    {
      if (used_indices[i] == -1)
        continue;

      Label_handle new_label;
      std::map<int, std::string>::iterator found
        = label_names.find (int(i));
      if (found != label_names.end())
        new_label = m_labels.add(found->second.c_str());
      else
      {
        std::ostringstream oss;
        oss << "label_" << i;
        new_label = m_labels.add(oss.str().c_str());
      }
    }
  }
  else
  {
    m_labels.add("ground");
    m_labels.add("vegetation");
    m_labels.add("roof");
    m_labels.add("facade");
  }

  update_comments_of_point_set_item();

  m_sowf = new Sum_of_weighted_features (m_labels, m_features);
  m_ethz = NULL;
#ifdef CGAL_LINKED_WITH_OPENCV
  m_random_forest = NULL;
#endif
#ifdef CGAL_LINKED_WITH_TENSORFLOW
  m_neural_network = NULL;
#endif

  // Compute neighborhood
#if 0
  typedef CGAL::Triangulation_vertex_base_with_info_3<int, Kernel>         Vb;
  typedef CGAL::Delaunay_triangulation_cell_base_3<Kernel>                 Cb;
  typedef CGAL::Triangulation_data_structure_3<Vb, Cb>                     Tds;
  typedef CGAL::Delaunay_triangulation_3<Kernel, Tds>                      Delaunay;

  Delaunay dt (boost::make_transform_iterator
               (m_points->point_set()->begin(),
                Point_set_with_cluster_info (m_points->point_set(),
                                             m_cluster_id)),
                boost::make_transform_iterator
               (m_points->point_set()->end(),
                Point_set_with_cluster_info (m_points->point_set(),
                                             m_cluster_id)));

  std::set<std::pair<int, int> > adjacencies;

  for (Delaunay::Finite_edges_iterator it = dt.finite_edges_begin();
       it != dt.finite_edges_end(); ++ it)
  {
    Delaunay::Vertex_handle v0 = it->first->vertex (it->second);
    Delaunay::Vertex_handle v1 = it->first->vertex (it->third);
    int a = v0->info();
    int b = v1->info();

    if (a == -1 || b == -1)
      continue;

    if (a > b)
      std::swap (a, b);

    adjacencies.insert (std::make_pair (a, b));
  }

  for (std::set<std::pair<int, int> >::iterator it = adjacencies.begin();
       it != adjacencies.end(); ++ it)
  {
    m_clusters[std::size_t(it->first)].neighbors->push_back (std::size_t(it->second));
    m_clusters[std::size_t(it->second)].neighbors->push_back (std::size_t(it->first));
  }
#endif
}


Cluster_classification::~Cluster_classification()
{
  if (m_sowf != NULL)
    delete m_sowf;
  if (m_ethz != NULL)
    delete m_ethz;
#ifdef CGAL_LINKED_WITH_OPENCV
  if (m_random_forest != NULL)
    delete m_random_forest;
#endif
#ifdef CGAL_LINKED_WITH_TENSORFLOW
  if (m_neural_network != NULL)
    delete m_neural_network;
#endif
  if (m_points != NULL)
  {
    for (Point_set::const_iterator it = m_points->point_set()->begin();
         it != m_points->point_set()->first_selected(); ++ it)
    {
      int c = m_cluster_id[*it];
      if (c != -1)
      {
        m_training[*it] = m_clusters[c].training();
        m_classif[*it] = m_clusters[c].label();
      }
      else
      {
        m_training[*it] = -1;
        m_classif[*it] = -1;
      }
    }

    // For LAS saving, convert classification info in the LAS standard
//    if (m_input_is_las)
    {
      Point_set::Property_map<unsigned char> las_classif
        = m_points->point_set()->add_property_map<unsigned char>("classification", 0).first;

      std::vector<unsigned char> label_indices;

      unsigned char custom = 19;
      for (std::size_t i = 0; i < m_labels.size(); ++ i)
      {
        if (m_labels[i]->name() == "ground")
          label_indices.push_back (2);
        else if (m_labels[i]->name() == "low_veget")
          label_indices.push_back (3);
        else if (m_labels[i]->name() == "med_veget" || m_labels[i]->name() == "vegetation")
          label_indices.push_back (4);
        else if (m_labels[i]->name() == "high_veget")
          label_indices.push_back (5);
        else if (m_labels[i]->name() == "building" || m_labels[i]->name() == "roof")
          label_indices.push_back (6);
        else if (m_labels[i]->name() == "noise")
          label_indices.push_back (7);
        else if (m_labels[i]->name() == "reserved" || m_labels[i]->name() == "facade")
          label_indices.push_back (8);
        else if (m_labels[i]->name() == "water")
          label_indices.push_back (9);
        else if (m_labels[i]->name() == "rail")
          label_indices.push_back (10);
        else if (m_labels[i]->name() == "road_surface")
          label_indices.push_back (11);
        else if (m_labels[i]->name() == "reserved_2")
          label_indices.push_back (12);
        else if (m_labels[i]->name() == "wire_guard")
          label_indices.push_back (13);
        else if (m_labels[i]->name() == "wire_conduct")
          label_indices.push_back (14);
        else if (m_labels[i]->name() == "trans_tower")
          label_indices.push_back (15);
        else if (m_labels[i]->name() == "wire_connect")
          label_indices.push_back (16);
        else if (m_labels[i]->name() == "bridge_deck")
          label_indices.push_back (17);
        else if (m_labels[i]->name() == "high_noise")
          label_indices.push_back (18);
        else
          label_indices.push_back (custom ++);
      }

      for (Point_set::const_iterator it = m_points->point_set()->begin();
           it != m_points->point_set()->end(); ++ it)
      {
        int c = m_classif[*it];
        unsigned char lc = 1; // unclassified in LAS standard
        if (c != -1)
          lc = label_indices[std::size_t(c)];

        las_classif[*it] = lc;

        int t = m_training[*it];
        unsigned char lt = 1; // unclassified in LAS standard
        if (t != -1)
          lt = label_indices[std::size_t(t)];

        m_training[*it] = int(lt);
      }

      m_points->point_set()->remove_property_map (m_classif);
    }


    reset_colors();
    erase_item();
  }
}


void Cluster_classification::backup_existing_colors_and_add_new()
{
  if (m_points->point_set()->has_colors())
  {
    m_color = m_points->point_set()->add_property_map<CGAL::Color>("real_color").first;
    for (Point_set::const_iterator it = m_points->point_set()->begin();
         it != m_points->point_set()->first_selected(); ++ it)
      m_color[*it] = CGAL::Color ((unsigned char)(255 * m_points->point_set()->red(*it)),
                                  (unsigned char)(255 * m_points->point_set()->green(*it)),
                                  (unsigned char)(255 * m_points->point_set()->blue(*it)));

    m_points->point_set()->remove_colors();
  }

  m_points->point_set()->add_colors();
}

void Cluster_classification::reset_colors()
{
  if (m_color == Point_set::Property_map<CGAL::Color>())
    m_points->point_set()->remove_colors();
  else
  {
    for (Point_set::const_iterator it = m_points->point_set()->begin();
         it != m_points->point_set()->first_selected(); ++ it)
      m_points->point_set()->set_color(*it, m_color[*it]);

    m_points->point_set()->remove_property_map(m_color);
  }
}

void Cluster_classification::change_color (int index, float* vmin, float* vmax)
{
  m_index_color = index;

  int index_color = real_index_color();

  // Colors
  static Color_ramp ramp;
  ramp.build_rainbow();
  reset_indices();

  if (index_color == -1) // item color
    m_points->point_set()->remove_colors();
  else
  {
    if (!m_points->point_set()->has_colors())
      m_points->point_set()->add_colors();

    if (index_color == 0) // real colors
    {

      for (Point_set::const_iterator it = m_points->point_set()->begin();
           it != m_points->point_set()->first_selected(); ++ it)
        m_points->point_set()->set_color(*it, m_color[*it]);
    }
    else if (index_color == 1) // classif
    {
      for (Point_set::const_iterator it = m_points->point_set()->begin();
           it != m_points->point_set()->first_selected(); ++ it)
      {
        QColor color (0, 0, 0);
        int cid = m_cluster_id[*it];
        if (cid != -1)
        {
          std::size_t c = m_clusters[cid].label();

          if (c != std::size_t(-1))
            color = label_qcolor (m_labels[std::size_t(c)]);
        }

        m_points->point_set()->set_color(*it, color);
      }
    }
    else if (index_color == 2) // training
    {
      for (Point_set::const_iterator it = m_points->point_set()->begin();
           it != m_points->point_set()->first_selected(); ++ it)
      {
        QColor color (0, 0, 0);
        int cid = m_cluster_id[*it];
        float div = 1;

        if (cid != -1)
        {
          int c = m_clusters[cid].training();
          int c2 = m_clusters[cid].label();

          if (c != -1)
            color = label_qcolor (m_labels[std::size_t(c)]);

          if (c != c2)
            div = 2;
        }
        m_points->point_set()->set_color(*it, color.red() / div, color.green() / div, color.blue() / div);
      }
    }
    else if (index_color == 3) // clusters
    {
      for (Point_set::const_iterator it = m_points->point_set()->begin();
           it != m_points->point_set()->first_selected(); ++ it)
      {
        int cid = m_cluster_id[*it];

        if (cid != -1)
        {
          srand(cid);
          m_points->point_set()->set_color(*it, 64 + rand() % 192, 64 + rand() % 192, 64 + rand() % 192);
        }
        else
        {
          m_points->point_set()->set_color(*it);
        }
      }
    }
    else
    {
      std::size_t corrected_index = index_color - 4;
      if (corrected_index < m_labels.size()) // Display label probabilities
      {
        if (m_label_probabilities.size() <= corrected_index ||
            m_label_probabilities[corrected_index].size() != m_clusters.size())
        {
          for (Point_set::const_iterator it = m_points->point_set()->begin();
               it != m_points->point_set()->first_selected(); ++ it)
            m_points->point_set()->set_color(*it);
        }
        else
        {
          for (Point_set::const_iterator it = m_points->point_set()->begin();
               it != m_points->point_set()->first_selected(); ++ it)
          {
            int cid = m_cluster_id[*it];
            if (cid != -1)
            {
              float v = (std::max) (0.f, (std::min)(1.f, m_label_probabilities[corrected_index][cid]));
              m_points->point_set()->set_color(*it, ramp.r(v) * 255, ramp.g(v) * 255, ramp.b(v) * 255);
            }
            else
              m_points->point_set()->set_color(*it);
          }
        }
      }
      else
      {
        corrected_index -= m_labels.size();
        if (corrected_index >= m_features.size())
        {
          std::cerr << "Error: trying to access feature " << corrected_index << " out of " << m_features.size() << std::endl;
          return;
        }

        Feature_handle feature = m_features[corrected_index];

        float min = (std::numeric_limits<float>::max)();
        float max = -(std::numeric_limits<float>::max)();

        if (vmin != NULL && vmax != NULL
            && *vmin != std::numeric_limits<float>::infinity()
            && *vmax != std::numeric_limits<float>::infinity())
        {
          min = *vmin;
          max = *vmax;
        }
        else
        {
          for (Point_set::const_iterator it = m_points->point_set()->begin();
               it != m_points->point_set()->first_selected(); ++ it)
          {
            int cid = m_cluster_id[*it];
            if (cid != -1)
            {
              if (feature->value(cid) > max)
                max = feature->value(cid);
              if (feature->value(cid) < min)
                min = feature->value(cid);
            }
          }
        }

        for (Point_set::const_iterator it = m_points->point_set()->begin();
             it != m_points->point_set()->first_selected(); ++ it)
        {
          int cid = m_cluster_id[*it];
          if (cid != -1)
          {
            float v = (feature->value(cid) - min) / (max - min);
            if (v < 0.f) v = 0.f;
            if (v > 1.f) v = 1.f;

            m_points->point_set()->set_color(*it, ramp.r(v) * 255, ramp.g(v) * 255, ramp.b(v) * 255);
          }
          else
            m_points->point_set()->set_color(*it);
        }

        if (vmin != NULL && vmax != NULL)
        {
          *vmin = min;
          *vmax = max;
        }
      }
    }
  }

  for (Point_set::const_iterator it = m_points->point_set()->first_selected();
       it != m_points->point_set()->end(); ++ it)
    m_points->point_set()->set_color(*it, 255, 0, 0);
}

int Cluster_classification::real_index_color() const
{
  int out = m_index_color;

  if (out == 0 && m_color == Point_set::Property_map<CGAL::Color>())
    out = -1;
  return out;
}

void Cluster_classification::reset_indices ()
{
  Point_set::Property_map<Point_set::Index> indices
    = m_points->point_set()->property_map<Point_set::Index>("index").first;

  m_points->point_set()->unselect_all();
  Point_set::Index idx;
  ++ idx;
  for (std::size_t i = 0; i < m_points->point_set()->size(); ++ i)
    *(indices.begin() + i) = idx ++;
}

void Cluster_classification::compute_features (std::size_t nb_scales, float voxel_size)
{
  CGAL_assertion (!(m_points->point_set()->empty()));

  reset_indices();

  std::cerr << "Computing pointwise features with " << nb_scales << " scale(s) and ";
  if (voxel_size == -1)
    std::cerr << "automatic voxel size" << std::endl;
  else
    std::cerr << "voxel size = " << voxel_size << std::endl;

  m_features.clear();

  Point_set::Vector_map normal_map;
  bool normals = m_points->point_set()->has_normal_map();
  if (normals)
    normal_map = m_points->point_set()->normal_map();

  bool colors = (m_color != Point_set::Property_map<CGAL::Color>());

  Point_set::Property_map<boost::uint8_t> echo_map;
  bool echo;
  boost::tie (echo_map, echo) = m_points->point_set()->template property_map<boost::uint8_t>("echo");
  if (!echo)
    boost::tie (echo_map, echo) = m_points->point_set()->template property_map<boost::uint8_t>("number_of_returns");

  Feature_set pointwise_features;

  Generator generator (*(m_points->point_set()), m_points->point_set()->point_map(), nb_scales, voxel_size);

  CGAL::Real_timer t;
  t.start();

#ifdef CGAL_LINKED_WITH_TBB
  pointwise_features.begin_parallel_additions();
#endif

  generator.generate_point_based_features(pointwise_features);
  if (normals)
    generator.generate_normal_based_features (pointwise_features, normal_map);
  if (colors)
    generator.generate_color_based_features (pointwise_features, m_color);
  if (echo)
    generator.generate_echo_based_features (pointwise_features, echo_map);

#ifdef CGAL_LINKED_WITH_TBB
  pointwise_features.end_parallel_additions();
#endif

  add_remaining_point_set_properties_as_features(pointwise_features);

  t.stop();
  std::cerr << pointwise_features.size() << " feature(s) computed in " << t.time() << " second(s)" << std::endl;
  t.reset();

  std::cerr << "Computing cluster features" << std::endl;
  t.start();

#ifdef CGAL_LINKED_WITH_TBB
  m_features.begin_parallel_additions();
#endif

  for (std::size_t i = 0; i < pointwise_features.size(); ++ i)
  {
    m_features.template add<Mean_of_feature> (m_clusters,
                                              pointwise_features[i]);
  }

#ifdef CGAL_LINKED_WITH_TBB
  m_features.end_parallel_additions();
  m_features.begin_parallel_additions();
#endif

  for (std::size_t i = 0; i < pointwise_features.size(); ++ i)
  {
    m_features.template add<Variance_of_feature> (m_clusters,
                                                  pointwise_features[i],
                                                  m_features[i]);
  }

  add_cluster_features();

#ifdef CGAL_LINKED_WITH_TBB
  m_features.end_parallel_additions();
#endif

  t.stop();
  std::cerr << m_features.size() << " feature(s) computed in " << t.time() << " second(s)" << std::endl;


  delete m_sowf;
  m_sowf = new Sum_of_weighted_features (m_labels, m_features);
  if (m_ethz != NULL)
  {
    delete m_ethz;
    m_ethz = NULL;
  }
#ifdef CGAL_LINKED_WITH_OPENCV
  if (m_random_forest != NULL)
  {
    delete m_random_forest;
    m_random_forest = NULL;
  }
#endif
#ifdef CGAL_LINKED_WITH_TENSORFLOW
  if (m_neural_network != NULL)
  {
    delete m_neural_network;
    m_neural_network = NULL;
  }
#endif

  std::cerr << "Features = " << m_features.size() << std::endl;
}

void Cluster_classification::select_random_region()
{
  std::size_t c = rand() % m_clusters.size();
  std::set<Point_set::Index> in_cluster;
  for (std::size_t i = 0; i < m_clusters[c].size(); ++ i)
    in_cluster.insert (m_clusters[c].index(i));

  std::vector<std::size_t> selected;
  std::vector<std::size_t> unselected;

  for (Point_set::const_iterator it = m_points->point_set()->begin();
       it != m_points->point_set()->end(); ++ it)
    if (in_cluster.find (*it) != in_cluster.end())
      selected.push_back (*it);
    else
      unselected.push_back (*it);

  for (std::size_t i = 0; i < unselected.size(); ++ i)
    *(m_points->point_set()->begin() + i) = unselected[i];
  for (std::size_t i = 0; i < selected.size(); ++ i)
    *(m_points->point_set()->begin() + (unselected.size() + i)) = selected[i];

  m_points->point_set()->set_first_selected
    (m_points->point_set()->begin() + unselected.size());

}

void Cluster_classification::add_remaining_point_set_properties_as_features(Feature_set& feature_set)
{
  const std::vector<std::string>& prop = m_points->point_set()->base().properties();

  for (std::size_t i = 0; i < prop.size(); ++ i)
  {
    if (prop[i] == "index" ||
        prop[i] == "point" ||
        prop[i] == "normal" ||
        prop[i] == "echo" ||
        prop[i] == "number_of_returns" ||
        prop[i] == "training" ||
        prop[i] == "label" ||
        prop[i] == "classification" ||
        prop[i] == "scan_direction_flag" ||
        prop[i] == "real_color" ||
        prop[i] == "shape" ||
        prop[i] == "red" || prop[i] == "green" || prop[i] == "blue" ||
        prop[i] == "r" || prop[i] == "g" || prop[i] == "b")
      continue;

    if (try_adding_simple_feature<boost::int8_t>(feature_set, prop[i]))
      continue;
    if (try_adding_simple_feature<boost::uint8_t>(feature_set, prop[i]))
      continue;
    if (try_adding_simple_feature<boost::int16_t>(feature_set, prop[i]))
      continue;
    if (try_adding_simple_feature<boost::uint16_t>(feature_set, prop[i]))
      continue;
    if (try_adding_simple_feature<boost::int32_t>(feature_set, prop[i]))
      continue;
    if (try_adding_simple_feature<boost::uint32_t>(feature_set, prop[i]))
      continue;
    if (try_adding_simple_feature<float>(feature_set, prop[i]))
      continue;
    if (try_adding_simple_feature<double>(feature_set, prop[i]))
      continue;
  }
}

void Cluster_classification::train(int classifier, const QMultipleInputDialog& dialog)
{
  if (m_features.size() == 0)
  {
    std::cerr << "Error: features not computed" << std::endl;
    return;
  }
  reset_indices();

  m_label_probabilities.clear();
  m_label_probabilities.resize (m_labels.size());
  for (std::size_t i = 0; i < m_label_probabilities.size(); ++ i)
    m_label_probabilities[i].resize (m_clusters.size(), -1);

  std::vector<std::size_t> nb_label (m_labels.size(), 0);
  std::size_t nb_total = 0;

  std::vector<int> training;
  training.reserve (m_clusters.size());
  for (std::size_t i = 0; i < m_clusters.size(); ++ i)
  {
    training.push_back (m_clusters[i].training());
    if (training.back() != -1)
    {
      nb_label[std::size_t(training.back())] ++;
      ++ nb_total;
    }
  }

  std::cerr << nb_total << " cluster(s) used for training ("
            << 100. * (nb_total / double(m_clusters.size())) << "% of the total):" << std::endl;
  for (std::size_t i = 0; i < m_labels.size(); ++ i)
    std::cerr << " * " << m_labels[i]->name() << ": " << nb_label[i] << " clusters(s)" << std::endl;

  std::vector<int> indices (m_clusters.size(), -1);

  if (classifier == CGAL_CLASSIFICATION_SOWF_NUMBER)
  {
    m_sowf->train<Concurrency_tag>(training, dialog.get<QSpinBox>("trials")->value());
    CGAL::Classification::classify<Concurrency_tag> (m_clusters,
                                                     m_labels, *m_sowf,
                                                     indices, m_label_probabilities);
  }
  else if (classifier == CGAL_CLASSIFICATION_ETHZ_NUMBER)
  {
    if (m_ethz != NULL)
      delete m_ethz;
    m_ethz = new ETHZ_random_forest (m_labels, m_features);
    m_ethz->train<Concurrency_tag>(training, true,
                                   dialog.get<QSpinBox>("num_trees")->value(),
                                   dialog.get<QSpinBox>("max_depth")->value());
    CGAL::Classification::classify<Concurrency_tag> (m_clusters,
                                                     m_labels, *m_ethz,
                                                     indices, m_label_probabilities);
  }
  else if (classifier == CGAL_CLASSIFICATION_OPENCV_NUMBER)
  {
#ifdef CGAL_LINKED_WITH_OPENCV
    if (m_random_forest != NULL)
      delete m_random_forest;
    m_random_forest = new Random_forest (m_labels, m_features,
                                         dialog.get<QSpinBox>("max_depth")->value(), 5, 15,
                                         dialog.get<QSpinBox>("num_trees")->value());
    m_random_forest->train (training);
    CGAL::Classification::classify<Concurrency_tag> (m_clusters,
                                                     m_labels, *m_random_forest,
                                                     indices, m_label_probabilities);
#endif
  }
  else if (classifier == CGAL_CLASSIFICATION_TENSORFLOW_NUMBER)
  {
#ifdef CGAL_LINKED_WITH_TENSORFLOW
    if (m_neural_network != NULL)
    {
      if (m_neural_network->initialized())
      {
        if (dialog.get<QCheckBox>("restart")->isChecked())
        {
          delete m_neural_network;
          m_neural_network = new Neural_network (m_labels, m_features);
        }
      }
      else
      {
        delete m_neural_network;
        m_neural_network = new Neural_network (m_labels, m_features);
      }
    }
    else
      m_neural_network = new Neural_network (m_labels, m_features);

    std::vector<std::size_t> hidden_layers;

    std::string hl_input = dialog.get<QLineEdit>("hidden_layers")->text().toStdString();
    if (hl_input != "")
    {
      std::istringstream iss(hl_input);
      int s;
      while (iss >> s)
        hidden_layers.push_back (std::size_t(s));
    }

    m_neural_network->train (training,
                             dialog.get<QCheckBox>("restart")->isChecked(),
                             dialog.get<QSpinBox>("trials")->value(),
                             dialog.get<DoubleEdit>("learning_rate")->value(),
                             dialog.get<QSpinBox>("batch_size")->value(),
                             hidden_layers);

    CGAL::Classification::classify<Concurrency_tag> (m_clusters,
                                                     m_labels, *m_neural_network,
                                                     indices, m_label_probabilities);
#endif
  }

  for (std::size_t i = 0; i < m_clusters.size(); ++ i)
    m_clusters[i].label() = indices[i];

  if (m_index_color == 1 || m_index_color == 2)
    change_color (m_index_color);
}

bool Cluster_classification::run (int method, int classifier,
                                  std::size_t subdivisions,
                                  double smoothing)
{
  if (m_features.size() == 0)
  {
    std::cerr << "Error: features not computed" << std::endl;
    return false;
  }
  reset_indices();

  if (classifier == CGAL_CLASSIFICATION_SOWF_NUMBER)
    run (method, *m_sowf, subdivisions, smoothing);
  else if (classifier == CGAL_CLASSIFICATION_ETHZ_NUMBER)
  {
    if (m_ethz == NULL)
    {
      std::cerr << "Error: ETHZ Random Forest must be trained or have a configuration loaded first" << std::endl;
      return false;
    }
    run (method, *m_ethz, subdivisions, smoothing);
  }
  else if (classifier == CGAL_CLASSIFICATION_OPENCV_NUMBER)
  {
#ifdef CGAL_LINKED_WITH_OPENCV
    if (m_random_forest == NULL)
    {
      std::cerr << "Error: OpenCV Random Forest must be trained or have a configuration loaded first" << std::endl;
      return false;
    }
    run (method, *m_random_forest, subdivisions, smoothing);
#endif
  }
  else if (classifier == CGAL_CLASSIFICATION_TENSORFLOW_NUMBER)
  {
#ifdef CGAL_LINKED_WITH_TENSORFLOW
    if (m_neural_network == NULL)
    {
      std::cerr << "Error: TensorFlow Neural Network must be trained or have a configuration loaded first" << std::endl;
      return false;
    }
    run (method, *m_neural_network, subdivisions, smoothing);
#endif
  }

  return true;
}
