#if defined (_MSC_VER) && !defined (_WIN64)
#pragma warning(disable:4244) // boost::number_distance::distance()
                              // converts 64 to 32 bits integers
#endif

#include "Point_set_item_classification.h"
#include "Color_ramp.h"

#include <CGAL/Timer.h>
#include <CGAL/Memory_sizer.h>

#include <CGAL/Three/Viewer_interface.h>

#include <set>
#include <stack>
#include <algorithm>
#include <boost/array.hpp>

Point_set_item_classification::Point_set_item_classification(Scene_points_with_normal_item* points)
  : m_points (points),
    m_generator (NULL)
{
  m_index_color = 1;

  reset_indices();

  Point_set::Property_map<int> cluster_id;
  bool cluster_found = false;
  boost::tie (cluster_id, cluster_found) = m_points->point_set()->property_map<int>("shape");
  if (cluster_found)
  {
    for (Point_set::const_iterator it = m_points->point_set()->begin();
         it != m_points->point_set()->first_selected(); ++ it)
    {
      int c = cluster_id[*it];
      if (c == -1)
        continue;
      if (std::size_t(c) >= m_clusters.size())
        m_clusters.resize (c + 1);
      m_clusters[c].inliers.push_back (*it);
    }
  }

  backup_existing_colors_and_add_new();
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

    if (std::size_t(current_idx) != used_indices.size()) // Empty indices -> reorder indices in point set
    {
      for (Point_set::const_iterator it = m_points->point_set()->begin();
           it != m_points->point_set()->first_selected(); ++ it)
      {
        if (training_found)
        {
          if (las_found && (m_training[*it] == 0 || m_training[*it] == 1)) // Unclassified class in LAS
            m_training[*it] = -1;
          else if (m_training[*it] != -1)
            m_training[*it] = used_indices[std::size_t(m_training[*it])];
        }
        if (classif_found)
        {
          if (las_found && (m_classif[*it] == 0 || m_classif[*it] == 1)) // Unclassified class in LAS
            m_classif[*it] = -1;
          else if (m_classif[*it] != -1)
            m_classif[*it] = used_indices[std::size_t(m_classif[*it])];
        }
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

      m_label_colors.push_back (this->get_new_label_color (new_label->name()));
    }
  }
  else
  {
    m_labels.add("ground");
    m_labels.add("vegetation");
    m_labels.add("roof");
    m_labels.add("facade");

    for (std::size_t i = 0; i < m_labels.size(); ++ i)
      m_label_colors.push_back (this->get_new_label_color (m_labels[i]->name()));
  }
  
  update_comments_of_point_set_item();

  m_sowf = new Sum_of_weighted_features (m_labels, m_features);
  m_ethz = new ETHZ_random_forest (m_labels, m_features);
#ifdef CGAL_LINKED_WITH_OPENCV
  m_random_forest = new Random_forest (m_labels, m_features);
#endif
}


Point_set_item_classification::~Point_set_item_classification()
{
  if (m_sowf != NULL)
    delete m_sowf;
  if (m_ethz != NULL)
    delete m_ethz;
#ifdef CGAL_LINKED_WITH_OPENCV
  if (m_random_forest != NULL)
    delete m_random_forest;
#endif
  if (m_generator != NULL)
    delete m_generator;
  if (m_points != NULL)
    {
      reset_colors();
      erase_item();
    }
}


void Point_set_item_classification::backup_existing_colors_and_add_new()
{
  if (m_points->point_set()->has_colors())
    {
      m_color = m_points->point_set()->add_property_map<Color>("real_color").first;
      for (Point_set::const_iterator it = m_points->point_set()->begin();
           it != m_points->point_set()->first_selected(); ++ it)
        m_color[*it] = {{ (unsigned char)(255 * m_points->point_set()->red(*it)),
                          (unsigned char)(255 * m_points->point_set()->green(*it)),
                          (unsigned char)(255 * m_points->point_set()->blue(*it)) }};

      m_points->point_set()->remove_colors();
    }
      
  m_red = m_points->point_set()->add_property_map<unsigned char>("red").first;
  m_green = m_points->point_set()->add_property_map<unsigned char>("green").first;
  m_blue = m_points->point_set()->add_property_map<unsigned char>("blue").first;
  for (Point_set::const_iterator it = m_points->point_set()->begin();
       it != m_points->point_set()->first_selected(); ++ it)
    {
      m_red[*it] = 0;
      m_green[*it] = 0;
      m_blue[*it] = 0;
    }
  m_points->point_set()->check_colors();
}

void Point_set_item_classification::reset_colors()
{
  if (m_color == Point_set::Property_map<Color>())
    m_points->point_set()->remove_colors();
  else
    {
      for (Point_set::const_iterator it = m_points->point_set()->begin();
           it != m_points->point_set()->first_selected(); ++ it)
        {
          m_red[*it] = m_color[*it][0];
          m_green[*it] = m_color[*it][1];
          m_blue[*it] = m_color[*it][2];
        }
      m_points->point_set()->remove_property_map(m_color);
    }
}

// Write point set to .PLY file
bool Point_set_item_classification::write_output(std::ostream& stream)
{
  if (m_features.size() == 0)
    return false;

  reset_indices();
  
  stream.precision (std::numeric_limits<double>::digits10 + 2);

  // std::vector<Color> colors;
  // for (std::size_t i = 0; i < m_labels.size(); ++ i)
  //   {
  //     Color c = {{ (unsigned char)(m_labels[i].second.red()),
  //                  (unsigned char)(m_labels[i].second.green()),
  //                  (unsigned char)(m_labels[i].second.blue()) }};
  //     colors.push_back (c);
  //   }
  
  //  m_psc->write_classification_to_ply (stream);
  return true;
}


void Point_set_item_classification::change_color (int index)
{
  m_index_color = index;

  int index_color = real_index_color();
    
  // Colors
  static Color_ramp ramp;
  ramp.build_rainbow();
  reset_indices();
  if (index_color == -1) // item color
    {
      for (Point_set::const_iterator it = m_points->point_set()->begin();
           it != m_points->point_set()->first_selected(); ++ it)
        {
          m_red[*it] = 0;
          m_green[*it] = 0;
          m_blue[*it] = 0;
        }
    }
  else if (index_color == 0) // real colors
    {

      for (Point_set::const_iterator it = m_points->point_set()->begin();
           it != m_points->point_set()->first_selected(); ++ it)
        {
          m_red[*it] = m_color[*it][0];
          m_green[*it] = m_color[*it][1];
          m_blue[*it] = m_color[*it][2];
        }
    }
  else if (index_color == 1) // classif
    {
      for (Point_set::const_iterator it = m_points->point_set()->begin();
           it != m_points->point_set()->first_selected(); ++ it)
        {
          QColor color (0, 0, 0);
          std::size_t c = m_classif[*it];
          
          if (c != std::size_t(-1))
            color = m_label_colors[c];

          m_red[*it] = color.red();
          m_green[*it] = color.green();
          m_blue[*it] = color.blue();
        }
    }
  else if (index_color == 2) // training
    {
      for (Point_set::const_iterator it = m_points->point_set()->begin();
           it != m_points->point_set()->first_selected(); ++ it)
        {
          QColor color (0, 0, 0);
          int c = m_training[*it];
          int c2 = m_classif[*it];
          
          if (c != -1)
            color = m_label_colors[std::size_t(c)];
          
          float div = 1;
          if (c != c2)
            div = 2;
          
          m_red[*it] = (color.red() / div);
          m_green[*it] = (color.green() / div);
          m_blue[*it] = (color.blue() / div);
        }
    }
  else
    {
      Feature_handle feature = m_features[index_color - 3];

      float max = 0.;
      for (Point_set::const_iterator it = m_points->point_set()->begin();
           it != m_points->point_set()->first_selected(); ++ it)
        if (feature->value(*it) > max)
          max = feature->value(*it);

      for (Point_set::const_iterator it = m_points->point_set()->begin();
           it != m_points->point_set()->first_selected(); ++ it)
        {
          float v = std::max (0.f, feature->value(*it) / max);
          m_red[*it] = (unsigned char)(ramp.r(v) * 255);
          m_green[*it] = (unsigned char)(ramp.g(v) * 255);
          m_blue[*it] = (unsigned char)(ramp.b(v) * 255);
        }
    }

  for (Point_set::const_iterator it = m_points->point_set()->first_selected();
       it != m_points->point_set()->end(); ++ it)
    {
      m_red[*it] = 255;
      m_green[*it] = 0;
      m_blue[*it] = 0;
    }

}

int Point_set_item_classification::real_index_color() const
{
  int out = m_index_color;
  
  if (out == 0 && m_color == Point_set::Property_map<Color>())
    out = -1;
  return out;
}

void Point_set_item_classification::reset_indices ()
{
  Point_set::Property_map<Point_set::Index> indices
    = m_points->point_set()->property_map<Point_set::Index>("index").first;

  m_points->point_set()->unselect_all();
  Point_set::Index idx;
  ++ idx;
  for (std::size_t i = 0; i < m_points->point_set()->size(); ++ i)
    *(indices.begin() + i) = idx ++;
}

void Point_set_item_classification::compute_features (std::size_t nb_scales)
{
  CGAL_assertion (!(m_points->point_set()->empty()));

  if (m_generator != NULL)
    delete m_generator;

  reset_indices();
  
  std::cerr << "Computing features with " << nb_scales << " scale(s)" << std::endl;
  m_features.clear();

  bool normals = m_points->point_set()->has_normal_map();
  bool colors = (m_color != Point_set::Property_map<Color>());
  Point_set::Property_map<boost::uint8_t> echo_map;
  bool echo;
  boost::tie (echo_map, echo) = m_points->point_set()->template property_map<boost::uint8_t>("echo");
  if (!echo)
    boost::tie (echo_map, echo) = m_points->point_set()->template property_map<boost::uint8_t>("number_of_returns");

  m_generator = new Generator (*(m_points->point_set()), m_points->point_set()->point_map(), nb_scales);

  CGAL::Real_timer t;
  t.start();
    
#ifdef CGAL_LINKED_WITH_TBB
  m_features.begin_parallel_additions();
#endif

  m_generator->generate_point_based_features(m_features);
  if (normals)
    m_generator->generate_normal_based_features (m_features, m_points->point_set()->normal_map());
  if (colors)
    m_generator->generate_color_based_features (m_features, m_color);
  if (echo)
    m_generator->generate_echo_based_features (m_features, echo_map);
  
  add_remaining_point_set_properties_as_features();

#ifdef CGAL_LINKED_WITH_TBB
  m_features.end_parallel_additions();
#endif

  delete m_sowf;
  m_sowf = new Sum_of_weighted_features (m_labels, m_features);
  delete m_ethz;
  m_ethz = new ETHZ_random_forest (m_labels, m_features);

#ifdef CGAL_LINKED_WITH_OPENCV
  delete m_random_forest;
  m_random_forest = new Random_forest (m_labels, m_features);
#endif

  t.stop();
  std::cerr << m_features.size() << " feature(s) computed in " << t.time() << " second(s)" << std::endl;
}

void Point_set_item_classification::select_random_region()
{
  m_points->point_set()->reset_indices();

  std::vector<std::size_t> selected;
  std::vector<std::size_t> unselected;

  // If no cluster, fallback mode on grid or region selection
  if (m_clusters.empty())
  {
    std::size_t scale = (rand() % m_generator->number_of_scales());

    bool use_grid = (rand() % 2);

    if (use_grid)
    {
      std::size_t x = (rand() % m_generator->grid(scale).width());
      std::size_t y = (rand() % m_generator->grid(scale).height());
      std::copy (m_generator->grid(scale).indices_begin(x,y),
                 m_generator->grid(scale).indices_end(x,y),
                 std::back_inserter (selected));
    }
    else
    {
      m_generator->neighborhood(0).sphere_neighbor_query (m_generator->radius_neighbors(scale))
        (*(m_points->point_set()->points().begin() + (rand() % m_points->point_set()->size())),
         std::back_inserter (selected));
    }

    if (selected.empty())
      return;

    std::sort (selected.begin(), selected.end());
    std::size_t current_idx = 0;
    for (Point_set::const_iterator it = m_points->point_set()->begin();
         it != m_points->point_set()->end(); ++ it)
      if (std::size_t(*it) == selected[current_idx])
        current_idx ++;
      else
        unselected.push_back (*it);
  
  }
  else
  {
    std::size_t c = rand() % m_clusters.size();
    std::set<Point_set::Index> in_cluster;
    for (std::size_t i = 0; i < m_clusters[c].size(); ++ i)
      in_cluster.insert (m_clusters[c][i]);

    for (Point_set::const_iterator it = m_points->point_set()->begin();
         it != m_points->point_set()->end(); ++ it)
      if (in_cluster.find (*it) != in_cluster.end())
        selected.push_back (*it);
      else
        unselected.push_back (*it);
  }
  
  for (std::size_t i = 0; i < unselected.size(); ++ i)
    *(m_points->point_set()->begin() + i) = unselected[i];
  for (std::size_t i = 0; i < selected.size(); ++ i)
    *(m_points->point_set()->begin() + (unselected.size() + i)) = selected[i];

  m_points->point_set()->set_first_selected
    (m_points->point_set()->begin() + unselected.size());

}

void Point_set_item_classification::add_remaining_point_set_properties_as_features()
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
        prop[i] == "real_color" ||
        prop[i] == "shape" ||
        prop[i] == "red" || prop[i] == "green" || prop[i] == "blue" ||
        prop[i] == "r" || prop[i] == "g" || prop[i] == "b")
      continue;

    if (try_adding_simple_feature<boost::int8_t>(prop[i]))
      continue;
    if (try_adding_simple_feature<boost::uint8_t>(prop[i]))
      continue;
    if (try_adding_simple_feature<boost::int16_t>(prop[i]))
      continue;
    if (try_adding_simple_feature<boost::uint16_t>(prop[i]))
      continue;
    if (try_adding_simple_feature<boost::int32_t>(prop[i]))
      continue;
    if (try_adding_simple_feature<boost::uint32_t>(prop[i]))
      continue;
    if (try_adding_simple_feature<float>(prop[i]))
      continue;
    if (try_adding_simple_feature<double>(prop[i]))
      continue;
  }
}

void Point_set_item_classification::train(int classifier, unsigned int nb_trials,
                                          std::size_t num_trees, std::size_t max_depth)
{
  if (m_features.size() == 0)
    {
      std::cerr << "Error: features not computed" << std::endl;
      return;
    }
  reset_indices();

  std::vector<int> training (m_points->point_set()->size(), -1);
  std::vector<int> indices (m_points->point_set()->size(), -1);

  std::vector<std::size_t> nb_label (m_labels.size(), 0);
  std::size_t nb_total = 0;
  
  for (Point_set::const_iterator it = m_points->point_set()->begin();
       it != m_points->point_set()->first_selected(); ++ it)
  {
    training[*it] = m_training[*it];
    if (training[*it] != -1)
    {
      nb_label[std::size_t(training[*it])] ++;
      ++ nb_total;
    }
  }

  std::cerr << nb_total << " point(s) used for training ("
            << 100. * (nb_total / double(m_points->point_set()->size())) << "% of the total):" << std::endl;
  for (std::size_t i = 0; i < m_labels.size(); ++ i)
    std::cerr << " * " << m_labels[i]->name() << ": " << nb_label[i] << " point(s)" << std::endl;

  if (classifier == 0)
    {
      m_sowf->train<Concurrency_tag>(training, nb_trials);
      CGAL::Classification::classify<Concurrency_tag> (*(m_points->point_set()),
                                                       m_labels, *m_sowf,
                                                       indices);
    }
  else if (classifier == 1)
  {
    m_ethz->train(training, true, num_trees, max_depth);
    CGAL::Classification::classify<Concurrency_tag> (*(m_points->point_set()),
                                                     m_labels, *m_ethz,
                                                     indices);
  }
  else
    {
#ifdef CGAL_LINKED_WITH_OPENCV
      if (m_random_forest != NULL)
        delete m_random_forest;
      m_random_forest = new Random_forest (m_labels, m_features,
                                           int(max_depth), 5, 15,
                                           int(num_trees));
      m_random_forest->train (training);
      CGAL::Classification::classify<Concurrency_tag> (*(m_points->point_set()),
                                                       m_labels, *m_random_forest,
                                                       indices);
#endif
    }
  for (Point_set::const_iterator it = m_points->point_set()->begin();
       it != m_points->point_set()->first_selected(); ++ it)
    m_classif[*it] = indices[*it];
  
  if (m_index_color == 1 || m_index_color == 2)
     change_color (m_index_color);
}

bool Point_set_item_classification::run (int method, int classifier,
                                         std::size_t subdivisions,
                                         double smoothing)
{
  if (m_features.size() == 0)
    {
      std::cerr << "Error: features not computed" << std::endl;
      return false;
    }
  reset_indices();

  if (classifier == 0)
    run (method, *m_sowf, subdivisions, smoothing);
  else if (classifier == 1)
    run (method, *m_ethz, subdivisions, smoothing);
#ifdef CGAL_LINKED_WITH_OPENCV
  else
    run (method, *m_random_forest, subdivisions, smoothing);
#endif
  
  return true;
}

