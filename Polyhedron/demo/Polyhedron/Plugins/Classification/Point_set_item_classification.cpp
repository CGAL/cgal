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

  backup_existing_colors_and_add_new();
  bool training_found = false;
  boost::tie (m_training, training_found) = m_points->point_set()->add_property_map<int>("training", -1);
  bool classif_found = false;
  boost::tie (m_classif, classif_found) = m_points->point_set()->add_property_map<int>("label", -1);
  
  training_found = !training_found; // add_property_map returns false if 
  classif_found = !classif_found;   // property was already there
  
  if (training_found || classif_found)
  {
    int lmax = 0;
    
    for (Point_set::const_iterator it = m_points->point_set()->begin();
         it != m_points->point_set()->first_selected(); ++ it)
    {
      if (training_found)
      {
        int l = m_training[*it];
        lmax = (std::max)(l, lmax);
      }
      if (classif_found)
      {
        int l = m_classif[*it];
        lmax = (std::max)(l, lmax);
      }
    }

    // Try to recover label names from PLY comments
    std::map<int, std::string> label_names;
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
    
    for (int i = 0; i <= lmax; ++ i)
    {
      std::map<int, std::string>::iterator found
        = label_names.find (i);
      if (found != label_names.end())
        m_labels.add(found->second.c_str());
      else
      {
        std::ostringstream oss;
        oss << "label_" << i;
        m_labels.add(oss.str().c_str());
      }
      
      CGAL::Classification::HSV_Color hsv;
      hsv[0] = 360. * (i / double(lmax + 1));
      hsv[1] = 76.;
      hsv[2] = 85.;
      Color rgb = CGAL::Classification::hsv_to_rgb(hsv);
      m_label_colors.push_back (QColor(rgb[0], rgb[1], rgb[2]));
    }
  }
  else
  {
    m_labels.add("ground");
    m_labels.add("vegetation");
    m_labels.add("roof");
    m_labels.add("facade");

    m_label_colors.push_back (QColor(245, 180, 0));
    m_label_colors.push_back (QColor(0, 255, 27));
    m_label_colors.push_back (QColor(255, 0, 170));
    m_label_colors.push_back (QColor(100, 0, 255));
  }
  
  update_comments_of_point_set_item();

  m_sowf = new Sum_of_weighted_features (m_labels, m_features);
#ifdef CGAL_LINKED_WITH_OPENCV
  m_random_forest = new Random_forest (m_labels, m_features);
#endif
}


Point_set_item_classification::~Point_set_item_classification()
{
  if (m_sowf != NULL)
    delete m_sowf;
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
  ramp.build_red();
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

  if (!normals && !colors && !echo)
    m_generator = new Generator (m_features, *(m_points->point_set()), m_points->point_set()->point_map(), nb_scales);
  else if (!normals && !colors && echo)
    m_generator = new Generator (m_features, *(m_points->point_set()), m_points->point_set()->point_map(), nb_scales,
                                 CGAL::Default(), CGAL::Default(), echo_map);
  else if (!normals && colors && !echo)
    m_generator = new Generator (m_features, *(m_points->point_set()), m_points->point_set()->point_map(), nb_scales,
                                 CGAL::Default(), m_color);
  else if (!normals && colors && echo)
    m_generator = new Generator (m_features, *(m_points->point_set()), m_points->point_set()->point_map(), nb_scales,
                                 CGAL::Default(), m_color, echo_map);
  else if (normals && !colors && !echo)
    m_generator = new Generator (m_features, *(m_points->point_set()), m_points->point_set()->point_map(), nb_scales,
                                 m_points->point_set()->normal_map());
  else if (normals && !colors && echo)
    m_generator = new Generator (m_features, *(m_points->point_set()), m_points->point_set()->point_map(), nb_scales,
                                 m_points->point_set()->normal_map(), CGAL::Default(), echo_map);
  else if (normals && colors && !echo)
    m_generator = new Generator (m_features, *(m_points->point_set()), m_points->point_set()->point_map(), nb_scales,
                                 m_points->point_set()->normal_map(), m_color);
  else
    m_generator = new Generator (m_features, *(m_points->point_set()), m_points->point_set()->point_map(), nb_scales,
                                 m_points->point_set()->normal_map(), m_color, echo_map);

  delete m_sowf;
  m_sowf = new Sum_of_weighted_features (m_labels, m_features);
#ifdef CGAL_LINKED_WITH_OPENCV
  delete m_random_forest;
  m_random_forest = new Random_forest (m_labels, m_features);
#endif
  std::cerr << "Features = " << m_features.size() << std::endl;
}



void Point_set_item_classification::train(int classifier, unsigned int nb_trials)
{
  if (m_features.size() == 0)
    {
      std::cerr << "Error: features not computed" << std::endl;
      return;
    }
  reset_indices();

  std::vector<int> training (m_points->point_set()->size(), -1);
  std::vector<std::size_t> indices (m_points->point_set()->size(), std::size_t(-1));

  for (Point_set::const_iterator it = m_points->point_set()->begin();
       it != m_points->point_set()->first_selected(); ++ it)
    training[*it] = m_training[*it];

  if (classifier == 0)
    {
      m_sowf->train<Concurrency_tag>(training, nb_trials);
      CGAL::Classification::classify<Concurrency_tag> (*(m_points->point_set()),
                                                       m_labels, *m_sowf,
                                                       indices);
    }
  else
    {
#ifdef CGAL_LINKED_WITH_OPENCV
      m_random_forest->train (training);
      CGAL::Classification::classify<Concurrency_tag> (*(m_points->point_set()),
                                                       m_labels, *m_random_forest,
                                                       indices);
#endif
    }
  for (Point_set::const_iterator it = m_points->point_set()->begin();
       it != m_points->point_set()->first_selected(); ++ it)
    m_classif[*it] = int(indices[*it]);
  
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
#ifdef CGAL_LINKED_WITH_OPENCV
  else
    run (method, *m_random_forest, subdivisions, smoothing);
#endif
  
  return true;
}

