#include "Point_set_item_classification.h"
#include "Color_ramp.h"

#include <CGAL/Timer.h>
#include <CGAL/Memory_sizer.h>

#include <CGAL/Three/Viewer_interface.h>

#include <set>
#include <stack>
#include <algorithm>
#include <boost/array.hpp>

#ifdef CGAL_LINKED_WITH_TBB
typedef CGAL::Parallel_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif

Point_set_item_classification::Point_set_item_classification(Scene_points_with_normal_item* points)
  : m_points (points),
    m_generator (NULL)
{
  m_nb_scales = 5;
  m_index_color = 1;
  m_nb_trials = 300;
  m_smoothing = 0.5;
  m_subdivisions = 16;

  reset_indices();
  
  backup_existing_colors_and_add_new();
  m_training = m_points->point_set()->add_property_map<std::size_t>("training", std::size_t(-1)).first;
  m_classif = m_points->point_set()->add_property_map<std::size_t>("label", std::size_t(-1)).first;

  m_labels.add("ground");
  m_labels.add("vegetation");
  m_labels.add("roof");
  m_labels.add("facade");
  
  m_label_colors.push_back (QColor(245, 180, 0));
  m_label_colors.push_back (QColor(0, 255, 27));
  m_label_colors.push_back (QColor(255, 0, 170));
  m_label_colors.push_back (QColor(100, 0, 255));

  m_predicate = new Predicate (m_labels, m_features);
}


Point_set_item_classification::~Point_set_item_classification()
{
  if (m_predicate != NULL)
    delete m_predicate;
  if (m_generator != NULL)
    delete m_generator;
  if (m_points != NULL)
    {
      reset_colors();
      m_points->point_set()->remove_property_map(m_training);
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
    {
      m_points->point_set()->remove_property_map(m_red);
      m_points->point_set()->remove_property_map(m_green);
      m_points->point_set()->remove_property_map(m_blue);
      m_points->point_set()->check_colors();
    }
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
          std::size_t c = m_training[*it];
          std::size_t c2 = m_classif[*it];
          
          if (c != std::size_t(-1))
            color = m_label_colors[c];
          
          double div = 1;
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

      double max = 0.;
      for (Point_set::const_iterator it = m_points->point_set()->begin();
           it != m_points->point_set()->first_selected(); ++ it)
        if (feature->value(*it) > max)
          max = feature->value(*it);

      for (Point_set::const_iterator it = m_points->point_set()->begin();
           it != m_points->point_set()->first_selected(); ++ it)
        {
          double v = std::max (0., feature->value(*it) / max);
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

void Point_set_item_classification::compute_features ()
{
  CGAL_assertion (!(m_points->point_set()->empty()));

  if (m_generator != NULL)
    delete m_generator;

  reset_indices();
  
  std::cerr << "Computing features with " << m_nb_scales << " scale(s)" << std::endl;
  m_features.clear();

  bool normals = m_points->point_set()->has_normal_map();
  bool colors = (m_color != Point_set::Property_map<Color>());
  Point_set::Property_map<boost::uint8_t> echo_map;
  bool echo;
  boost::tie (echo_map, echo) = m_points->point_set()->template property_map<boost::uint8_t>("echo");

  if (!normals && !colors && !echo)
    m_generator = new Generator (m_features, m_nb_scales, *(m_points->point_set()), m_points->point_set()->point_map());
  else if (!normals && !colors && echo)
    m_generator = new Generator (m_features, m_nb_scales, *(m_points->point_set()), m_points->point_set()->point_map(),
                                 CGAL::Default(), CGAL::Default(), echo_map);
  else if (!normals && colors && !echo)
    m_generator = new Generator (m_features, m_nb_scales, *(m_points->point_set()), m_points->point_set()->point_map(),
                                 CGAL::Default(), m_color);
  else if (!normals && colors && echo)
    m_generator = new Generator (m_features, m_nb_scales, *(m_points->point_set()), m_points->point_set()->point_map(),
                                 CGAL::Default(), m_color, echo_map);
  else if (normals && !colors && !echo)
    m_generator = new Generator (m_features, m_nb_scales, *(m_points->point_set()), m_points->point_set()->point_map(),
                                 m_points->point_set()->normal_map());
  else if (normals && !colors && echo)
    m_generator = new Generator (m_features, m_nb_scales, *(m_points->point_set()), m_points->point_set()->point_map(),
                                 m_points->point_set()->normal_map(), CGAL::Default(), echo_map);
  else if (normals && colors && !echo)
    m_generator = new Generator (m_features, m_nb_scales, *(m_points->point_set()), m_points->point_set()->point_map(),
                                 m_points->point_set()->normal_map(), m_color);
  else
    m_generator = new Generator (m_features, m_nb_scales, *(m_points->point_set()), m_points->point_set()->point_map(),
                                 m_points->point_set()->normal_map(), m_color, echo_map);

  delete m_predicate;
  m_predicate = new Predicate (m_labels, m_features);
  std::cerr << "Features = " << m_features.size() << std::endl;
}



void Point_set_item_classification::train()
{
  if (m_features.size() == 0)
    {
      std::cerr << "Error: features not computed" << std::endl;
      return;
    }
  reset_indices();

  std::vector<std::size_t> indices (m_points->point_set()->size(), std::size_t(-1));

  for (Point_set::const_iterator it = m_points->point_set()->begin();
       it != m_points->point_set()->first_selected(); ++ it)
    indices[*it] = m_training[*it];

  m_predicate->train<Concurrency_tag>(indices, m_nb_trials);
  CGAL::Classification::classify<Concurrency_tag> (*(m_points->point_set()),
                                                   m_labels, *m_predicate,
                                                   indices);
  for (Point_set::const_iterator it = m_points->point_set()->begin();
       it != m_points->point_set()->first_selected(); ++ it)
    m_classif[*it] = indices[*it];
  
  if (m_index_color == 1 || m_index_color == 2)
     change_color (m_index_color);
}

bool Point_set_item_classification::run (int method)
{
  if (m_features.size() == 0)
    {
      std::cerr << "Error: features not computed" << std::endl;
      return false;
    }
  reset_indices();

  std::vector<std::size_t> indices;

  if (method == 0)
    CGAL::Classification::classify<Concurrency_tag> (*(m_points->point_set()),
                                                     m_labels, *m_predicate,
                                                     indices);
  else if (method == 1)
    CGAL::Classification::classify_with_local_smoothing<Concurrency_tag>
      (*(m_points->point_set()), m_points->point_set()->point_map(), m_labels, *m_predicate,
       m_generator->neighborhood().range_neighbor_query(m_generator->radius_neighbors()),
       indices);
  else if (method == 2)
    CGAL::Classification::classify_with_graphcut<Concurrency_tag>
      (*(m_points->point_set()), m_points->point_set()->point_map(),
       m_points->point_set()->point_map(), m_labels, *m_predicate,
       m_generator->neighborhood().k_neighbor_query(12),
       m_smoothing, m_subdivisions, indices);
  
  for (Point_set::const_iterator it = m_points->point_set()->begin();
       it != m_points->point_set()->first_selected(); ++ it)
    m_classif[*it] = indices[*it];
  
  if (m_index_color == 1 || m_index_color == 2)
    change_color (m_index_color);
  
  return true;
}
