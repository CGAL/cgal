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
    m_psc (NULL),
    m_trainer (NULL)
{
  m_nb_scales = 5;
  m_index_color = 1;
  m_nb_trials = 300;
  m_smoothing = 0.5;
  m_subdivisions = 16;

  reset_indices();
  
  m_psc = new PSC(*(m_points->point_set()), m_points->point_set()->point_map());

  backup_existing_colors_and_add_new();
  
  m_trainer = new Trainer(*m_psc);

  Label_handle ground = m_psc->add_label("ground");
  Label_handle vegetation = m_psc->add_label("vegetation");
  Label_handle roof = m_psc->add_label("roof");
  Label_handle facade = m_psc->add_label("facade");
  m_labels.push_back (std::make_pair(ground, QColor(245, 180, 0)));
  m_labels.push_back (std::make_pair(vegetation, QColor(0, 255, 27)));
  m_labels.push_back (std::make_pair(roof, QColor(255, 0, 170)));
  m_labels.push_back (std::make_pair(facade, QColor(100, 0, 255)));

}


Point_set_item_classification::~Point_set_item_classification()
{
  if (m_psc != NULL)
    delete m_psc;
  if (m_trainer != NULL)
    delete m_trainer;
  if (m_points != NULL)
    reset_colors();
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
bool Point_set_item_classification::write_ply_point_set(std::ostream& stream)
{
  if (m_psc->number_of_features() == 0)
    return false;

  reset_indices();
  
  stream.precision (std::numeric_limits<double>::digits10 + 2);

  std::vector<Color> colors;
  for (std::size_t i = 0; i < m_labels.size(); ++ i)
    {
      Color c = {{ (unsigned char)(m_labels[i].second.red()),
                   (unsigned char)(m_labels[i].second.green()),
                   (unsigned char)(m_labels[i].second.blue()) }};
      colors.push_back (c);
    }
  
  m_psc->write_classification_to_ply (stream);
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
      std::map<Label_handle, QColor> map_colors;
      for (std::size_t i = 0; i < m_labels.size(); ++ i)
        map_colors.insert (m_labels[i]);
          
      for (Point_set::const_iterator it = m_points->point_set()->begin();
           it != m_points->point_set()->first_selected(); ++ it)
        {
          QColor color (0, 0, 0);
          Label_handle c = m_psc->label_of(*it);
          
          if (c != Label_handle())
            color = map_colors[c];

          m_red[*it] = color.red();
          m_green[*it] = color.green();
          m_blue[*it] = color.blue();
        }
    }
  else if (index_color == 2) // training
    {
      std::map<Label_handle, QColor> map_colors;
      for (std::size_t i = 0; i < m_labels.size(); ++ i)
        map_colors.insert (m_labels[i]);
          
      for (Point_set::const_iterator it = m_points->point_set()->begin();
           it != m_points->point_set()->first_selected(); ++ it)
        {
          QColor color (0, 0, 0);
          Label_handle c = m_trainer->training_label_of(*it);
          Label_handle c2 = m_psc->label_of(*it);
          
          if (c != Label_handle())
            color = map_colors[c];
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
      Feature_handle att = m_psc->feature(index_color - 3);
      double weight = att->weight();
      att->set_weight(att->max);
      for (Point_set::const_iterator it = m_points->point_set()->begin();
           it != m_points->point_set()->first_selected(); ++ it)
        {
          m_red[*it] = (unsigned char)(ramp.r(att->normalized(*it)) * 255);
          m_green[*it] = (unsigned char)(ramp.g(att->normalized(*it)) * 255);
          m_blue[*it] = (unsigned char)(ramp.b(att->normalized(*it)) * 255);
        }
      att->set_weight(weight);
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
  CGAL_assertion (m_psc != NULL);
  m_psc->clear_features();
  reset_indices();
  
  std::cerr << "Computing features with " << m_nb_scales << " scale(s)" << std::endl;
  if (m_psc->number_of_features() != 0)
    m_psc->clear();

  bool normals = m_points->point_set()->has_normal_map();
  bool colors = (m_color != Point_set::Property_map<Color>());
  Point_set::Property_map<boost::uint8_t> echo_map;
  bool echo;
  boost::tie (echo_map, echo) = m_points->point_set()->template property_map<boost::uint8_t>("echo");

  if (!normals && !colors && !echo)
    m_psc->generate_features (m_nb_scales);
  else if (!normals && !colors && echo)
    m_psc->generate_features (m_nb_scales, CGAL::Default(), CGAL::Default(), echo_map);
  else if (!normals && colors && !echo)
    m_psc->generate_features (m_nb_scales, CGAL::Default(), m_color);
  else if (!normals && colors && echo)
    m_psc->generate_features (m_nb_scales, CGAL::Default(), m_color, echo_map);
  else if (normals && !colors && !echo)
    m_psc->generate_features (m_nb_scales, m_points->point_set()->normal_map());
  else if (normals && !colors && echo)
    m_psc->generate_features (m_nb_scales, m_points->point_set()->normal_map(), CGAL::Default(), echo_map);
  else if (normals && colors && !echo)
    m_psc->generate_features (m_nb_scales, m_points->point_set()->normal_map(), m_color);
  else
    m_psc->generate_features (m_nb_scales, m_points->point_set()->normal_map(), m_color, echo_map);

}



void Point_set_item_classification::train()
{
  if (m_psc->number_of_features() == 0)
    {
      std::cerr << "Error: features not computed" << std::endl;
      return;
    }
  reset_indices();
  
  m_trainer->train(m_nb_trials);
  m_psc->run();
  m_psc->info();
  if (m_index_color == 1 || m_index_color == 2)
    change_color (m_index_color);
}

bool Point_set_item_classification::run (int method)
{
  if (m_psc->number_of_features() == 0)
    {
      std::cerr << "Error: features not computed" << std::endl;
      return false;
    }
  reset_indices();

  if (method == 0)
    {
      m_psc->run();
      m_psc->info();
    }
  else if (method == 1)
    m_psc->run_with_local_smoothing (m_psc->neighborhood().range_neighbor_query(m_psc->radius_neighbors()));
  else if (method == 2)
    m_psc->run_with_graphcut (m_psc->neighborhood().k_neighbor_query(12), m_smoothing, m_subdivisions);
  
  if (m_index_color == 1 || m_index_color == 2)
    change_color (m_index_color);
  
  return true;
}
