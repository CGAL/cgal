#ifndef CGAL_DATA_CLASSIFICATION_HELPER_H
#define CGAL_DATA_CLASSIFICATION_HELPER_H

#include <CGAL/property_map.h>
#include <CGAL/Point_set_classification.h>
#include <CGAL/Data_classification/Planimetric_grid.h>
#include <CGAL/Data_classification/Local_eigen_analysis.h>
#include <CGAL/Data_classification/Attribute.h>
#include <CGAL/Data_classification/Attribute_color.h>
#include <CGAL/Data_classification/Attribute_distance_to_plane.h>
#include <CGAL/Data_classification/Attribute_echo_scatter.h>
#include <CGAL/Data_classification/Attribute_elevation.h>
#include <CGAL/Data_classification/Attribute_vertical_dispersion.h>
#include <CGAL/Data_classification/Attribute_verticality.h>
#include <CGAL/Data_classification/Attributes_eigen.h>
#include <CGAL/Data_classification/Type.h>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/lexical_cast.hpp>

namespace CGAL {

namespace Data_classification {
  

template <typename Kernel,
          typename RandomAccessIterator,
          typename PointMap,
          typename DiagonalizeTraits = CGAL::Default_diagonalize_traits<double,3> >
class Helper
{
  
public:
  typedef typename Kernel::Iso_cuboid_3                       Iso_cuboid_3;
  typedef CGAL::Point_set_classification
  <Kernel, RandomAccessIterator, PointMap>                    Point_set_classification;
  typedef CGAL::Data_classification::Planimetric_grid
  <Kernel, RandomAccessIterator, PointMap>                    Planimetric_grid;
  typedef CGAL::Data_classification::Neighborhood
  <Kernel, RandomAccessIterator, PointMap>                    Neighborhood;
  typedef CGAL::Data_classification::Local_eigen_analysis
  <Kernel, RandomAccessIterator, PointMap, DiagonalizeTraits> Local_eigen_analysis;
  
  typedef CGAL::Data_classification::Attribute_handle         Attribute_handle;
  typedef Attribute_anisotropy
  <Kernel, RandomAccessIterator, PointMap, DiagonalizeTraits> Anisotropy;
  typedef Attribute_distance_to_plane
  <Kernel, RandomAccessIterator, PointMap, DiagonalizeTraits> Distance_to_plane;
  typedef Attribute_eigentropy
  <Kernel, RandomAccessIterator, PointMap, DiagonalizeTraits> Eigentropy;
  typedef Attribute_elevation
  <Kernel, RandomAccessIterator, PointMap>                    Elevation;
  typedef Attribute_linearity
  <Kernel, RandomAccessIterator, PointMap, DiagonalizeTraits> Linearity;
  typedef Attribute_omnivariance
  <Kernel, RandomAccessIterator, PointMap, DiagonalizeTraits> Omnivariance;
  typedef Attribute_planarity
  <Kernel, RandomAccessIterator, PointMap, DiagonalizeTraits> Planarity;
  typedef Attribute_sphericity
  <Kernel, RandomAccessIterator, PointMap, DiagonalizeTraits> Sphericity;
  typedef Attribute_sum_eigenvalues
  <Kernel, RandomAccessIterator, PointMap, DiagonalizeTraits> Sum_eigen;
  typedef Attribute_surface_variation
  <Kernel, RandomAccessIterator, PointMap, DiagonalizeTraits> Surface_variation;
  typedef Attribute_vertical_dispersion
  <Kernel, RandomAccessIterator, PointMap>                    Dispersion;
  typedef Attribute_verticality
  <Kernel, RandomAccessIterator, PointMap, DiagonalizeTraits> Verticality;

  typedef CGAL::Data_classification::Type                     Type;
  typedef CGAL::Data_classification::Type_handle              Type_handle;
  
private:

  double m_grid_resolution;
  double m_radius_neighbors;
  double m_radius_dtm;
  Iso_cuboid_3 m_bbox;
  Planimetric_grid* m_grid;
  Neighborhood* m_neighborhood;
  Local_eigen_analysis* m_eigen;

public:

  Helper()
    : m_grid (NULL), m_neighborhood (NULL), m_eigen (NULL)
  {
  }

  Helper (RandomAccessIterator begin, RandomAccessIterator end, PointMap point_map,
          double grid_resolution,
          double radius_neighbors = -1.,
          double radius_dtm = -1.)
    : m_grid_resolution (grid_resolution)
    , m_radius_neighbors (radius_neighbors)
    , m_radius_dtm (radius_dtm)
  {
    if (m_radius_neighbors < 0.)
      m_radius_neighbors = 5. * m_grid_resolution;
    if (m_radius_dtm < 0.)
      m_radius_dtm = 5. * radius_neighbors;

    m_bbox = CGAL::bounding_box
      (boost::make_transform_iterator (begin, CGAL::Property_map_to_unary_function<PointMap>(point_map)),
       boost::make_transform_iterator (end, CGAL::Property_map_to_unary_function<PointMap>(point_map)));
    m_neighborhood = new Neighborhood (begin, end, point_map);
    m_grid = new Planimetric_grid (begin, end, point_map, m_bbox, m_grid_resolution);
    m_eigen = new Local_eigen_analysis (begin, end, point_map, *m_neighborhood, m_radius_neighbors);
  }

  
  template<typename VectorMap = CGAL::Empty_property_map<RandomAccessIterator, typename Kernel::Vector_3>,
           typename ColorMap = CGAL::Empty_property_map<RandomAccessIterator, RGB_Color>,
           typename EchoMap  = CGAL::Empty_property_map<RandomAccessIterator, std::size_t> >
  Helper (Point_set_classification& psc, const char* filename, 
          RandomAccessIterator begin, RandomAccessIterator end, 
          PointMap point_map,
          VectorMap normal_map = VectorMap(),
          ColorMap color_map = ColorMap(),
          EchoMap echo_map = EchoMap())
    : m_grid (NULL), m_neighborhood (NULL), m_eigen (NULL)
  {
    load (psc, filename, begin, end, point_map, normal_map, color_map, echo_map);
  }

  virtual ~Helper()
  {
    clear();
  }

  const Iso_cuboid_3& bbox() const { return m_bbox; }
  const Neighborhood& neighborhood() const { return *m_neighborhood; }
  const Planimetric_grid& grid() const { return *m_grid; }
  const Local_eigen_analysis& eigen() const { return *m_eigen; }

  
  void clear()
  {
    if (m_grid != NULL)
      {
        delete m_grid;
        delete m_neighborhood;
        delete m_eigen;
      }
  }

  template<typename VectorMap = CGAL::Empty_property_map<RandomAccessIterator, typename Kernel::Vector_3>,
           typename ColorMap = CGAL::Empty_property_map<RandomAccessIterator, RGB_Color>,
           typename EchoMap  = CGAL::Empty_property_map<RandomAccessIterator, std::size_t> >
  void generate_attributes(Point_set_classification& psc,
                           RandomAccessIterator begin, RandomAccessIterator end,
                           PointMap point_map,
                           VectorMap normal_map = VectorMap(),
                           ColorMap color_map = ColorMap(),
                           EchoMap echo_map = EchoMap())
  {
    generate_point_based_attributes (psc, begin, end, point_map);
    generate_normal_based_attributes (psc, begin, end, normal_map);
    generate_color_based_attributes (psc, begin, end, color_map);
    generate_echo_based_attributes (psc, begin, end, echo_map);
  }

  
  void generate_point_based_attributes (Point_set_classification& psc,
                                        RandomAccessIterator begin, RandomAccessIterator end,
                                        PointMap point_map)
  {
    psc.add_attribute (Attribute_handle (new Anisotropy(begin, end, *m_eigen)));
    psc.add_attribute (Attribute_handle (new Distance_to_plane(begin, end, point_map, *m_eigen)));
    psc.add_attribute (Attribute_handle (new Eigentropy(begin, end, *m_eigen)));
    psc.add_attribute (Attribute_handle (new Elevation(begin, end, point_map,
                                     m_bbox, *m_grid, m_grid_resolution,
                                     m_radius_neighbors, m_radius_dtm)));
    psc.add_attribute (Attribute_handle (new Linearity(begin, end, *m_eigen)));
    psc.add_attribute (Attribute_handle (new Omnivariance(begin, end, *m_eigen)));
    psc.add_attribute (Attribute_handle (new Planarity(begin, end, *m_eigen)));
    psc.add_attribute (Attribute_handle (new Sphericity(begin, end, *m_eigen)));
    psc.add_attribute (Attribute_handle (new Sum_eigen(begin, end, *m_eigen)));
    psc.add_attribute (Attribute_handle (new Surface_variation(begin, end, *m_eigen)));
    psc.add_attribute (Attribute_handle (new Dispersion(begin, end, point_map,
                                      *m_grid, m_grid_resolution, m_radius_neighbors)));
  }

  template <typename VectorMap>
  void generate_normal_based_attributes(Point_set_classification& psc,
                                        RandomAccessIterator begin, RandomAccessIterator end,
                                        VectorMap normal_map)
  {
    psc.add_attribute (Attribute_handle (new Verticality(begin, end, normal_map)));
  }

  void generate_normal_based_attributes(Point_set_classification& psc,
                                        RandomAccessIterator begin, RandomAccessIterator end,
                                        const CGAL::Empty_property_map<RandomAccessIterator, typename Kernel::Vector_3>&
                                        = CGAL::Empty_property_map<RandomAccessIterator, typename Kernel::Vector_3>())
  {
    psc.add_attribute (Attribute_handle (new Verticality(begin, end, *m_eigen)));
  }

  template <typename ColorMap>
  void generate_color_based_attributes(Point_set_classification& psc,
                                       RandomAccessIterator begin, RandomAccessIterator end,
                                       ColorMap color_map)
  {

    typedef Attribute_hsv<Kernel, RandomAccessIterator, ColorMap> Hsv;

    for (std::size_t i = 0; i <= 8; ++ i)
      psc.add_attribute (Attribute_handle (new Hsv(begin, end, color_map, 0, 45 * i, 22.5)));
    for (std::size_t i = 0; i <= 4; ++ i)
      psc.add_attribute (Attribute_handle (new Hsv(begin, end, color_map, 1, 25 * i, 12.5)));
    for (std::size_t i = 0; i <= 4; ++ i)
      psc.add_attribute (Attribute_handle (new Hsv(begin, end, color_map, 2, 25 * i, 12.5)));
  }

  void generate_color_based_attributes(const Point_set_classification&,
                                       RandomAccessIterator, RandomAccessIterator,
                                       const CGAL::Empty_property_map<RandomAccessIterator, RGB_Color>&)
  {
  }

  template <typename EchoMap>
  void generate_echo_based_attributes(Point_set_classification& psc,
                                      RandomAccessIterator begin, RandomAccessIterator end,
                                      EchoMap echo_map)
  {
    typedef Attribute_echo_scatter<Kernel, RandomAccessIterator, PointMap, EchoMap> Echo_scatter;
    psc.add_attribute (Attribute_handle (new Echo_scatter(begin, end, echo_map, *m_grid,
                                                          m_grid_resolution, m_radius_neighbors)));
  }

  void generate_echo_based_attributes(const Point_set_classification&,
                                      RandomAccessIterator, RandomAccessIterator,
                                      const CGAL::Empty_property_map<RandomAccessIterator, std::size_t>&)
  {
  }


  void save (const char* filename, Point_set_classification& psc)
  {
    boost::property_tree::ptree tree;

    tree.put("classification.parameters.grid_resolution", m_grid_resolution);
    tree.put("classification.parameters.radius_neighbors", m_radius_neighbors);
    tree.put("classification.parameters.radius_dtm", m_radius_dtm);


    for (std::size_t i = 0; i < psc.number_of_attributes(); ++ i)
      {
        Attribute_handle att = psc.get_attribute(i);
        if (att->weight == 0)
          continue;
        boost::property_tree::ptree ptr;
        ptr.put("id", att->id());
        ptr.put("weight", att->weight);
        tree.add_child("classification.attributes.attribute", ptr);
      }

    for (std::size_t i = 0; i < psc.number_of_classification_types(); ++ i)
      {
        Type_handle type = psc.get_classification_type(i);
        boost::property_tree::ptree ptr;
        ptr.put("id", type->id());
        for (std::size_t j = 0; j < psc.number_of_attributes(); ++ j)
          {
            Attribute_handle att = psc.get_attribute(j);
            if (att->weight == 0)
              continue;
            boost::property_tree::ptree ptr2;
            ptr2.put("id", att->id());
            Type::Attribute_effect effect = type->attribute_effect(att);
            if (effect == Type::Attribute_effect::PENALIZED_ATT)
              ptr2.put("effect", "penalized");
            else if (effect == Type::Attribute_effect::NEUTRAL_ATT)
              ptr2.put("effect", "neutral");
            else if (effect == Type::Attribute_effect::FAVORED_ATT)
              ptr2.put("effect", "favored");
            ptr.add_child("attribute", ptr2);
          }
        tree.add_child("classification.types.type", ptr);
      }

    // Write property tree to XML file
    boost::property_tree::xml_writer_settings<std::string> settings(' ', 3);
    boost::property_tree::write_xml(filename, tree, std::locale(), settings);
  }

  
  template<typename VectorMap = CGAL::Empty_property_map<RandomAccessIterator, typename Kernel::Vector_3>,
           typename ColorMap = CGAL::Empty_property_map<RandomAccessIterator, RGB_Color>,
           typename EchoMap  = CGAL::Empty_property_map<RandomAccessIterator, std::size_t> >
  bool load (Point_set_classification& psc, const char* filename, 
             RandomAccessIterator begin, RandomAccessIterator end, 
             PointMap point_map,
             VectorMap normal_map = VectorMap(),
             ColorMap color_map = ColorMap(),
             EchoMap echo_map = EchoMap())

  {
    typedef Attribute_echo_scatter<Kernel, RandomAccessIterator, PointMap, EchoMap> Echo_scatter;
    typedef Attribute_hsv<Kernel, RandomAccessIterator, ColorMap> Hsv;
    
    clear();
    
    boost::property_tree::ptree tree;

    boost::property_tree::read_xml(filename, tree);

    m_grid_resolution = tree.get<double>("classification.parameters.grid_resolution");
    m_radius_neighbors = tree.get<double>("classification.parameters.radius_neighbors");
    m_radius_dtm = tree.get<double>("classification.parameters.radius_dtm");

    m_bbox = CGAL::bounding_box
      (boost::make_transform_iterator (begin, CGAL::Property_map_to_unary_function<PointMap>(point_map)),
       boost::make_transform_iterator (end, CGAL::Property_map_to_unary_function<PointMap>(point_map)));
    m_neighborhood = new Neighborhood (begin, end, point_map);
    m_grid = new Planimetric_grid (begin, end, point_map, m_bbox, m_grid_resolution);
    m_eigen = new Local_eigen_analysis (begin, end, point_map, *m_neighborhood, m_radius_neighbors);

    std::map<std::string, Attribute_handle> att_map;
    BOOST_FOREACH(boost::property_tree::ptree::value_type &v, tree.get_child("classification.attributes"))
      {
        std::string id = v.second.get<std::string>("id");
        double weight = v.second.get<double>("weight");

        // Generate the right attribute if possible
        if (id == "anisotropy")
          psc.add_attribute (Attribute_handle (new Anisotropy(begin, end, *m_eigen)));
        else if (id == "distance_to_plane")
          psc.add_attribute (Attribute_handle (new Distance_to_plane(begin, end, point_map, *m_eigen)));
        else if (id == "eigentropy")
          psc.add_attribute (Attribute_handle (new Eigentropy(begin, end, *m_eigen)));
        else if (id == "elevation")
          psc.add_attribute (Attribute_handle (new Elevation(begin, end, point_map,
                                                             m_bbox, *m_grid, m_grid_resolution,
                                                             m_radius_neighbors, m_radius_dtm)));
        else if (id == "linearity")
          psc.add_attribute (Attribute_handle (new Linearity(begin, end, *m_eigen)));
        else if (id == "omnivariance")
          psc.add_attribute (Attribute_handle (new Omnivariance(begin, end, *m_eigen)));
        else if (id == "planarity")
          psc.add_attribute (Attribute_handle (new Planarity(begin, end, *m_eigen)));
        else if (id == "sphericity")
          psc.add_attribute (Attribute_handle (new Sphericity(begin, end, *m_eigen)));
        else if (id == "sum_eigen")
          psc.add_attribute (Attribute_handle (new Sum_eigen(begin, end, *m_eigen)));
        else if (id == "surface_variation")
          psc.add_attribute (Attribute_handle (new Surface_variation(begin, end, *m_eigen)));
        else if (id == "vertical_dispersion")
          psc.add_attribute (Attribute_handle (new Dispersion(begin, end, point_map,
                                                              *m_grid, m_grid_resolution, m_radius_neighbors)));
        else if (id == "verticality")
          {
            if (boost::is_convertible<VectorMap,
                typename CGAL::Empty_property_map<RandomAccessIterator, typename Kernel::Vector_3> >::value)
              psc.add_attribute (Attribute_handle (new Verticality(begin, end, *m_eigen)));
            else
              psc.add_attribute (Attribute_handle (new Verticality(begin, end, normal_map)));
          }
        else if (id == "echo_scatter")
          {
            if (boost::is_convertible<EchoMap,
                typename CGAL::Empty_property_map<RandomAccessIterator, std::size_t> >::value)
              {
                std::cerr << "Warning: echo_scatter required but no echo map given." << std::endl;
                continue;
              }
            psc.add_attribute (Attribute_handle (new Echo_scatter(begin, end, echo_map, *m_grid,
                                                                  m_grid_resolution, m_radius_neighbors)));

          }
        else if (boost::starts_with(id.c_str(), "hue")
                 || boost::starts_with(id.c_str(), "saturation")
                 || boost::starts_with(id.c_str(), "value"))
          {
            if (boost::is_convertible<ColorMap,
                typename CGAL::Empty_property_map<RandomAccessIterator, RGB_Color> >::value)
              {
                std::cerr << "Warning: color attribute required but no color map given." << std::endl;
                continue;
              }
            if (boost::starts_with(id.c_str(), "hue"))
              {
                double value = boost::lexical_cast<int>(id.c_str() + 4);
                psc.add_attribute (Attribute_handle (new Hsv(begin, end, color_map, 0, value, 22.5)));
              }
            else if (boost::starts_with(id.c_str(), "saturation"))
              {
                double value = boost::lexical_cast<int>(id.c_str() + 11);
                psc.add_attribute (Attribute_handle (new Hsv(begin, end, color_map, 1, value, 12.5)));
              }
            else if (boost::starts_with(id.c_str(), "value"))
              {
                double value = boost::lexical_cast<int>(id.c_str() + 6);
                psc.add_attribute (Attribute_handle (new Hsv(begin, end, color_map, 2, value, 12.5)));
              }
          }
        else
          {
            std::cerr << "Warning: unknown attribute \"" << id << "\"" << std::endl;
            continue;
          }

        Attribute_handle att = psc.get_attribute (psc.number_of_attributes() - 1);
        att->weight = weight;
        att_map[id] = att;
      }

    BOOST_FOREACH(boost::property_tree::ptree::value_type &v, tree.get_child("classification.types"))
      {
        std::string type_id = v.second.get<std::string>("id");

        Type_handle new_type = psc.add_classification_type (type_id.c_str());
        
        BOOST_FOREACH(boost::property_tree::ptree::value_type &v2, v.second)
          {
            if (v2.first == "id")
              continue;
            std::string att_id = v2.second.get<std::string>("id");
            std::map<std::string, Attribute_handle>::iterator it = att_map.find(att_id);
            if (it == att_map.end())
              continue;
            Attribute_handle att = it->second;
            std::string effect = v2.second.get<std::string>("effect");
            if (effect == "penalized")
              new_type->set_attribute_effect (att, Type::PENALIZED_ATT);
            else if (effect == "neutral")
              new_type->set_attribute_effect (att, Type::NEUTRAL_ATT);
            else
              new_type->set_attribute_effect (att, Type::FAVORED_ATT);
          }
      }
    
    return true;
  }
  
};


} // namespace Data_classification

} // namespace CGAL


#endif // CGAL_DATA_CLASSIFICATION_HELPER_H
