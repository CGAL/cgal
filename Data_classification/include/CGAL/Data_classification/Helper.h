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
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include <CGAL/Timer.h>

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
  typedef typename Kernel::Point_3                            Point;
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

  struct Scale
  {
    Neighborhood* neighborhood;
    Planimetric_grid* grid;
    Local_eigen_analysis* eigen;
    std::vector<Attribute_handle> attributes;
    double voxel_size;
    
    Scale (RandomAccessIterator begin, RandomAccessIterator end, PointMap point_map,
           const Iso_cuboid_3& bbox, double voxel_size)
      : voxel_size (voxel_size)
    {
      CGAL::Timer t;
      t.start();
      if (voxel_size < 0.)
        neighborhood = new Neighborhood (begin, end, point_map);
      else
        neighborhood = new Neighborhood (begin, end, point_map, voxel_size);
      t.stop();
      
      if (voxel_size < 0.)
        std::cerr << "Neighborhood computed in " << t.time() << " second(s)" << std::endl;
      else
        std::cerr << "Neighborhood with voxel size " << voxel_size
                  << " computed in " << t.time() << " second(s)" << std::endl;
      t.reset();
      t.start();
      double range;
      eigen = new Local_eigen_analysis (begin, end, point_map, *neighborhood, (std::size_t)6, range);
      if (this->voxel_size < 0)
        this->voxel_size = range / 3;
      t.stop();
      std::cerr << "Eigen values computed in " << t.time() << " second(s)" << std::endl;
      std::cerr << "Range = " << range << std::endl;
      t.reset();
      t.start();
      
      grid = new Planimetric_grid (begin, end, point_map, bbox, this->voxel_size);
      t.stop();
      std::cerr << "Planimetric grid computed in " << t.time() << " second(s)" << std::endl;
      t.reset();
    }
    ~Scale()
    {
      delete neighborhood;
      delete grid;
      delete eigen;
    }

    double grid_resolution() const { return voxel_size; }
    double radius_neighbors() const { return voxel_size * 5; }
    double radius_dtm() const { return voxel_size * 100; }
    
  };

  Iso_cuboid_3 m_bbox;
  std::vector<Scale*> m_scales;

public:

  Helper()
  {
  }

  Helper (RandomAccessIterator begin, RandomAccessIterator end, PointMap point_map,
          std::size_t nb_scales = 1)
  {
    m_bbox = CGAL::bounding_box
      (boost::make_transform_iterator (begin, CGAL::Property_map_to_unary_function<PointMap>(point_map)),
       boost::make_transform_iterator (end, CGAL::Property_map_to_unary_function<PointMap>(point_map)));

    CGAL::Timer t; t.start();

    m_scales.reserve (nb_scales);
    double voxel_size = - 1.;
    for (std::size_t i = 0; i < nb_scales; ++ i)
      {
        m_scales.push_back (new Scale (begin, end, point_map, m_bbox, voxel_size));
        if (i == 0)
          voxel_size = m_scales[0]->voxel_size;
        voxel_size *= 2;
      }
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
  {
    load (psc, filename, begin, end, point_map, normal_map, color_map, echo_map);
  }

  virtual ~Helper()
  {
    clear();
  }

  const Iso_cuboid_3& bbox() const { return m_bbox; }
  const Neighborhood& neighborhood() const { return (*m_scales[0]->neighborhood); }
  const Planimetric_grid& grid() const { return *(m_scales[0]->grid); }
  const Local_eigen_analysis& eigen() const { return *(m_scales[0]->eigen); }

  void info() const
  {
    std::cerr << m_scales.size() << " scale(s) used:" << std::endl;
    for (std::size_t i = 0; i < m_scales.size(); ++ i)
      {
        std::size_t nb_useful = 0;
        for (std::size_t j = 0; j < m_scales[i]->attributes.size(); ++ j)
          if (m_scales[i]->attributes[j]->weight != 0.)
            nb_useful ++;
        std::cerr << " * scale " << i << " with size " << m_scales[i]->voxel_size
                  << ", " << nb_useful << " useful attribute(s)";
        if (nb_useful != 0) std::cerr << ":" << std::endl;
        else std::cerr << std::endl;
        for (std::size_t j = 0; j < m_scales[i]->attributes.size(); ++ j)
          if (m_scales[i]->attributes[j]->weight != 0.)
            std::cerr << "   - " << m_scales[i]->attributes[j]->id()
                      << " (weight = " << m_scales[i]->attributes[j]->weight << ")" << std::endl;
      }
  }

  void clear()
  {
    for (std::size_t i = 0; i < m_scales.size(); ++ i)
      delete m_scales[i];
    m_scales.clear();
  }

  template <typename Attribute_type>
  void generate_multiscale_attribute_variant_0 (Point_set_classification& psc,
                                                RandomAccessIterator begin,
                                                RandomAccessIterator end,
                                                const PointMap& = PointMap())
  {
    for (std::size_t i = 0; i < m_scales.size(); ++ i)
      {
        psc.add_attribute (Attribute_handle (new Attribute_type(begin, end, *(m_scales[i]->eigen))));
        m_scales[i]->attributes.push_back (psc.get_attribute (psc.number_of_attributes() - 1));
      }
  }

  template <typename Attribute_type>
  void generate_multiscale_attribute_variant_1 (Point_set_classification& psc,
                                                RandomAccessIterator begin,
                                                RandomAccessIterator end,
                                                PointMap point_map)
  {
    for (std::size_t i = 0; i < m_scales.size(); ++ i)
      {
        psc.add_attribute (Attribute_handle (new Attribute_type(begin, end, point_map, *(m_scales[i]->eigen))));
        m_scales[i]->attributes.push_back (psc.get_attribute (psc.number_of_attributes() - 1));
      }
  }

  template <typename Attribute_type>
  void generate_multiscale_attribute_variant_2 (Point_set_classification& psc,
                                                RandomAccessIterator begin,
                                                RandomAccessIterator end,
                                                PointMap point_map)
  {
    for (std::size_t i = 0; i < m_scales.size(); ++ i)
      {
        psc.add_attribute (Attribute_handle (new Attribute_type(begin, end, point_map,
                                                                  *(m_scales[i]->grid),
                                                                  m_scales[i]->grid_resolution(),
                                                                  m_scales[i]->radius_neighbors())));
        m_scales[i]->attributes.push_back (psc.get_attribute (psc.number_of_attributes() - 1));
      }
  }

  template <typename Attribute_type>
  void generate_multiscale_attribute_variant_3 (Point_set_classification& psc,
                                                RandomAccessIterator begin,
                                                RandomAccessIterator end,
                                                PointMap point_map)
  {
    for (std::size_t i = 0; i < m_scales.size(); ++ i)
      {
        psc.add_attribute (Attribute_handle (new Attribute_type(begin, end, point_map,
                                                                m_bbox, *(m_scales[i]->grid),
                                                                m_scales[i]->grid_resolution(),
                                                                m_scales[i]->radius_neighbors(),
                                                                m_scales[i]->radius_dtm())));
        m_scales[i]->attributes.push_back (psc.get_attribute (psc.number_of_attributes() - 1));
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
    CGAL::Timer teigen, tpoint;
    teigen.start();
    generate_multiscale_attribute_variant_0<Anisotropy> (psc, begin, end, point_map);
    generate_multiscale_attribute_variant_1<Distance_to_plane> (psc, begin, end, point_map);
    generate_multiscale_attribute_variant_0<Eigentropy> (psc, begin, end, point_map);
    generate_multiscale_attribute_variant_0<Linearity> (psc, begin, end, point_map);
    generate_multiscale_attribute_variant_0<Omnivariance> (psc, begin, end, point_map);
    generate_multiscale_attribute_variant_0<Planarity> (psc, begin, end, point_map);
    generate_multiscale_attribute_variant_0<Sphericity> (psc, begin, end, point_map);
    generate_multiscale_attribute_variant_0<Sum_eigen> (psc, begin, end, point_map);
    generate_multiscale_attribute_variant_0<Surface_variation> (psc, begin, end, point_map);
    teigen.stop();
    tpoint.start();
    generate_multiscale_attribute_variant_2<Dispersion> (psc, begin, end, point_map);
    generate_multiscale_attribute_variant_3<Elevation> (psc, begin, end, point_map);
    tpoint.stop();
    
    std::cerr << "Point based attributes computed in " << tpoint.time() << " second(s)" << std::endl;
    std::cerr << "Eigen based attributes computed in " << teigen.time() << " second(s)" << std::endl;
  }

  template <typename VectorMap>
  void generate_normal_based_attributes(Point_set_classification& psc,
                                        RandomAccessIterator begin, RandomAccessIterator end,
                                        VectorMap normal_map)
  {
    CGAL::Timer t; t.start();
    psc.add_attribute (Attribute_handle (new Verticality(begin, end, normal_map)));
    m_scales[0]->attributes.push_back (psc.get_attribute (psc.number_of_attributes() - 1));
    t.stop();
    std::cerr << "Normal based attributes computed in " << t.time() << " second(s)" << std::endl;
  }

  void generate_normal_based_attributes(Point_set_classification& psc,
                                        RandomAccessIterator begin, RandomAccessIterator end,
                                        const CGAL::Empty_property_map<RandomAccessIterator, typename Kernel::Vector_3>&
                                        = CGAL::Empty_property_map<RandomAccessIterator, typename Kernel::Vector_3>())
  {
    CGAL::Timer t; t.start();
    generate_multiscale_attribute_variant_0<Verticality> (psc, begin, end);
    std::cerr << "Normal based attributes computed in " << t.time() << " second(s)" << std::endl;
  }

  template <typename ColorMap>
  void generate_color_based_attributes(Point_set_classification& psc,
                                       RandomAccessIterator begin, RandomAccessIterator end,
                                       ColorMap color_map)
  {

    typedef Attribute_hsv<Kernel, RandomAccessIterator, ColorMap> Hsv;
    CGAL::Timer t; t.start();
    for (std::size_t i = 0; i <= 8; ++ i)
      {
        psc.add_attribute (Attribute_handle (new Hsv(begin, end, color_map, 0, 45 * i, 22.5)));
        m_scales[0]->attributes.push_back (psc.get_attribute (psc.number_of_attributes() - 1));
      }
    for (std::size_t i = 0; i <= 4; ++ i)
      {
        psc.add_attribute (Attribute_handle (new Hsv(begin, end, color_map, 1, 25 * i, 12.5)));
        m_scales[0]->attributes.push_back (psc.get_attribute (psc.number_of_attributes() - 1));
      }
    for (std::size_t i = 0; i <= 4; ++ i)
      {
        psc.add_attribute (Attribute_handle (new Hsv(begin, end, color_map, 2, 25 * i, 12.5)));
        m_scales[0]->attributes.push_back (psc.get_attribute (psc.number_of_attributes() - 1));
      }
    t.stop();
    std::cerr << "Color based attributes computed in " << t.time() << " second(s)" << std::endl;
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
    CGAL::Timer t; t.start();
    for (std::size_t i = 0; i < m_scales.size(); ++ i)
      {
        psc.add_attribute (Attribute_handle (new Echo_scatter(begin, end, echo_map, *(m_scales[i]->grid),
                                                              m_scales[i]->grid_resolution(),
                                                              m_scales[i]->radius_neighbors())));
        m_scales[i]->attributes.push_back (psc.get_attribute (psc.number_of_attributes() - 1));
      }
    t.stop();
    std::cerr << "Echo based attributes computed in " << t.time() << " second(s)" << std::endl;
  }

  void generate_echo_based_attributes(const Point_set_classification&,
                                      RandomAccessIterator, RandomAccessIterator,
                                      const CGAL::Empty_property_map<RandomAccessIterator, std::size_t>&)
  {
  }

  void get_map_scale (std::map<Attribute_handle, std::size_t>& map_scale)
  {
    for (std::size_t i = 0; i < m_scales.size(); ++ i)
      for (std::size_t j = 0; j < m_scales[i]->attributes.size(); ++ j)
        map_scale[m_scales[i]->attributes[j]] = i;
  }

  std::size_t scale_of_attribute (Attribute_handle att)
  {
    for (std::size_t i = 0; i < m_scales.size(); ++ i)
      for (std::size_t j = 0; j < m_scales[i]->attributes.size(); ++ j)
        if (m_scales[i]->attributes[j] == att)
          return i;
    return (std::size_t)(-1);    
  }

  std::string name_att (Attribute_handle att,
                        std::map<Attribute_handle, std::size_t>& map_scale)
  {
    std::ostringstream oss;
    oss << att->id() << "_" << map_scale[att];
    return oss.str();
  }

  void save (const char* filename, Point_set_classification& psc)
  {
    boost::property_tree::ptree tree;

    // tree.put("classification.parameters.grid_resolution", m_grid_resolution);
    // tree.put("classification.parameters.radius_neighbors", m_radius_neighbors);
    tree.put("classification.parameters.voxel_size", m_scales[0]->voxel_size);

    std::map<Attribute_handle, std::size_t> map_scale;
    get_map_scale (map_scale);

    for (std::size_t i = 0; i < psc.number_of_attributes(); ++ i)
      {
        Attribute_handle att = psc.get_attribute(i);
        if (att->weight == 0)
          continue;
        boost::property_tree::ptree ptr;
        
        ptr.put("id", name_att (att, map_scale));
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
            ptr2.put("id", name_att (att, map_scale));
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
    
    m_bbox = CGAL::bounding_box
      (boost::make_transform_iterator (begin, CGAL::Property_map_to_unary_function<PointMap>(point_map)),
       boost::make_transform_iterator (end, CGAL::Property_map_to_unary_function<PointMap>(point_map)));

    boost::property_tree::ptree tree;
    boost::property_tree::read_xml(filename, tree);

    double voxel_size = tree.get<double>("classification.parameters.voxel_size");
    
    m_scales.push_back (new Scale (begin, end, point_map, m_bbox, voxel_size));

    std::map<std::string, Attribute_handle> att_map;
    BOOST_FOREACH(boost::property_tree::ptree::value_type &v, tree.get_child("classification.attributes"))
      {
        std::string full_id = v.second.get<std::string>("id");

        std::vector<std::string> splitted_id;
        boost::split(splitted_id, full_id, boost::is_any_of("_"));
        std::string id = splitted_id[0];
        for (std::size_t i = 1; i < splitted_id.size() - 1; ++ i)
          id = id + "_" + splitted_id[i];
        double scale = std::atof (splitted_id.back().c_str());

        while (m_scales.size() <= scale)
          {
            voxel_size *= 2;
            m_scales.push_back (new Scale (begin, end, point_map, m_bbox, voxel_size));
          }
        
        double weight = v.second.get<double>("weight");

        // Generate the right attribute if possible
        if (id == "anisotropy")
          psc.add_attribute (Attribute_handle (new Anisotropy(begin, end, *(m_scales[scale]->eigen))));
        else if (id == "distance_to_plane")
          psc.add_attribute (Attribute_handle (new Distance_to_plane(begin, end, point_map, *(m_scales[scale]->eigen))));
        else if (id == "eigentropy")
          psc.add_attribute (Attribute_handle (new Eigentropy(begin, end, *(m_scales[scale]->eigen))));
        else if (id == "elevation")
          psc.add_attribute (Attribute_handle (new Elevation(begin, end, point_map,
                                                             m_bbox, *(m_scales[scale]->grid),
                                                             m_scales[scale]->grid_resolution(),
                                                             m_scales[scale]->radius_neighbors(),
                                                             m_scales[scale]->radius_dtm())));
        else if (id == "linearity")
          psc.add_attribute (Attribute_handle (new Linearity(begin, end, *(m_scales[scale]->eigen))));
        else if (id == "omnivariance")
          psc.add_attribute (Attribute_handle (new Omnivariance(begin, end, *(m_scales[scale]->eigen))));
        else if (id == "planarity")
          psc.add_attribute (Attribute_handle (new Planarity(begin, end, *(m_scales[scale]->eigen))));
        else if (id == "sphericity")
          psc.add_attribute (Attribute_handle (new Sphericity(begin, end, *(m_scales[scale]->eigen))));
        else if (id == "sum_eigen")
          psc.add_attribute (Attribute_handle (new Sum_eigen(begin, end, *(m_scales[scale]->eigen))));
        else if (id == "surface_variation")
          psc.add_attribute (Attribute_handle (new Surface_variation(begin, end, *(m_scales[scale]->eigen))));
        else if (id == "vertical_dispersion")
          psc.add_attribute (Attribute_handle (new Dispersion(begin, end, point_map,
                                                              *(m_scales[scale]->grid),
                                                              m_scales[scale]->grid_resolution(),
                                                              m_scales[scale]->radius_neighbors())));
        else if (id == "verticality")
          {
            if (boost::is_convertible<VectorMap,
                typename CGAL::Empty_property_map<RandomAccessIterator, typename Kernel::Vector_3> >::value)
              psc.add_attribute (Attribute_handle (new Verticality(begin, end, *(m_scales[scale]->eigen))));
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
            psc.add_attribute (Attribute_handle (new Echo_scatter(begin, end, echo_map, *(m_scales[scale]->grid),
                                                                  m_scales[scale]->grid_resolution(),
                                                                  m_scales[scale]->radius_neighbors())));

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
        m_scales[scale]->attributes.push_back (att);
        att->weight = weight;
        att_map[full_id] = att;
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



  void write_ply (std::ostream& stream,
                  RandomAccessIterator begin, RandomAccessIterator end,
                  PointMap point_map, Point_set_classification& psc,
                  std::vector<RGB_Color>* colors = NULL)
  {
    stream << "ply" << std::endl
           << "format ascii 1.0" << std::endl
           << "element vertex " << (end - begin) << std::endl
           << "property double x" << std::endl
           << "property double y" << std::endl
           << "property double z" << std::endl
           << "property uchar red" << std::endl
           << "property uchar green" << std::endl
           << "property uchar blue" << std::endl
           << "property int label" << std::endl;

    bool delete_colors = false;
    if (colors == NULL)
      {
        colors = new std::vector<RGB_Color>();
        delete_colors = true;
      }
    
    std::map<Type_handle, std::size_t> map_types;
    for (std::size_t i = 0; i < psc.number_of_classification_types(); ++ i)
      {
        map_types.insert (std::make_pair (psc.get_classification_type(i), i));
        stream << "comment label " << i << " is " << psc.get_classification_type(i)->id() << std::endl;
      }
    map_types.insert (std::make_pair (Type_handle(), psc.number_of_classification_types()));
    
    stream << "end_header" << std::endl;

    while (colors->size() < psc.number_of_classification_types())
      {
        RGB_Color c = {{ (unsigned char)(192 + rand() % 60),
                         (unsigned char)(192 + rand() % 60),
                         (unsigned char)(192 + rand() % 60) }};
        colors->push_back (c);
      }
    RGB_Color c = {{ 0, 0, 0 }};
    colors->push_back (c);

    std::size_t i = 0;
    for (RandomAccessIterator it = begin; it != end; ++ it)
      {
        Type_handle t = psc.classification_type_of(i);
        std::size_t idx = map_types[t];
        
        stream << get(point_map, *it) << " "
               << (int)((*colors)[idx][0]) << " "
               << (int)((*colors)[idx][1]) << " "
               << (int)((*colors)[idx][2]) << " "
               << idx << std::endl;
        ++ i;
      }

    if (delete_colors)
      delete colors;
  }
  
};


} // namespace Data_classification

} // namespace CGAL


#endif // CGAL_DATA_CLASSIFICATION_HELPER_H
