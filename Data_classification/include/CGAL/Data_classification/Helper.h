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
    typedef Attribute_color<Kernel, RandomAccessIterator, ColorMap> Color;
    psc.add_attribute (Attribute_handle (new Color(begin, end, color_map)));
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

  
};


} // namespace Data_classification

} // namespace CGAL


#endif // CGAL_DATA_CLASSIFICATION_HELPER_H
