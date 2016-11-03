// Copyright (c) 2016  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_CLASSIFICATION_HELPER_H
#define CGAL_CLASSIFICATION_HELPER_H

#include <CGAL/property_map.h>
#include <CGAL/Classifier.h>
#include <CGAL/Classification/Point_set_neighborhood.h>
#include <CGAL/Classification/Planimetric_grid.h>
#include <CGAL/Classification/Local_eigen_analysis.h>
#include <CGAL/Classification/Attribute.h>
#include <CGAL/Classification/Attribute_color.h>
#include <CGAL/Classification/Attribute_distance_to_plane.h>
#include <CGAL/Classification/Attribute_echo_scatter.h>
#include <CGAL/Classification/Attribute_elevation.h>
#include <CGAL/Classification/Attribute_vertical_dispersion.h>
#include <CGAL/Classification/Attribute_verticality.h>
#include <CGAL/Classification/Attributes_eigen.h>
#include <CGAL/Classification/Type.h>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include <CGAL/Timer.h>
#include <CGAL/demangle.h>

namespace CGAL {

namespace Classification {
  

/*!
\ingroup PkgClassification

\brief Class designed to help the user using the classification
algorithm.

The classification algorithm is designed to be as flexible as
possible, letting the user free to use its own data structures,
attributes and classification types.

Nevertheless, \cgal provides a predefined framework that should work
correctly on common urban point sets. Using this class, classification
can be performed without having to instanciate all attributes and data
structures separately.

This class also provides functions to save and load a specific
configuration (attribute effects and weights), as well as a function
to write a classified point set in a PLY format with colors and
labels.

\tparam Kernel The geometric kernel used
\tparam RandomAccessIterator Iterator over the input
\tparam PointMap is a model of `ReadablePropertyMap` with value type `Point_3<Kernel>`.
\tparam DiagonalizeTraits Solver used for matrix diagonalization.
*/
template <typename Kernel,
          typename RandomAccessIterator,
          typename PointMap,
          typename DiagonalizeTraits = CGAL::Default_diagonalize_traits<double,3> >
class Helper
{
  
public:
  typedef typename Kernel::Iso_cuboid_3                       Iso_cuboid_3;

  typedef CGAL::Classifier
  <RandomAccessIterator, PointMap>                            Classifier;
  typedef CGAL::Classification::Planimetric_grid
  <Kernel, RandomAccessIterator, PointMap>                    Planimetric_grid;
  typedef CGAL::Classification::Point_set_neighborhood
  <Kernel, RandomAccessIterator, PointMap>                    Neighborhood;
  typedef CGAL::Classification::Local_eigen_analysis
  <Kernel, RandomAccessIterator, PointMap, DiagonalizeTraits> Local_eigen_analysis;
  
  typedef CGAL::Classification::Attribute_handle         Attribute_handle;
  typedef CGAL::Classification::Type                     Type;
  typedef CGAL::Classification::Type_handle              Type_handle;
  

  /// \cond SKIP_IN_MANUAL
  typedef typename Kernel::Point_3                            Point;
  
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
  /// \endcond
  
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
      eigen = new Local_eigen_analysis (begin, end, point_map, neighborhood->k_neighbor_query(6), range);
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
  /// \cond SKIP_IN_MANUAL
  Helper()
  {
  }
  /// \endcond
  

  /*!
    \brief Constructs an helper object.

    This triggers the instanciation of all necessary data structures
    (neighborhood, eigen analysis, etc.) on the specified number of
    scales.

    \param begin Iterator to the first input object
    \param end Past-the-end iterator
    \param point_map Property map to access the input points

    \param nb_scales Number of scales. Default only uses the first
    scale (unaltered point set). Increasing this value generates
    recursively simplified point set and computes the associated
    attributes. Using `nb_scales=5` usually provides satisfying
    results. Using more than 10 scales is not recommended as it does
    not help generating better results.
  */
  Helper (RandomAccessIterator begin, RandomAccessIterator end, PointMap point_map,
          std::size_t nb_scales = 1)
  {
    m_bbox = CGAL::bounding_box
      (boost::make_transform_iterator (begin, CGAL::Property_map_to_unary_function<PointMap>(point_map)),
       boost::make_transform_iterator (end, CGAL::Property_map_to_unary_function<PointMap>(point_map)));

    CGAL::Timer t; t.start();

    m_scales.reserve (nb_scales);
    double voxel_size = - 1.;

    m_scales.push_back (new Scale (begin, end, point_map, m_bbox, voxel_size));
    voxel_size = m_scales[0]->grid_resolution();
    
    for (std::size_t i = 1; i < nb_scales; ++ i)
      {
        voxel_size *= 2;
        m_scales.push_back (new Scale (begin, end, point_map, m_bbox, voxel_size));
      }
  }


  /*!
    \brief Constructs an helper object by loading a configuration file.

    All data structures, attributes and types specified in the input
    file `filename` are instanciated if possible (in particular,
    property maps needed should be provided).

    \tparam VectorMap is a model of `ReadablePropertyMap` with value type `Vector_3<Kernel>`.
    \tparam ColorMap is a model of `ReadablePropertyMap` with value type `CGAL::Classification::RGB_Color`.
    \tparam EchoMap is a model of `ReadablePropertyMap` with value type `std::size_`.
    \param filename Name of the output file
    \param psc Classification object where to store attributes and types
    \param begin Iterator to the first input object
    \param end Past-the-end iterator
    \param point_map Property map to access the input points
    \param normal_map Property map to access the normal vectors of the input points (if any).
    \param color_map Property map to access the colors of the input points (if any).
    \param echo_map Property map to access the echo values of the input points (if any).
  */
  template<typename VectorMap = CGAL::Empty_property_map<RandomAccessIterator, typename Kernel::Vector_3>,
           typename ColorMap = CGAL::Empty_property_map<RandomAccessIterator, RGB_Color>,
           typename EchoMap  = CGAL::Empty_property_map<RandomAccessIterator, std::size_t> >
  Helper (const char* filename, Classifier& psc, 
          RandomAccessIterator begin, RandomAccessIterator end, 
          PointMap point_map,
          VectorMap normal_map = VectorMap(),
          ColorMap color_map = ColorMap(),
          EchoMap echo_map = EchoMap())
  {
    load (filename, psc, begin, end, point_map, normal_map, color_map, echo_map);
  }

  /// \cond SKIP_IN_MANUAL
  virtual ~Helper()
  {
    clear();
  }
  /// \endcond

  /*!
    \brief Returns the bounding box of the input point set.
  */
  const Iso_cuboid_3& bbox() const { return m_bbox; }
  /*!
    \brief Returns the neighborhood structure at scale `scale`.
  */
  const Neighborhood& neighborhood(std::size_t scale = 0) const { return (*m_scales[scale]->neighborhood); }
  /*!
    \brief Returns the planimetric grid structure at scale `scale`.
  */
  const Planimetric_grid& grid(std::size_t scale = 0) const { return *(m_scales[scale]->grid); }
  /*!
    \brief Returns the local eigen analysis structure at scale `scale`.
  */
  const Local_eigen_analysis& eigen(std::size_t scale = 0) const { return *(m_scales[scale]->eigen); }
  /*!
    \brief Returns the grid resolution at scale `scale`.
  */
  double grid_resolution(std::size_t scale = 0) const { return m_scales[scale]->grid_resolution(); }
  /*!
    \brief Returns the radius used for neighborhood queries at scale `scale`.
  */
  double radius_neighbors(std::size_t scale = 0) const { return m_scales[scale]->radius_neighbors(); }
  /*!
    \brief Returns the radius used for digital terrain modeling at scale `scale`.
  */
  double radius_dtm(std::size_t scale = 0) const { return m_scales[scale]->radius_dtm(); }

  /// \cond SKIP_IN_MANUAL  
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
  /// \endcond

  /*!
    \brief Clears all computed data structures.
  */
  void clear()
  {
    for (std::size_t i = 0; i < m_scales.size(); ++ i)
      delete m_scales[i];
    m_scales.clear();
  }

  /*!
    \brief Generate all possible attributes from an input range.

    This method calls `generate_point_based_attributes()`,
    `generate_normal_based_attributes()`,
    `generate_color_based_attributes()` and
    `generate_echo_based_attributes()`.

    If a property map is left to its default
    `CGAL::Empty_property_map` type, the corresponding attributes are
    not computed (this method can thus be called even if a point set
    does not have associated normal, color or echo properties).

    \tparam VectorMap is a model of `ReadablePropertyMap` with value type `Vector_3<Kernel>`.
    \tparam ColorMap is a model of `ReadablePropertyMap` with value type `CGAL::Classification::RGB_Color`.
    \tparam EchoMap is a model of `ReadablePropertyMap` with value type `std::size_`.
    \param psc The classification object where to store the attributes
    \param begin Iterator to the first input object
    \param end Past-the-end iterator
    \param point_map Property map to access the input points
    \param normal_map Property map to access the normal vectors of the input points (if any).
    \param color_map Property map to access the colors of the input points (if any).
    \param echo_map Property map to access the echo values of the input points (if any).
  */
  template<typename VectorMap = CGAL::Empty_property_map<RandomAccessIterator, typename Kernel::Vector_3>,
           typename ColorMap = CGAL::Empty_property_map<RandomAccessIterator, RGB_Color>,
           typename EchoMap  = CGAL::Empty_property_map<RandomAccessIterator, std::size_t> >
  void generate_attributes(Classifier& psc,
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

  
  /*!
    \brief Generate all points attributes from an input range.

    Generate, for all precomputed scales, the following attributes:

    - `CGAL::Classification::Attribute_anisotropy`
    - `CGAL::Classification::Attribute_distance_to_plane`
    - `CGAL::Classification::Attribute_eigentropy`
    - `CGAL::Classification::Attribute_elevation`
    - `CGAL::Classification::Attribute_linearity`
    - `CGAL::Classification::Attribute_omnivariance`
    - `CGAL::Classification::Attribute_planarity`
    - `CGAL::Classification::Attribute_sphericity`
    - `CGAL::Classification::Attribute_sum_eigenvalues`
    - `CGAL::Classification::Attribute_surface_variation`
    - `CGAL::Classification::Attribute_vertical_dispersion`

    \param psc The classification object where to store the attributes
    \param begin Iterator to the first input object
    \param end Past-the-end iterator
    \param point_map Property map to access the input points
  */
  void generate_point_based_attributes (Classifier& psc,
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

  /*!
    \brief Generate all normal attributes from an input range.

    Generate, for all precomputed scales, the following attribute:

    - `CGAL::Classification::Attribute_verticality`

    If the normal map is left to its default type
    `CGAL::Empty_property_map`, then the verticality attributes are
    still computed by using an approximation of the normal vector
    provided by the corresponding `Local_eigen_analysis` object.

    \tparam VectorMap Property map to access the normal vectors of the input points (if any).
    \param psc The classification object where to store the attributes
    \param begin Iterator to the first input object
    \param end Past-the-end iterator
    \param normal_map Property map to access the normal vectors of the input points (if any).
  */
  #ifdef DOXYGEN_RUNNING
  template<typename VectorMap = CGAL::Empty_property_map<RandomAccessIterator, typename Kernel::Vector_3> >
  #else
  template <typename VectorMap>
  #endif
  void generate_normal_based_attributes(Classifier& psc,
                                        RandomAccessIterator begin, RandomAccessIterator end,
#ifdef DOXYGEN_RUNNING
                                        VectorMap normal_map = VectorMap())
#else
                                        VectorMap normal_map)
#endif
  {
    CGAL::Timer t; t.start();
    psc.add_attribute (Attribute_handle (new Verticality(begin, end, normal_map)));
    m_scales[0]->attributes.push_back (psc.get_attribute (psc.number_of_attributes() - 1));
    t.stop();
    std::cerr << "Normal based attributes computed in " << t.time() << " second(s)" << std::endl;
  }

  /// \cond SKIP_IN_MANUAL
  void generate_normal_based_attributes(Classifier& psc,
                                        RandomAccessIterator begin, RandomAccessIterator end,
                                        const CGAL::Empty_property_map<RandomAccessIterator, typename Kernel::Vector_3>&
                                        = CGAL::Empty_property_map<RandomAccessIterator, typename Kernel::Vector_3>())
  {
    CGAL::Timer t; t.start();
    generate_multiscale_attribute_variant_0<Verticality> (psc, begin, end);
    std::cerr << "Normal based attributes computed in " << t.time() << " second(s)" << std::endl;
  }
  /// \endcond

  /*!
    \brief Generate a set of color attributes from an input range.

    Generate the following attributes:

    - 9 attributes `CGAL::Classification::Attribute_hsv` on
      channel 0 (hue) with mean ranging from 0° to 360° and standard
      deviation of 22.5.

    - 5 attributes `CGAL::Classification::Attribute_hsv` on
      channel 1 (saturation) with mean ranging from 0 to 100 and standard
      deviation of 12.5

    - 5 attributes `CGAL::Classification::Attribute_hsv` on
      channel 2 (value) with mean ranging from 0 to 100 and standard
      deviation of 12.5

    This decomposition allows to handle all the color spectrum with a
    usually sufficiently accurate precision.

    \tparam ColorMap is a model of `ReadablePropertyMap` with value type `CGAL::Classification::RGB_Color`.
    \param psc The classification object where to store the attributes
    \param begin Iterator to the first input object
    \param end Past-the-end iterator
    \param color_map Property map to access the colors of the input points.
  */
  template <typename ColorMap>
  void generate_color_based_attributes(Classifier& psc,
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

  /// \cond SKIP_IN_MANUAL
  void generate_color_based_attributes(const Classifier&,
                                       RandomAccessIterator, RandomAccessIterator,
                                       const CGAL::Empty_property_map<RandomAccessIterator, RGB_Color>&)
  {
  }
  /// \endcond

  /*!
    \brief Generate all echo attributes from an input range.

    Generate, for all precomputed scales, the following attribute:

    - `CGAL::Classification::Attribute_echo_scatter`

    \tparam EchoMap Property map to access the echo values of the input points (if any).
    \param psc The classification object where to store the attributes
    \param begin Iterator to the first input object
    \param end Past-the-end iterator
    \param echo_map Property map to access the echo values of the input points (if any).
  */
  template <typename EchoMap>
  void generate_echo_based_attributes(Classifier& psc,
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

  /// \cond SKIP_IN_MANUAL
  void generate_echo_based_attributes(const Classifier&,
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
  /// \endcond


  /*!
    \brief Saves the current configuration in the file `filename`.

    This allows to easily save and recover a specific classification
    configuration, that is to say:

    - The smallest voxel size defined
    - The attributes and their respective weights
    - The classification types and the effect the attributes have on them

    The output file is written in an XML format that is readable by
    the `load()` method and the constructor that takes a file name
    as parameter.

    \param filename Name of the output file
    \param psc Classification object whose attributes and types must be saved
  */
  void save (const char* filename, Classifier& psc)
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

  
  /*!
    \brief Load a configuration from the file `filename`.

    All data structures, attributes and types specified in the input
    file `filename` are instanciated if possible (in particular,
    property maps needed should be provided).

    The input file is written in an XML format written by the `save()`
    method.

    \tparam VectorMap is a model of `ReadablePropertyMap` with value type `Vector_3<Kernel>`.
    \tparam ColorMap is a model of `ReadablePropertyMap` with value type `CGAL::Classification::RGB_Color`.
    \tparam EchoMap is a model of `ReadablePropertyMap` with value type `std::size_`.
    \param filename Name of the output file
    \param psc Classification object where to store attributes and types
    \param begin Iterator to the first input object
    \param end Past-the-end iterator
    \param point_map Property map to access the input points
    \param normal_map Property map to access the normal vectors of the input points (if any).
    \param color_map Property map to access the colors of the input points (if any).
    \param echo_map Property map to access the echo values of the input points (if any).
  */
  template<typename VectorMap = CGAL::Empty_property_map<RandomAccessIterator, typename Kernel::Vector_3>,
           typename ColorMap = CGAL::Empty_property_map<RandomAccessIterator, RGB_Color>,
           typename EchoMap  = CGAL::Empty_property_map<RandomAccessIterator, std::size_t> >
  bool load (const char* filename, Classifier& psc, 
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

    CGAL::Timer t;
    std::map<std::string, Attribute_handle> att_map;
    BOOST_FOREACH(boost::property_tree::ptree::value_type &v, tree.get_child("classification.attributes"))
      {
        std::string full_id = v.second.get<std::string>("id");

        std::vector<std::string> splitted_id;
        boost::split(splitted_id, full_id, boost::is_any_of("_"));
        std::string id = splitted_id[0];
        for (std::size_t i = 1; i < splitted_id.size() - 1; ++ i)
          id = id + "_" + splitted_id[i];
        std::size_t scale = std::atoi (splitted_id.back().c_str());

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
          {
            t.start();
            psc.add_attribute (Attribute_handle (new Elevation(begin, end, point_map,
                                                               *(m_scales[scale]->grid),
                                                               m_scales[scale]->grid_resolution(),
                                                               m_scales[scale]->radius_dtm())));
            t.stop();
          }
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
    std::cerr << "Elevation took " << t.time() << " second(s)" << std::endl;

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


  /*!
    \brief Writes a classification in a colored and labeled PLY format.

    The input point set is written in a PLY format with the addition
    of several PLY properties:

    - a property `label` to indicate which classification type is
    assigned to the point. The types are indexed from 0 to N (the
    correspondancy is given as comments in the PLY header).

    - 3 properties `red`, `green` and `blue` to associate each label
    to a color (this is useful to visualize the classification in a
    viewer that supports PLY colors)

    \param stream The output stream where to write the content
    \param begin Iterator to the first input object
    \param end Past-the-end iterator
    \param point_map Property map to access the input points
    \param psc The classification object to write from
    \param colors A set of colors to be used to represent the
    different classification types. If none is given, random colors
    are picked.
  */
  void write_ply (std::ostream& stream,
                  RandomAccessIterator begin, RandomAccessIterator end,
                  PointMap point_map, Classifier& psc,
                  std::vector<RGB_Color>* colors = NULL)
  {
    stream << "ply" << std::endl
           << "format ascii 1.0" << std::endl
           << "comment Generated by the CGAL library www.cgal.org" << std::endl
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
    stream << "comment label -1 is (unclassified)" << std::endl;
    for (std::size_t i = 0; i < psc.number_of_classification_types(); ++ i)
      {
        map_types.insert (std::make_pair (psc.get_classification_type(i), i));
        stream << "comment label " << i << " is " << psc.get_classification_type(i)->id() << std::endl;
      }
    map_types.insert (std::make_pair (Type_handle(), psc.number_of_classification_types()));
    
    stream << "end_header" << std::endl;

    while (colors->size() < psc.number_of_classification_types())
      {
        RGB_Color c = {{ (unsigned char)(64 + rand() % 128),
                         (unsigned char)(64 + rand() % 128),
                         (unsigned char)(64 + rand() % 128) }};
        colors->push_back (c);
      }

    std::size_t i = 0;
    for (RandomAccessIterator it = begin; it != end; ++ it)
      {
        Type_handle t = psc.classification_type_of(i);
        std::size_t idx = map_types[t];

        if (idx == psc.number_of_classification_types())
          stream << get(point_map, *it) << " 0 0 0 -1" << std::endl;
        else
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

private:
    /// \cond SKIP_IN_MANUAL
  template <typename Attribute_type>
  void generate_multiscale_attribute_variant_0 (Classifier& psc,
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
  void generate_multiscale_attribute_variant_1 (Classifier& psc,
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
  void generate_multiscale_attribute_variant_2 (Classifier& psc,
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
  void generate_multiscale_attribute_variant_3 (Classifier& psc,
                                                RandomAccessIterator begin,
                                                RandomAccessIterator end,
                                                PointMap point_map)
  {
    for (std::size_t i = 0; i < m_scales.size(); ++ i)
      {
        psc.add_attribute (Attribute_handle (new Attribute_type(begin, end, point_map,
                                                                *(m_scales[i]->grid),
                                                                m_scales[i]->grid_resolution(),
                                                                m_scales[i]->radius_dtm())));
        m_scales[i]->attributes.push_back (psc.get_attribute (psc.number_of_attributes() - 1));
      }
  }

  /// \endcond
};


} // namespace Classification

} // namespace CGAL


#endif // CGAL_CLASSIFICATION_HELPER_H
