// Copyright (c) 2017 GeometryFactory Sarl (France).
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
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_CLASSIFICATION_POINT_SET_FEATURE_GENERATOR_H
#define CGAL_CLASSIFICATION_POINT_SET_FEATURE_GENERATOR_H

#include <CGAL/license/Classification.h>

#include <CGAL/property_map.h>
#include <CGAL/Classification/Point_set_neighborhood.h>
#include <CGAL/Classification/Planimetric_grid.h>
#include <CGAL/Classification/Local_eigen_analysis.h>
#include <CGAL/Classification/Feature_base.h>
#include <CGAL/Classification/Feature/Distance_to_plane.h>
#include <CGAL/Classification/Feature/Echo_scatter.h>
#include <CGAL/Classification/Feature/Elevation.h>
#include <CGAL/Classification/Feature/Vertical_dispersion.h>
#include <CGAL/Classification/Feature/Verticality.h>
#include <CGAL/Classification/Feature/Eigenvalue.h>
#include <CGAL/Classification/Feature/Color_channel.h>

// Experimental feature, not used officially
#ifdef CGAL_CLASSIFICATION_USE_GRADIENT_OF_FEATURE
#include <CGAL/Classification/Feature/Gradient_of_feature.h>
#endif

#include <CGAL/Classification/Label.h>
#include <CGAL/Classification/internal/verbosity.h>
#include <CGAL/Classification/Feature_set.h>
#include <CGAL/bounding_box.h>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include <CGAL/Real_timer.h>
#include <CGAL/demangle.h>


namespace CGAL {

namespace Classification {

/*!
  \ingroup PkgClassificationPointSet

  \brief Generates a set of generic features for point set
  classification.

  This class takes care of computing all necessary data structures and
  of generating a set of generic features at multiple scales to
  increase the reliability of the classification.

  \tparam GeomTraits model of \cgal Kernel.
  \tparam PointRange model of `ConstRange`. Its iterator type is
  `RandomAccessIterator` and its value type is the key type of
  `PointMap`.
  \tparam PointMap model of `ReadablePropertyMap` whose key
  type is the value type of the iterator of `PointRange` and value type
  is `GeomTraits::Point_3`.
  \tparam ConcurrencyTag enables sequential versus parallel
  computation of `CGAL::Classification::Local_eigen_analysis`
  objects. Possible values are `Parallel_tag` (default value is %CGAL
  is linked with TBB) or `Sequential_tag` (default value otherwise).
  \tparam DiagonalizeTraits model of `DiagonalizeTraits` used for
  matrix diagonalization. It can be omitted: if Eigen 3 (or greater)
  is available and `CGAL_EIGEN3_ENABLED` is defined: in that case, an
  overload using `Eigen_diagonalize_traits` is provided.
*/
template <typename GeomTraits,
          typename PointRange,
          typename PointMap,
#if defined(DOXYGEN_RUNNING)
          typename ConcurrencyTag,
#elif defined(CGAL_LINKED_WITH_TBB)
          typename ConcurrencyTag = CGAL::Parallel_tag,
#else
          typename ConcurrencyTag = CGAL::Sequential_tag,
#endif
#if defined(DOXYGEN_RUNNING)
          typename DiagonalizeTraits>
#else
          typename DiagonalizeTraits = CGAL::Default_diagonalize_traits<float,3> >
#endif
class Point_set_feature_generator
{
  
public:
  typedef typename GeomTraits::Iso_cuboid_3             Iso_cuboid_3;

  /// \cond SKIP_IN_MANUAL
  typedef typename PointRange::const_iterator Iterator;
  typedef typename PointMap::value_type       Point;
  /// \endcond
  
  typedef Classification::Planimetric_grid
  <GeomTraits, PointRange, PointMap>                    Planimetric_grid;
  typedef Classification::Point_set_neighborhood
  <GeomTraits, PointRange, PointMap>                    Neighborhood;
  typedef Classification::Local_eigen_analysis           Local_eigen_analysis;

  /// \cond SKIP_IN_MANUAL
  typedef Classification::Feature_handle                 Feature_handle;
  typedef Classification::Label                          Label;
  typedef Classification::Label_handle                   Label_handle;

  typedef Classification::Feature::Distance_to_plane
  <PointRange, PointMap>                                 Distance_to_plane;
  typedef Classification::Feature::Elevation
  <GeomTraits, PointRange, PointMap>                    Elevation;
  typedef Classification::Feature::Vertical_dispersion
  <GeomTraits, PointRange, PointMap>                    Dispersion;
  typedef Classification::Feature::Verticality
  <GeomTraits>                                          Verticality;
  typedef Classification::Feature::Eigenvalue           Eigenvalue;
  
  typedef typename Neighborhood::K_neighbor_query       Neighbor_query;

#ifdef CGAL_CLASSIFICATION_USE_GRADIENT_OF_FEATURE
  typedef Classification::Feature::Gradient_of_feature
  <PointRange, PointMap, Neighbor_query>                Gradient_of_feature;
#endif
  
  typedef typename Classification::RGB_Color RGB_Color;
  /// \endcond
    
private:

  struct Scale
  {
    Neighborhood* neighborhood;
    Planimetric_grid* grid;
    Local_eigen_analysis* eigen;
    float voxel_size;
    
    Scale (const PointRange& input, PointMap point_map,
           const Iso_cuboid_3& bbox, float voxel_size,
           Planimetric_grid* lower_grid = NULL)
      : voxel_size (voxel_size)
    {
      CGAL::Real_timer t;
      t.start();
      if (lower_grid == NULL)
        neighborhood = new Neighborhood (input, point_map);
      else
        neighborhood = new Neighborhood (input, point_map, voxel_size);
      t.stop();
      
      if (voxel_size < 0.)
        CGAL_CLASSIFICATION_CERR << "Neighborhood computed in " << t.time() << " second(s)" << std::endl;
      else
        CGAL_CLASSIFICATION_CERR << "Neighborhood with voxel size " << voxel_size
                                 << " computed in " << t.time() << " second(s)" << std::endl;
      t.reset();
      t.start();
      
      eigen = new Local_eigen_analysis
        (Local_eigen_analysis::create_from_point_set
         (input, point_map, neighborhood->k_neighbor_query(12), ConcurrencyTag(), DiagonalizeTraits()));
      
      float range = eigen->mean_range();
      if (this->voxel_size < 0)
        this->voxel_size = range;
      t.stop();
      CGAL_CLASSIFICATION_CERR << "Eigen values computed in " << t.time() << " second(s)" << std::endl;
      CGAL_CLASSIFICATION_CERR << "Range = " << range << std::endl;
      t.reset();
      t.start();

      if (lower_grid == NULL)
        grid = new Planimetric_grid (input, point_map, bbox, this->voxel_size);
      else
        grid = new Planimetric_grid(lower_grid);
      t.stop();
      CGAL_CLASSIFICATION_CERR << "Planimetric grid computed in " << t.time() << " second(s)" << std::endl;
      t.reset();
    }
    ~Scale()
    {
      if (neighborhood != NULL)
        delete neighborhood;
      if (grid != NULL)
        delete grid;
      delete eigen;
    }

    void reduce_memory_footprint(bool delete_neighborhood)
    {
      delete grid;
      grid = NULL;
      if (delete_neighborhood)
      {
        delete neighborhood;
        neighborhood = NULL;
      }
    }

    float grid_resolution() const { return voxel_size; }
    float radius_neighbors() const { return voxel_size * 5; }
    float radius_dtm() const { return voxel_size * 100; }
    
  };

  Iso_cuboid_3 m_bbox;
  std::vector<Scale*> m_scales;

  const PointRange& m_input;
  PointMap m_point_map;
  
public:
  
  /// \name Constructor
  /// @{
  
  /*!
    \brief Initializes a feature generator from an input range.

    If not provided by the user, The size of the smallest scale is
    automatically estimated using a method equivalent to
    `CGAL::compute_average_spacing()` using 6 neighbors. The data
    structures needed (`Neighborhood`, `Planimetric_grid` and
    `Local_eigen_analysis`) are computed at `nb_scales` recursively
    larger scales. 

    \param input point range.
    \param point_map property map to access the input points.
    \param nb_scales number of scales to compute.
    \param voxel_size smallest scale used as a voxel size for the
    planimetric grid (if the default value -1 is used, its value is
    automatically estimated).
  */
  Point_set_feature_generator(const PointRange& input,
                              PointMap point_map,
                              std::size_t nb_scales,
                              float voxel_size = -1.f)
    : m_input (input), m_point_map (point_map)
  {
    m_bbox = CGAL::bounding_box
      (boost::make_transform_iterator (m_input.begin(), CGAL::Property_map_to_unary_function<PointMap>(m_point_map)),
       boost::make_transform_iterator (m_input.end(), CGAL::Property_map_to_unary_function<PointMap>(m_point_map)));

    CGAL::Real_timer t; t.start();
    
    m_scales.reserve (nb_scales);
    
    m_scales.push_back (new Scale (m_input, m_point_map, m_bbox, voxel_size));

    if (voxel_size == -1.f)
      voxel_size = m_scales[0]->grid_resolution();
    
    for (std::size_t i = 1; i < nb_scales; ++ i)
    {
      voxel_size *= 2;
      m_scales.push_back (new Scale (m_input, m_point_map, m_bbox, voxel_size, m_scales[i-1]->grid));
    }
    t.stop();
    CGAL_CLASSIFICATION_CERR << "Scales computed in " << t.time() << " second(s)" << std::endl;
    t.reset();
  }

  /// @}
  
  /// \cond SKIP_IN_MANUAL
  
#ifndef CGAL_NO_DEPRECATED_CODE
  // deprecated
  template <typename VectorMap = Default,
            typename ColorMap = Default,
            typename EchoMap = Default>
  CGAL_DEPRECATED_MSG("you are using a deprecated constructor of CGAL::Classification::Point_set_feature_generator, please update your code")
  Point_set_feature_generator(Feature_set& features,
                              const PointRange& input,
                              PointMap point_map,
                              std::size_t nb_scales,
                              VectorMap normal_map = VectorMap(),
                              ColorMap color_map = ColorMap(),
                              EchoMap echo_map = EchoMap(),
                              float voxel_size = -1.f)
    : m_input (input), m_point_map (point_map)
  {
    m_bbox = CGAL::bounding_box
      (boost::make_transform_iterator (m_input.begin(), CGAL::Property_map_to_unary_function<PointMap>(m_point_map)),
       boost::make_transform_iterator (m_input.end(), CGAL::Property_map_to_unary_function<PointMap>(m_point_map)));

    CGAL::Real_timer t; t.start();
    
    m_scales.reserve (nb_scales);
    
    m_scales.push_back (new Scale (m_input, m_point_map, m_bbox, voxel_size));

    if (voxel_size == -1.f)
      voxel_size = m_scales[0]->grid_resolution();
    
    for (std::size_t i = 1; i < nb_scales; ++ i)
    {
      voxel_size *= 2;
      m_scales.push_back (new Scale (m_input, m_point_map, m_bbox, voxel_size, m_scales[i-1]->grid));
    }
    t.stop();
    CGAL_CLASSIFICATION_CERR << "Scales computed in " << t.time() << " second(s)" << std::endl;
    t.reset();

    typedef typename Default::Get<VectorMap, typename GeomTraits::Vector_3 >::type
      Vmap;
    typedef typename Default::Get<ColorMap, RGB_Color >::type
      Cmap;
    typedef typename Default::Get<EchoMap, std::size_t >::type
      Emap;

    generate_point_based_features (features);
    generate_normal_based_features (features, get_parameter<Vmap>(normal_map));
    generate_color_based_features (features, get_parameter<Cmap>(color_map));
    generate_echo_based_features (features, get_parameter<Emap>(echo_map));
  }

  // Functions to remove when deprecated constructor is removed
  void generate_normal_based_features(const CGAL::Constant_property_map<Iterator, typename GeomTraits::Vector_3>&) { }
  void generate_color_based_features(const CGAL::Constant_property_map<Iterator, RGB_Color>&) { }
  void generate_echo_based_features(const CGAL::Constant_property_map<Iterator, std::size_t>&) { }
#endif
  
  virtual ~Point_set_feature_generator()
  {
    clear();
  }

  void reduce_memory_footprint()
  {
    for (std::size_t i = 0; i < m_scales.size(); ++ i)
    {
      m_scales[i]->reduce_memory_footprint(i > 0);
    }
  }
  /// \endcond

  /// \name Feature Generation
  /// @{


  /*!
    \brief Generate geometric features based on point position information.

    At each scale, the following features are generated:

    - `CGAL::Classification::Feature::Eigenvalue` with indices 0, 1 and 2
    - `CGAL::Classification::Feature::Distance_to_plane`
    - `CGAL::Classification::Feature::Elevation`
    - `CGAL::Classification::Feature::Vertical_dispersion`
    - The version of `CGAL::Classification::Feature::Verticality` based on eigenvalues

    \param features the feature set where the features are instantiated.
   */
  void generate_point_based_features (Feature_set& features)
  {
    for (int j = 0; j < 3; ++ j)
      for (std::size_t i = 0; i < m_scales.size(); ++ i)
        features.add_with_scale_id<Eigenvalue> (i, m_input, eigen(i), (unsigned int)(j));
    for (std::size_t i = 0; i < m_scales.size(); ++ i)
      features.add_with_scale_id<Distance_to_plane> (i, m_input, m_point_map, eigen(i));
    for (std::size_t i = 0; i < m_scales.size(); ++ i)
      features.add_with_scale_id<Dispersion> (i, m_input, m_point_map, grid(i), radius_neighbors(i));
    for (std::size_t i = 0; i < m_scales.size(); ++ i)
      features.add_with_scale_id<Elevation> (i, m_input, m_point_map, grid(i), radius_dtm(i));
    for (std::size_t i = 0; i < m_scales.size(); ++ i)
      features.add_with_scale_id<Verticality> (i, m_input, eigen(i));
  }

  /*!
    \brief Generate geometric features based on normal vector information.

    Generates the version of `CGAL::Classification::Feature::Verticality` based on normal vectors.

    \tparam VectorMap model of `ReadablePropertyMap` whose key type is
    the value type of the iterator of `PointRange` and value type is
    `GeomTraits::Vector_3`.

    \param features the feature set where the features are instantiated.
    \param normal_map property map to access the normal vectors of the input points (if any).

   */
  template <typename VectorMap>
  void generate_normal_based_features(Feature_set& features, const VectorMap& normal_map)
  {
    features.add<Verticality> (m_input, normal_map);
  }

  /*!
    \brief Generate geometric features based on point color information.

    Generates `CGAL::Classification::Feature::Color_channel` with
    channels `HUE`, `SATURATION` and `VALUE`.

    \tparam ColorMap model of `ReadablePropertyMap`  whose key type is
    the value type of the iterator of `PointRange` and value type is
    `CGAL::Classification::RGB_Color`.

    \param features the feature set where the features are instantiated.
    \param color_map property map to access the colors of the input points (if any).
   */
  template <typename ColorMap>
  void generate_color_based_features(Feature_set& features, const ColorMap& color_map)
  {
    typedef Feature::Color_channel<GeomTraits, PointRange, ColorMap> Color_channel;
    for (std::size_t i = 0; i < 3; ++ i)
      features.add<Color_channel> (m_input, color_map, typename Color_channel::Channel(i));
  }

  /*!
    \brief Generate geometric features based on echo information.

    At each scale, generates `CGAL::Classification::Feature::Echo_scatter`.

    \tparam EchoMap model of `ReadablePropertyMap` whose key type is
    the value type of the iterator of `PointRange` and value type is
    `std::size_t`.

    \param features the feature set where the features are instantiated.
    \param echo_map property map to access the echo values of the input points (if any).
   */
  template <typename EchoMap>
  void generate_echo_based_features(Feature_set& features, const EchoMap& echo_map)
  {
    typedef Feature::Echo_scatter<GeomTraits, PointRange, PointMap, EchoMap> Echo_scatter;
    for (std::size_t i = 0; i < m_scales.size(); ++ i)
      features.add_with_scale_id<Echo_scatter> (i, m_input, echo_map, grid(i), radius_neighbors(i));
  }

  /// @}

  /// \name Data Structures Access
  /// @{

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

  /// @}

  /// \name Parameters
  /// @{
  
  /*!
    \brief Returns the number of scales that were computed.
  */
  std::size_t number_of_scales() const { return m_scales.size(); }
  
  /*!
    \brief Returns the grid resolution at scale `scale`. This
    resolution is the length and width of a cell of the
    `Planimetric_grid` defined at this scale.
  */
  float grid_resolution(std::size_t scale = 0) const { return m_scales[scale]->grid_resolution(); }
  /*!

    \brief Returns the radius used for neighborhood queries at scale
    `scale`. This radius is the smallest radius that is relevant from
    a geometric point of view at this scale (that is to say that
    encloses a few cells of `Planimetric_grid`).
  */
  float radius_neighbors(std::size_t scale = 0) const { return m_scales[scale]->radius_neighbors(); }
  /*!
    \brief Returns the radius used for digital terrain modeling at
    scale `scale`. This radius represents the minimum size of a
    building at this scale.
  */
  float radius_dtm(std::size_t scale = 0) const { return m_scales[scale]->radius_dtm(); }

  /// @}

    
private:

  void clear()
  {
    for (std::size_t i = 0; i < m_scales.size(); ++ i)
      delete m_scales[i];
    m_scales.clear();
  }

#ifdef CGAL_CLASSIFICATION_USE_GRADIENT_OF_FEATURE
  void generate_gradient_features(Feature_set& features)
  {
    std::size_t size = features->size();

    for (std::size_t i = 0; i < size; ++ i)
    {
      for (int j = m_scales.size() - 1; j >= 0; -- j)
      {
        std::ostringstream oss;
        oss << "_" << j;
        if ((*features)[i]->name().find (oss.str()))
        {
          const Neighbor_query& neighbor_query = neighborhood(std::size_t(j)).k_neighbor_query(6);
          features->template add<Gradient_of_feature> (m_input, m_point_map, (*features)[i], neighbor_query);
          break;
        }
      }
    }
  }
#endif

  template <typename T>
  const T& get_parameter (const T& t)
  {
    return t;
  }

  template <typename T>
  Constant_property_map<Iterator, T>
  get_parameter (const Default&)
  {
    return Constant_property_map<Iterator, T>();
  }

};


} // namespace Classification
  
} // namespace CGAL



#endif // CGAL_CLASSIFICATION_POINT_SET_FEATURE_GENERATOR_H
