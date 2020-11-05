// Copyright (c) 2017 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
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
#include <CGAL/Classification/Feature/Height_below.h>
#include <CGAL/Classification/Feature/Height_above.h>
#include <CGAL/Classification/Feature/Vertical_range.h>

#include <CGAL/Classification/internal/verbosity.h>
#include <CGAL/Classification/Feature_set.h>
#include <CGAL/bounding_box.h>

#include <CGAL/Real_timer.h>


namespace CGAL {

namespace Classification {

/*!
  \ingroup PkgClassificationPointSet

  \brief Generates a set of generic features for point set
  classification.

  This class takes care of computing and storing all necessary data
  structures and of generating a set of generic features at multiple
  scales to increase the reliability of the classification.

  \warning The generated features use data structures that are stored
  inside the generator. For this reason, the generator should be
  instantiated _within the same scope_ as the feature set and should
  not be deleted before the feature set.

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
#else
          typename ConcurrencyTag = CGAL::Parallel_if_available_tag,
#endif
#if defined(DOXYGEN_RUNNING)
          typename DiagonalizeTraits>
#else
          typename DiagonalizeTraits = CGAL::Default_diagonalize_traits<float,3> >
#endif
class Point_set_feature_generator
{

public:
  using Iso_cuboid_3 = typename GeomTraits::Iso_cuboid_3;

  /// \cond SKIP_IN_MANUAL
  using Iterator = typename PointRange::const_iterator;
  using Point = typename PointMap::value_type;
  /// \endcond

  using Planimetric_grid = Classification::Planimetric_grid<GeomTraits, PointRange, PointMap>;
  using Neighborhood = Classification::Point_set_neighborhood<GeomTraits, PointRange, PointMap>;
  using Local_eigen_analysis = Classification::Local_eigen_analysis;

  /// \cond SKIP_IN_MANUAL
  using Feature_handle = Classification::Feature_handle;
  using Distance_to_plane = Classification::Feature::Distance_to_plane<PointRange, PointMap>;
  using Elevation = Classification::Feature::Elevation<GeomTraits, PointRange, PointMap>;
  using Height_below = Classification::Feature::Height_below<GeomTraits, PointRange, PointMap>;
  using Height_above = Classification::Feature::Height_above<GeomTraits, PointRange, PointMap>;
  using Vertical_range = Classification::Feature::Vertical_range<GeomTraits, PointRange, PointMap>;
  using Dispersion = Classification::Feature::Vertical_dispersion<GeomTraits, PointRange, PointMap>;
  using Verticality = Classification::Feature::Verticality<GeomTraits>;
  using Eigenvalue = Classification::Feature::Eigenvalue;
  using Neighbor_query = typename Neighborhood::K_neighbor_query;
  /// \endcond

private:

  struct Scale
  {
    std::unique_ptr<Neighborhood> neighborhood;
    std::unique_ptr<Planimetric_grid> grid;
    std::unique_ptr<Local_eigen_analysis> eigen;
    float voxel_size;

    Scale (const PointRange& input, PointMap point_map,
           const Iso_cuboid_3& bbox, float voxel_size,
           const std::unique_ptr<Planimetric_grid>& lower_grid
           = std::unique_ptr<Planimetric_grid>())
      : voxel_size (voxel_size)
    {
      CGAL::Real_timer t;
      t.start();
      if (!lower_grid)
        neighborhood = std::make_unique<Neighborhood> (input, point_map, ConcurrencyTag());
      else
        neighborhood = std::make_unique<Neighborhood> (input, point_map, voxel_size, ConcurrencyTag());
      t.stop();

      if (!lower_grid)
        CGAL_CLASSIFICATION_CERR << "Neighborhood computed in " << t.time() << " second(s)" << std::endl;
      else
        CGAL_CLASSIFICATION_CERR << "Neighborhood with voxel size " << voxel_size
                                 << " computed in " << t.time() << " second(s)" << std::endl;
      t.reset();
      t.start();

      eigen = std::make_unique<Local_eigen_analysis>
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

      if (!lower_grid)
        grid = std::make_unique<Planimetric_grid> (input, point_map, bbox, this->voxel_size);
      else
        grid = std::make_unique<Planimetric_grid>(lower_grid.get());
      t.stop();
      CGAL_CLASSIFICATION_CERR << "Planimetric grid computed in " << t.time() << " second(s)" << std::endl;
      t.reset();
    }

    float grid_resolution() const { return voxel_size; }
    float radius_neighbors() const { return voxel_size * 3; }
    float radius_dtm() const { return voxel_size * 10; }

  };

  Iso_cuboid_3 m_bbox;
  std::vector<std::unique_ptr<Scale> > m_scales;

  const PointRange& m_input;
  PointMap m_point_map;

public:

  /// \name Constructor
  /// @{

  /*!
    \brief initializes a feature generator from an input range.

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
      (CGAL::make_transform_iterator_from_property_map (m_input.begin(), m_point_map),
       CGAL::make_transform_iterator_from_property_map (m_input.end(), m_point_map));

    CGAL::Real_timer t; t.start();

    m_scales.reserve (nb_scales);

    m_scales.emplace_back (std::make_unique<Scale> (m_input, m_point_map, m_bbox, voxel_size));

    if (voxel_size == -1.f)
      voxel_size = m_scales[0]->grid_resolution();

    for (std::size_t i = 1; i < nb_scales; ++ i)
    {
      voxel_size *= 2;
      m_scales.push_back (std::make_unique<Scale> (m_input, m_point_map, m_bbox, voxel_size, m_scales[i-1]->grid));
    }
    t.stop();
    CGAL_CLASSIFICATION_CERR << "Scales computed in " << t.time() << " second(s)" << std::endl;
    t.reset();
  }

  /// @}

  /// \name Feature Generation
  /// @{


  /*!
    \brief generates geometric features based on point position information.

    This is a meta-function that calls the following functions:

    - `generate_eigen_features()`
    - `generate_dispersion_features()`
    - `generate_elevation_features()`
    - The version of `generate_normal_based_features()` without a normal map

    \param features the feature set where the features are instantiated.
   */
  void generate_point_based_features (Feature_set& features)
  {
    generate_eigen_features (features);
    generate_dispersion_features (features);
    generate_elevation_features (features);
    generate_normal_based_features (features);
  }

  /*!
    \brief generates geometric eigen features.

    At each scale, features
    `CGAL::Classification::Feature::Eigenvalue` with indices 0, 1 and
    2 are generated.

    \param features the feature set where the features are instantiated.
   */
  void generate_eigen_features (Feature_set& features)
  {
    for (int j = 0; j < 3; ++ j)
      for (std::size_t i = 0; i < m_scales.size(); ++ i)
        features.add_with_scale_id<Eigenvalue> (i, m_input, eigen(i), (unsigned int)(j));
  }

  /*!
    \brief generates geometric features based on local dispersion information.

    At each scale, the following features are generated:

    - `CGAL::Classification::Feature::Distance_to_plane`
    - `CGAL::Classification::Feature::Vertical_dispersion`

    \param features the feature set where the features are instantiated.
   */
  void generate_dispersion_features (Feature_set& features)
  {
    for (std::size_t i = 0; i < m_scales.size(); ++ i)
      features.add_with_scale_id<Distance_to_plane> (i, m_input, m_point_map, eigen(i));
    for (std::size_t i = 0; i < m_scales.size(); ++ i)
      features.add_with_scale_id<Dispersion> (i, m_input, m_point_map, grid(i), radius_neighbors(i));
  }

  /*!
    \brief generates geometric features based on elevation information.

    At each scale, the following features are generated:

    - `CGAL::Classification::Feature::Elevation`
    - `CGAL::Classification::Feature::Height_above`
    - `CGAL::Classification::Feature::Height_below`
    - `CGAL::Classification::Feature::Vertical_range`

    \param features the feature set where the features are instantiated.
   */
  void generate_elevation_features (Feature_set& features)
  {
    for (std::size_t i = 0; i < m_scales.size(); ++ i)
      features.add_with_scale_id<Elevation> (i, m_input, m_point_map, grid(i), radius_dtm(i));
    for (std::size_t i = 0; i < m_scales.size(); ++ i)
      features.add_with_scale_id<Height_below> (i, m_input, m_point_map, grid(i));
    for (std::size_t i = 0; i < m_scales.size(); ++ i)
      features.add_with_scale_id<Height_above> (i, m_input, m_point_map, grid(i));
    for (std::size_t i = 0; i < m_scales.size(); ++ i)
      features.add_with_scale_id<Vertical_range> (i, m_input, m_point_map, grid(i));
  }


  /*!
    \brief generates geometric features based on normal analysis information.

    At each scale, the version of
    `CGAL::Classification::Feature::Verticality` based on eigenvalue
    is generated.

    \param features the feature set where the features are instantiated.
   */
  void generate_normal_based_features (Feature_set& features)
  {
    for (std::size_t i = 0; i < m_scales.size(); ++ i)
      features.add_with_scale_id<Verticality> (i, m_input, eigen(i));
  }

  /*!
    \brief generates geometric features based on normal vector information.

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
    \brief generates geometric features based on point color information.

    Generates `CGAL::Classification::Feature::Color_channel` with
    channels `HUE`, `SATURATION` and `VALUE`.

    \tparam ColorMap model of `ReadablePropertyMap`  whose key type is
    the value type of the iterator of `PointRange` and value type is
    `CGAL::Color`.

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
    \brief generates geometric features based on echo information.

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
    \brief returns the bounding box of the input point set.
  */
  const Iso_cuboid_3& bbox() const { return m_bbox; }
  /*!
    \brief returns the neighborhood structure at scale `scale`.
  */
  const Neighborhood& neighborhood(std::size_t scale = 0) const { return (*m_scales[scale]->neighborhood); }
  /*!
    \brief returns the planimetric grid structure at scale `scale`.
  */
  const Planimetric_grid& grid(std::size_t scale = 0) const { return *(m_scales[scale]->grid); }
  /*!
    \brief returns the local eigen analysis structure at scale `scale`.
  */
  const Local_eigen_analysis& eigen(std::size_t scale = 0) const { return *(m_scales[scale]->eigen); }

  /// @}

  /// \name Parameters
  /// @{

  /*!
    \brief returns the number of scales that were computed.
  */
  std::size_t number_of_scales() const { return m_scales.size(); }

  /*!
    \brief returns the grid resolution at scale `scale`. This
    resolution is the length and width of a cell of the
    `Planimetric_grid` defined at this scale.
  */
  float grid_resolution(std::size_t scale = 0) const { return m_scales[scale]->grid_resolution(); }
  /*!

    \brief returns the radius used for neighborhood queries at scale
    `scale`. This radius is the smallest radius that is relevant from
    a geometric point of view at this scale (that is to say that
    encloses a few cells of `Planimetric_grid`).
  */
  float radius_neighbors(std::size_t scale = 0) const { return m_scales[scale]->radius_neighbors(); }
  /*!
    \brief returns the radius used for digital terrain modeling at
    scale `scale`. This radius represents the minimum size of a
    building at this scale.
  */
  float radius_dtm(std::size_t scale = 0) const { return m_scales[scale]->radius_dtm(); }

  /// @}

};


} // namespace Classification

} // namespace CGAL



#endif // CGAL_CLASSIFICATION_POINT_SET_FEATURE_GENERATOR_H
