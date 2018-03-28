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
#include <CGAL/Classification/Feature/Hsv.h>
#include <CGAL/Classification/Feature/Distance_to_plane.h>
#include <CGAL/Classification/Feature/Echo_scatter.h>
#include <CGAL/Classification/Feature/Elevation.h>
#include <CGAL/Classification/Feature/Vertical_dispersion.h>
#include <CGAL/Classification/Feature/Verticality.h>
#include <CGAL/Classification/Feature/Eigen.h>

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

#ifdef CGAL_LINKED_WITH_TBB
#include <tbb/task_group.h>
#include <tbb/mutex.h>
#endif // CGAL_LINKED_WITH_TBB

namespace CGAL {

namespace Classification {

/*!
  \ingroup PkgClassificationDataStructures

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
  algorithm. Possible values are `Parallel_tag` (default value is %CGAL
  is linked with TBB) or `Sequential_tag` (default value otherwise).
  \tparam DiagonalizeTraits model of `DiagonalizeTraits` used for
  matrix diagonalization. It can be omitted: if Eigen 3 (or greater)
  is available and `CGAL_EIGEN3_ENABLED` is defined then an overload
  using `Eigen_diagonalize_traits` is provided. Otherwise, the
  internal implementation `Diagonalize_traits` is used.

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
  
  typedef Classification::Feature::Anisotropy            Anisotropy;
  typedef Classification::Feature::Distance_to_plane
  <PointRange, PointMap>                                 Distance_to_plane;
  typedef Classification::Feature::Eigentropy            Eigentropy;
  typedef Classification::Feature::Elevation
  <GeomTraits, PointRange, PointMap>                    Elevation;
  typedef Classification::Feature::Linearity             Linearity;
  typedef Classification::Feature::Omnivariance          Omnivariance;
  typedef Classification::Feature::Planarity             Planarity;
  typedef Classification::Feature::Sphericity            Sphericity;
  typedef Classification::Feature::Sum_eigenvalues       Sum_eigen;
  typedef Classification::Feature::Surface_variation     Surface_variation;
  typedef Classification::Feature::Vertical_dispersion
  <GeomTraits, PointRange, PointMap>                    Dispersion;
  typedef Classification::Feature::Verticality
  <GeomTraits>                                          Verticality;
  
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
      if (voxel_size < 0.)
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
        (input, point_map, neighborhood->k_neighbor_query(6), ConcurrencyTag(), DiagonalizeTraits());
      
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

#ifdef CGAL_LINKED_WITH_TBB
  tbb::task_group* m_tasks;
#endif
  
  struct Feature_adder
  {
    mutable Point_set_feature_generator* generator;
    std::size_t scale;

    Feature_adder (Point_set_feature_generator* generator, std::size_t scale)
      : generator (generator), scale (scale) { }

    virtual ~Feature_adder() { }

    virtual void operator()() const = 0;
  };
  friend Feature_adder;

  const PointRange& m_input;
  PointMap m_point_map;
  std::vector<Feature_adder*> m_adders;
  Feature_set* m_features;
  
public:

  
  /*!
    \brief Generates all possible features from an input range.

    The size of the smallest scale is automatically estimated using a
    method equivalent to `CGAL::compute_average_spacing()` using 6
    neighbors. The data structures needed (`Neighborhood`,
    `Planimetric_grid` and `Local_eigen_analysis`) are computed at
    `nb_scales` recursively larger scales. At each scale, the
    following features are generated:

    - `CGAL::Classification::Feature::Anisotropy`
    - `CGAL::Classification::Feature::Distance_to_plane`
    - `CGAL::Classification::Feature::Eigentropy`
    - `CGAL::Classification::Feature::Elevation`
    - `CGAL::Classification::Feature::Linearity`
    - `CGAL::Classification::Feature::Omnivariance`
    - `CGAL::Classification::Feature::Planarity`
    - `CGAL::Classification::Feature::Sphericity`
    - `CGAL::Classification::Feature::Sum_eigenvalues`
    - `CGAL::Classification::Feature::Surface_variation`
    - `CGAL::Classification::Feature::Vertical_dispersion`
    - The version of `CGAL::Classification::Feature::Verticality` based on eigenvalues

    If normal vectors are provided (if `VectorMap` is different from
    `CGAL::Default`), the following feature is generated at each
    scale:

    - The version of `CGAL::Classification::Feature::Verticality` based on normal vectors

    If colors are provided (if `ColorMap` is different from
    `CGAL::Default`), the following features are generated at each
    scale:

    - 9 features `CGAL::Classification::Feature::Hsv` on
    channel 0 (hue) with mean ranging from 0° to 360° and standard
    deviation of 22.5.

    - 5 features `CGAL::Classification::Feature::Hsv` on
    channel 1 (saturation) with mean ranging from 0 to 100 and standard
    deviation of 12.5.

    - 5 features `CGAL::Classification::Feature::Hsv` on channel 2
    (value) with mean ranging from 0 to 100 and standard deviation
    of 12.5.

    If echo numbers are provided (if `EchoMap` is different from
    `CGAL::Default`), the following feature is computed at each
    scale:

    - `CGAL::Classification::Feature::Echo_scatter`

    \tparam VectorMap model of `ReadablePropertyMap` whose key type is
    the value type of the iterator of `PointRange` and value type is
    `GeomTraits::Vector_3`.
    \tparam ColorMap model of `ReadablePropertyMap`  whose key type is
    the value type of the iterator of `PointRange` and value type is
    `CGAL::Classification::RGB_Color`.
    \tparam EchoMap model of `ReadablePropertyMap` whose key type is
    the value type of the iterator of `PointRange` and value type is
    `std::size_t`.
    \param features the feature set where the features are instantiated.
    \param input point range.
    \param point_map property map to access the input points.
    \param nb_scales number of scales to compute.
    \param normal_map property map to access the normal vectors of the input points (if any).
    \param color_map property map to access the colors of the input points (if any).
    \param echo_map property map to access the echo values of the input points (if any).
  */
  template <typename VectorMap = Default,
            typename ColorMap = Default,
            typename EchoMap = Default>
  Point_set_feature_generator(Feature_set& features,
                              const PointRange& input,
                              PointMap point_map,
                              std::size_t nb_scales,
                              VectorMap normal_map = VectorMap(),
                              ColorMap color_map = ColorMap(),
                              EchoMap echo_map = EchoMap()
#ifndef DOXYGEN_RUNNING
                              , float voxel_size = -1.f // Undocumented way of changing base voxel size
#endif
                              )
    : m_input (input), m_point_map (point_map), m_features (&features)
  {
    m_bbox = CGAL::bounding_box
      (boost::make_transform_iterator (m_input.begin(), CGAL::Property_map_to_unary_function<PointMap>(m_point_map)),
       boost::make_transform_iterator (m_input.end(), CGAL::Property_map_to_unary_function<PointMap>(m_point_map)));

    typedef typename Default::Get<VectorMap, typename GeomTraits::Vector_3 >::type
      Vmap;
    typedef typename Default::Get<ColorMap, RGB_Color >::type
      Cmap;
    typedef typename Default::Get<EchoMap, std::size_t >::type
      Emap;

    generate_features_impl (nb_scales, voxel_size,
                            get_parameter<Vmap>(normal_map),
                            get_parameter<Cmap>(color_map),
                            get_parameter<Emap>(echo_map));

    m_features->sort_features_by_name();
  }

  
  /// \cond SKIP_IN_MANUAL
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


    
private:

  void clear()
  {
    for (std::size_t i = 0; i < m_scales.size(); ++ i)
      delete m_scales[i];
    m_scales.clear();
  }

  void generate_point_based_features ()
  {

    generate_multiscale_feature_variant_0<Anisotropy> ();
    generate_multiscale_feature_variant_0<Eigentropy> ();
    generate_multiscale_feature_variant_0<Linearity> ();
    generate_multiscale_feature_variant_0<Omnivariance> ();
    generate_multiscale_feature_variant_0<Planarity> ();
    generate_multiscale_feature_variant_0<Sphericity> ();
    generate_multiscale_feature_variant_0<Sum_eigen> ();
    generate_multiscale_feature_variant_0<Surface_variation> ();

    generate_multiscale_feature_variant_1<Distance_to_plane> ();
    generate_multiscale_feature_variant_2<Dispersion> ();
    generate_multiscale_feature_variant_3<Elevation> ();
  }

  template <typename FeatureAdder>
  void launch_feature_computation (FeatureAdder* adder)
  {
    m_adders.push_back (adder);
    
#ifndef CGAL_LINKED_WITH_TBB
    CGAL_static_assertion_msg (!(boost::is_convertible<ConcurrencyTag, Parallel_tag>::value),
                               "Parallel_tag is enabled but TBB is unavailable.");
#else
    if (boost::is_convertible<ConcurrencyTag,Parallel_tag>::value)
    {
      m_tasks->run (*adder);
    }
    else
#endif
    {
      (*adder)();
    }
  }

  template <typename VectorMap>
  struct Feature_adder_verticality : public Feature_adder
  {
    using Feature_adder::generator;
    using Feature_adder::scale;
    VectorMap normal_map;

    Feature_adder_verticality (Point_set_feature_generator* generator, VectorMap normal_map, std::size_t scale)
      : Feature_adder (generator, scale), normal_map (normal_map) { }
    
    void operator() () const
    {
      Feature_handle fh = generator->m_features->template add<Verticality> (generator->m_input, normal_map);
      std::ostringstream oss;
      oss << fh->name() << "_" << scale;
      fh->set_name (oss.str());
    }
  };

  template <typename VectorMap>
  void generate_normal_based_features(VectorMap normal_map)
  {
    launch_feature_computation (new Feature_adder_verticality<VectorMap> (this, normal_map, 0));
  }

  void generate_normal_based_features(const CGAL::Default_property_map<Iterator, typename GeomTraits::Vector_3>&)
  {
    generate_multiscale_feature_variant_0<Verticality> ();
  }

  template <typename ColorMap>
  struct Feature_adder_color : public Feature_adder
  {
    typedef Classification::Feature::Hsv<GeomTraits, PointRange, ColorMap> Hsv;
    
    using Feature_adder::generator;
    using Feature_adder::scale;
    ColorMap color_map;
    std::size_t channel;
    float mean;
    float sd;

    Feature_adder_color (Point_set_feature_generator* generator, ColorMap color_map, std::size_t scale,
                         std::size_t channel, float mean, float sd)
      : Feature_adder (generator, scale), color_map (color_map),
        channel (channel), mean (mean), sd (sd) { }
    
    void operator() () const
    {
      Feature_handle fh = generator->m_features->template add<Hsv> (generator->m_input, color_map,
                                                                    (typename Hsv::Channel)(channel), mean, sd);
      std::ostringstream oss;
      oss << fh->name() << "_" << scale;
      fh->set_name (oss.str());
    }
  };

  template <typename ColorMap>
  void generate_color_based_features(ColorMap color_map)
  {
    for (std::size_t i = 0; i <= 8; ++ i)
      launch_feature_computation (new Feature_adder_color<ColorMap> (this, color_map, 0,
                                                                     0, 45.f * float(i), 22.5f));

    for (std::size_t i = 0; i <= 4; ++ i)
      launch_feature_computation (new Feature_adder_color<ColorMap> (this, color_map, 0,
                                                                     1, 25.f * float(i), 12.5f));
    
    for (std::size_t i = 0; i <= 4; ++ i)
      launch_feature_computation (new Feature_adder_color<ColorMap> (this, color_map, 0,
                                                                     2, 25.f * float(i), 12.5f));
  }

  void generate_color_based_features(const CGAL::Default_property_map<Iterator, RGB_Color>&)
  {
  }


  template <typename EchoMap>
  struct Feature_adder_echo : public Feature_adder
  {
    typedef Classification::Feature::Echo_scatter<GeomTraits, PointRange, PointMap, EchoMap> Echo_scatter;
    
    using Feature_adder::generator;
    using Feature_adder::scale;
    EchoMap echo_map;

    Feature_adder_echo (Point_set_feature_generator* generator, EchoMap echo_map, std::size_t scale)
      : Feature_adder (generator, scale), echo_map (echo_map) { }
    
    void operator() () const
    {
      Feature_handle fh = generator->m_features->template add<Echo_scatter> (generator->m_input,
                                                                             echo_map,
                                                                             generator->grid(scale),
                                                                             generator->radius_neighbors(scale));
      std::ostringstream oss;
      oss << fh->name() << "_" << scale;
      fh->set_name (oss.str());
    }
  };


  template <typename EchoMap>
  void generate_echo_based_features(EchoMap echo_map)
  {
    for (std::size_t i = 0; i < m_scales.size(); ++ i)
      launch_feature_computation (new Feature_adder_echo<EchoMap> (this, echo_map, i));
  }

  void generate_echo_based_features(const CGAL::Default_property_map<Iterator, std::size_t>&)
  {
  }

  void generate_gradient_features()
  {
#ifdef CGAL_CLASSIFICATION_USE_GRADIENT_OF_FEATURE
    std::size_t size = m_features->size();

    for (std::size_t i = 0; i < size; ++ i)
    {
      for (int j = m_scales.size() - 1; j >= 0; -- j)
      {
        std::ostringstream oss;
        oss << "_" << j;
        if ((*m_features)[i]->name().find (oss.str()))
        {
          const Neighbor_query& neighbor_query = neighborhood(std::size_t(j)).k_neighbor_query(6);
          m_features->template add<Gradient_of_feature> (m_input, m_point_map, (*m_features)[i], neighbor_query);
          break;
        }
      }
    }
#endif
  }


  template <typename T>
  const T& get_parameter (const T& t)
  {
    return t;
  }

  template <typename T>
  Default_property_map<Iterator, T>
  get_parameter (const Default&)
  {
    return Default_property_map<Iterator, T>();
  }

  template<typename VectorMap, typename ColorMap, typename EchoMap>
  void generate_features_impl (std::size_t nb_scales, float voxel_size,
                               VectorMap normal_map,
                               ColorMap color_map,
                               EchoMap echo_map)
  {
    CGAL::Real_timer t; t.start();
    
    m_scales.reserve (nb_scales);
    
    m_scales.push_back (new Scale (m_input, m_point_map, m_bbox, -1.f));

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

    t.start();
    
#ifdef CGAL_LINKED_WITH_TBB
    m_tasks = new tbb::task_group;
#endif
    
    generate_point_based_features ();
    generate_normal_based_features (normal_map);
    generate_color_based_features (color_map);
    generate_echo_based_features (echo_map);

#ifdef CGAL_LINKED_WITH_TBB
    m_tasks->wait();
    delete m_tasks;
#endif

    generate_gradient_features();
    
    t.stop();
    CGAL_CLASSIFICATION_CERR << "Features computed in " << t.time() << " second(s)" << std::endl;
    for (std::size_t i = 0; i < m_adders.size(); ++ i)
      delete m_adders[i];
  }

  template <typename Feature_type>
  struct Feature_adder_variant_0 : public Feature_adder
  {
    using Feature_adder::generator;
    using Feature_adder::scale;

    Feature_adder_variant_0 (Point_set_feature_generator* generator, std::size_t scale)
      : Feature_adder (generator, scale) { }
    
    void operator() () const
    {
      Feature_handle fh = generator->m_features->template add<Feature_type> (generator->m_input, generator->eigen(scale));
      std::ostringstream oss;
      oss << fh->name() << "_" << scale;
      fh->set_name (oss.str());
    }
  };
  
  template <typename Feature_type>
  void generate_multiscale_feature_variant_0 ()
  {
    for (std::size_t i = 0; i < m_scales.size(); ++ i)
      launch_feature_computation (new Feature_adder_variant_0<Feature_type> (this, i));
  }

  template <typename Feature_type>
  struct Feature_adder_variant_1 : public Feature_adder
  {
    using Feature_adder::generator;
    using Feature_adder::scale;
    PointMap point_map;

    Feature_adder_variant_1 (Point_set_feature_generator* generator, PointMap point_map, std::size_t scale)
      : Feature_adder (generator, scale), point_map (point_map) { }
    
    void operator() () const
    {
      Feature_handle fh = generator->m_features->template add<Feature_type> (generator->m_input, point_map,
                                                                             generator->eigen(scale));
      std::ostringstream oss;
      oss << fh->name() << "_" << scale;
      fh->set_name (oss.str());
    }
  };
  
  template <typename Feature_type>
  void generate_multiscale_feature_variant_1 ()
  {
    for (std::size_t i = 0; i < m_scales.size(); ++ i)
      launch_feature_computation (new Feature_adder_variant_1<Feature_type> (this, m_point_map, i));
  }

  template <typename Feature_type>
  struct Feature_adder_variant_2 : public Feature_adder
  {
    using Feature_adder::generator;
    using Feature_adder::scale;
    PointMap point_map;

    Feature_adder_variant_2 (Point_set_feature_generator* generator, PointMap point_map, std::size_t scale)
      : Feature_adder (generator, scale), point_map (point_map) { }
    
    void operator() () const
    {
      Feature_handle fh = generator->m_features->template add<Feature_type>
        (generator->m_input, point_map,
         generator->grid(scale),
         generator->radius_neighbors(scale));
      std::ostringstream oss;
      oss << fh->name() << "_" << scale;
      fh->set_name (oss.str());
    }
  };
  
  template <typename Feature_type>
  void generate_multiscale_feature_variant_2 ()
  {
    for (std::size_t i = 0; i < m_scales.size(); ++ i)
      launch_feature_computation (new Feature_adder_variant_2<Feature_type> (this, m_point_map, i));
  }

  template <typename Feature_type>
  struct Feature_adder_variant_3 : public Feature_adder
  {
    using Feature_adder::generator;
    using Feature_adder::scale;
    PointMap point_map;

    Feature_adder_variant_3 (Point_set_feature_generator* generator, PointMap point_map, std::size_t scale)
      : Feature_adder (generator, scale), point_map (point_map) { }
    
    void operator() () const
    {
      Feature_handle fh = generator->m_features->template add<Feature_type> (generator->m_input,
                                                                             point_map,
                                                                             generator->grid(scale),
                                                                             generator->radius_dtm(scale));
      std::ostringstream oss;
      oss << fh->name() << "_" << scale;
      fh->set_name (oss.str());
    }
  };
  
  template <typename Feature_type>
  void generate_multiscale_feature_variant_3 ()
  {
    for (std::size_t i = 0; i < m_scales.size(); ++ i)
      launch_feature_computation (new Feature_adder_variant_3<Feature_type> (this, m_point_map, i));
  }

  
};


} // namespace Classification
  
} // namespace CGAL



#endif // CGAL_CLASSIFICATION_POINT_SET_FEATURE_GENERATOR_H
