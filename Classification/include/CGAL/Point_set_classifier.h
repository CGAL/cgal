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
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_POINT_SET_CLASSIFIER_H
#define CGAL_POINT_SET_CLASSIFIER_H

#include <CGAL/property_map.h>
#include <CGAL/Classifier.h>
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
#include <CGAL/Classification/Label.h>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include <CGAL/Real_timer.h>
#include <CGAL/demangle.h>

#ifdef CGAL_LINKED_WITH_TBB
#include <tbb/parallel_invoke.h>
#include <tbb/blocked_range.h>
#include <tbb/scalable_allocator.h>
#endif // CGAL_LINKED_WITH_TBB


namespace CGAL {

/*!
  \ingroup PkgClassification

  \brief Classifies a point set based on a set of features and a set
  of labels.

  This class specializes `Classifier` to point sets. It takes care of
  generating necessary data structures and automatically generate a
  set of generic features. Features can be generated at multiple
  scales to increase the reliability of the classification.

  \tparam Geom_traits model of \cgal Kernel.
  \tparam PointRange model of `ConstRange`. Its iterator type is
  `RandomAccessIterator`.
  \tparam PointMap model of `ReadablePropertyMap` whose key
  type is the value type of the iterator of `PointRange` and value type
  is `Geom_traits::Point_3`.
  \tparam DiagonalizeTraits model of `DiagonalizeTraits` used
  for matrix diagonalization.

*/
template <typename Geom_traits,
          typename PointRange,
          typename PointMap,
          typename DiagonalizeTraits = CGAL::Default_diagonalize_traits<double,3> >
class Point_set_classifier : public Classifier<PointRange, PointMap>
{
  
public:
  typedef typename Geom_traits::Iso_cuboid_3             Iso_cuboid_3;

  /// \cond SKIP_IN_MANUAL
  typedef Classifier<PointRange, PointMap> Base;
  typedef typename PointRange::const_iterator Iterator;
  using Base::m_input;
  using Base::m_item_map;
  /// \endcond
  
  typedef Classification::Planimetric_grid
  <Geom_traits, PointRange, PointMap>                    Planimetric_grid;
  typedef Classification::Point_set_neighborhood
  <Geom_traits, PointRange, PointMap>                    Neighborhood;
  typedef Classification::Local_eigen_analysis
  <Geom_traits, PointRange, PointMap, DiagonalizeTraits> Local_eigen_analysis;

  /// \cond SKIP_IN_MANUAL
  typedef Classification::Feature_handle                 Feature_handle;
  typedef Classification::Label                          Label;
  typedef Classification::Label_handle                   Label_handle;

  typedef typename Geom_traits::Point_3                  Point;
  
  typedef Classification::Feature::Anisotropy
  <Geom_traits, PointRange, PointMap, DiagonalizeTraits> Anisotropy;
  typedef Classification::Feature::Distance_to_plane
  <Geom_traits, PointRange, PointMap, DiagonalizeTraits> Distance_to_plane;
  typedef Classification::Feature::Eigentropy
  <Geom_traits, PointRange, PointMap, DiagonalizeTraits> Eigentropy;
  typedef Classification::Feature::Elevation
  <Geom_traits, PointRange, PointMap>                    Elevation;
  typedef Classification::Feature::Linearity
  <Geom_traits, PointRange, PointMap, DiagonalizeTraits> Linearity;
  typedef Classification::Feature::Omnivariance
  <Geom_traits, PointRange, PointMap, DiagonalizeTraits> Omnivariance;
  typedef Classification::Feature::Planarity
  <Geom_traits, PointRange, PointMap, DiagonalizeTraits> Planarity;
  typedef Classification::Feature::Sphericity
  <Geom_traits, PointRange, PointMap, DiagonalizeTraits> Sphericity;
  typedef Classification::Feature::Sum_eigenvalues
  <Geom_traits, PointRange, PointMap, DiagonalizeTraits> Sum_eigen;
  typedef Classification::Feature::Surface_variation
  <Geom_traits, PointRange, PointMap, DiagonalizeTraits> Surface_variation;
  typedef Classification::Feature::Vertical_dispersion
  <Geom_traits, PointRange, PointMap>                    Dispersion;
  typedef Classification::Feature::Verticality
  <Geom_traits, PointRange, PointMap, DiagonalizeTraits> Verticality;
  typedef typename Classification::RGB_Color RGB_Color;
  /// \endcond
    
private:

  struct Scale
  {
    Neighborhood* neighborhood;
    Planimetric_grid* grid;
    Local_eigen_analysis* eigen;
    double voxel_size;
    std::vector<Feature_handle> features;
    
    Scale (const PointRange& input, PointMap point_map,
           const Iso_cuboid_3& bbox, double voxel_size)
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
      
      eigen = new Local_eigen_analysis (input, point_map, neighborhood->k_neighbor_query(6));
      double range = eigen->mean_range();
      if (this->voxel_size < 0)
        this->voxel_size = range;
      t.stop();
      CGAL_CLASSIFICATION_CERR << "Eigen values computed in " << t.time() << " second(s)" << std::endl;
      CGAL_CLASSIFICATION_CERR << "Range = " << range << std::endl;
      t.reset();
      t.start();
      
      grid = new Planimetric_grid (input, point_map, bbox, this->voxel_size);
      t.stop();
      CGAL_CLASSIFICATION_CERR << "Planimetric grid computed in " << t.time() << " second(s)" << std::endl;
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

  
  /// \name Constructor
  /// @{
  
  /*! 
    \brief Initializes a classification object.

    \param input input range.

    \param point_map property map to access the input points.
  */
  Point_set_classifier(const PointRange& input, PointMap point_map) : Base (input, point_map)
  {
    m_bbox = CGAL::bounding_box
      (boost::make_transform_iterator (m_input.begin(), CGAL::Property_map_to_unary_function<PointMap>(m_item_map)),
       boost::make_transform_iterator (m_input.end(), CGAL::Property_map_to_unary_function<PointMap>(m_item_map)));
  }

  /// @}
  
  /// \cond SKIP_IN_MANUAL
  virtual ~Point_set_classifier()
  {
    clear();
  }
  /// \endcond

  /// \name Features
  /// @{

  
  /*!
    \brief Generate all possible features from an input range.

    The size of the smallest scale is automatically estimated and the
    data structures needed (`Neighborhood`, `Planimetric_grid` and
    `Local_eigen_analysis`) are computed at `nb_scales` recursively
    larger scales. At each scale, the following features are
    generated:

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
    - `CGAL::Classification::Feature::Vertical_dispersion` based on eigenvalues

    If normal vectors are provided (if `VectorMap` is different from
    `CGAL::Default`), the following feature is generated at each
    scale:

    - `CGAL::Classification::Feature::Vertical_dispersion` based on normal vectors

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
    `Geom_traits::Vector_3`.
    \tparam ColorMap model of `ReadablePropertyMap`  whose key type is
    the value type of the iterator of `PointRange` and value type is
    `CGAL::Classification::RGB_Color`.
    \tparam EchoMap model of `ReadablePropertyMap` whose key type is
    the value type of the iterator of `PointRange` and value type is
    `std::size_t`.
    \param nb_scales number of scales to compute.
    \param normal_map property map to access the normal vectors of the input points (if any).
    \param color_map property map to access the colors of the input points (if any).
    \param echo_map property map to access the echo values of the input points (if any).
  */
  template <typename VectorMap = Default,
            typename ColorMap = Default,
            typename EchoMap = Default>
  void generate_features (std::size_t nb_scales,
                            VectorMap normal_map = VectorMap(),
                            ColorMap color_map = ColorMap(),
                            EchoMap echo_map = EchoMap())
  {
    typedef typename Default::Get<VectorMap, typename Geom_traits::Vector_3 >::type
      Vmap;
    typedef typename Default::Get<ColorMap, RGB_Color >::type
      Cmap;
    typedef typename Default::Get<EchoMap, std::size_t >::type
      Emap;

    generate_features_impl (nb_scales,
                              get_parameter<Vmap>(normal_map),
                              get_parameter<Cmap>(color_map),
                              get_parameter<Emap>(echo_map));
  }

  /// @}

  /// \name Data Structures and Parameters
  /// @{
  
  /*!
    \brief Returns the bounding box of the input point set.
  */
  const Iso_cuboid_3& bbox() const { return m_bbox; }
  /*!
    \brief Returns the neighborhood structure at scale `scale`.

    \note `generate_features()` must have been called before calling
    this method.
  */
  const Neighborhood& neighborhood(std::size_t scale = 0) const { return (*m_scales[scale]->neighborhood); }
  /*!
    \brief Returns the planimetric grid structure at scale `scale`.

    \note `generate_features()` must have been called before calling
    this method.
  */
  const Planimetric_grid& grid(std::size_t scale = 0) const { return *(m_scales[scale]->grid); }
  /*!
    \brief Returns the local eigen analysis structure at scale `scale`.

    \note `generate_features()` must have been called before calling
    this method.
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

    \note `generate_features()` must have been called before calling
    this method.
  */
  double grid_resolution(std::size_t scale = 0) const { return m_scales[scale]->grid_resolution(); }
  /*!

    \brief Returns the radius used for neighborhood queries at scale
    `scale`. This radius is the smallest radius that is relevant from
    a geometric point of view at this scale (that is to say that
    encloses a few cells of `Planimetric_grid`).

    \note `generate_features()` must have been called before calling
    this method.
  */
  double radius_neighbors(std::size_t scale = 0) const { return m_scales[scale]->radius_neighbors(); }
  /*!
    \brief Returns the radius used for digital terrain modeling at
    scale `scale`. This radius represents the minimum size of a
    building at this scale.

    \note `generate_features()` must have been called before calling
    this method.
  */
  double radius_dtm(std::size_t scale = 0) const { return m_scales[scale]->radius_dtm(); }

  /*!
    \brief Clears all computed data structures.
  */
  void clear()
  {
    for (std::size_t i = 0; i < m_scales.size(); ++ i)
      delete m_scales[i];
    m_scales.clear();
    
    this->clear_labels();
    this->clear_features();
  }

  /// @}
  
  /// @}

  /// \cond SKIP_IN_MANUAL  
  void info() const
  {
    CGAL_CLASSIFICATION_CERR << m_scales.size() << " scale(s) used:" << std::endl;
    for (std::size_t i = 0; i < m_scales.size(); ++ i)
      {
        std::size_t nb_useful = 0;
        for (std::size_t j = 0; j < m_scales[i]->features.size(); ++ j)
          if (m_scales[i]->features[j]->weight() != 0.)
            nb_useful ++;
        CGAL_CLASSIFICATION_CERR << " * scale " << i << " with size " << m_scales[i]->voxel_size
                  << ", " << nb_useful << " useful feature(s)";
        if (nb_useful != 0) CGAL_CLASSIFICATION_CERR << ":" << std::endl;
        else CGAL_CLASSIFICATION_CERR << std::endl;
        for (std::size_t j = 0; j < m_scales[i]->features.size(); ++ j)
          if (m_scales[i]->features[j]->weight() != 0.)
            CGAL_CLASSIFICATION_CERR << "   - " << m_scales[i]->features[j]->name()
                      << " (weight = " << m_scales[i]->features[j]->weight() << ")" << std::endl;
      }
  }

  void generate_point_based_features ()
  {
    CGAL::Real_timer teigen, tpoint;
    teigen.start();
    generate_multiscale_feature_variant_0<Anisotropy> ();
    generate_multiscale_feature_variant_0<Eigentropy> ();
    generate_multiscale_feature_variant_0<Linearity> ();
    generate_multiscale_feature_variant_0<Omnivariance> ();
    generate_multiscale_feature_variant_0<Planarity> ();
    generate_multiscale_feature_variant_0<Sphericity> ();
    generate_multiscale_feature_variant_0<Sum_eigen> ();
    generate_multiscale_feature_variant_0<Surface_variation> ();

    teigen.stop();
    tpoint.start();
    generate_multiscale_feature_variant_1<Distance_to_plane> ();
    generate_multiscale_feature_variant_2<Dispersion> ();
    generate_multiscale_feature_variant_3<Elevation> ();
    tpoint.stop();
    
    CGAL_CLASSIFICATION_CERR << "Point based features computed in " << tpoint.time() << " second(s)" << std::endl;
    CGAL_CLASSIFICATION_CERR << "Eigen based features computed in " << teigen.time() << " second(s)" << std::endl;
  }


  template <typename VectorMap>
  void generate_normal_based_features(VectorMap normal_map)
  {
    CGAL::Real_timer t; t.start();
    this->template add_feature<Verticality> (normal_map);
    m_scales[0]->features.push_back (this->feature (this->number_of_features() - 1));
    t.stop();
    CGAL_CLASSIFICATION_CERR << "Normal based features computed in " << t.time() << " second(s)" << std::endl;
  }

  void generate_normal_based_features(const CGAL::Default_property_map<Iterator, typename Geom_traits::Vector_3>&)
  {
    CGAL::Real_timer t; t.start();
    generate_multiscale_feature_variant_0<Verticality> ();
    CGAL_CLASSIFICATION_CERR << "Normal based features computed in " << t.time() << " second(s)" << std::endl;
  }
  template <typename ColorMap>
  void generate_color_based_features(ColorMap color_map)
  {

    typedef Classification::Feature::Hsv<Geom_traits, PointRange, ColorMap> Hsv;
    CGAL::Real_timer t; t.start();
    for (std::size_t i = 0; i <= 8; ++ i)
      {
        this->template add_feature<Hsv> (color_map, 0, 45 * i, 22.5);
        m_scales[0]->features.push_back (this->feature (this->number_of_features() - 1));
      }
    for (std::size_t i = 0; i <= 4; ++ i)
      {
        this->template add_feature<Hsv> (color_map, 1, 25 * i, 12.5);
        m_scales[0]->features.push_back (this->feature (this->number_of_features() - 1));
      }
    for (std::size_t i = 0; i <= 4; ++ i)
      {
        this->template add_feature<Hsv> (color_map, 2, 25 * i, 12.5);
        m_scales[0]->features.push_back (this->feature (this->number_of_features() - 1));
      }
    t.stop();
    CGAL_CLASSIFICATION_CERR << "Color based features computed in " << t.time() << " second(s)" << std::endl;
  }

  void generate_color_based_features(const CGAL::Default_property_map<Iterator, RGB_Color>&)
  {
  }

  template <typename EchoMap>
  void generate_echo_based_features(EchoMap echo_map)
  {
    typedef Classification::Feature::Echo_scatter<Geom_traits, PointRange, PointMap, EchoMap> Echo_scatter;
    CGAL::Real_timer t; t.start();
    for (std::size_t i = 0; i < m_scales.size(); ++ i)
      {
        this->template add_feature<Echo_scatter>(echo_map, *(m_scales[i]->grid),
                                          m_scales[i]->grid_resolution(),
                                          m_scales[i]->radius_neighbors());
        m_scales[i]->features.push_back (this->feature (this->number_of_features() - 1));
      }
    t.stop();
    CGAL_CLASSIFICATION_CERR << "Echo based features computed in " << t.time() << " second(s)" << std::endl;
  }

  void generate_echo_based_features(const CGAL::Default_property_map<Iterator, std::size_t>&)
  {
  }

  void get_map_scale (std::map<Feature_handle, std::size_t>& map_scale)
  {
    for (std::size_t i = 0; i < m_scales.size(); ++ i)
      for (std::size_t j = 0; j < m_scales[i]->features.size(); ++ j)
        map_scale[m_scales[i]->features[j]] = i;
  }

  std::size_t scale_of_feature (Feature_handle att)
  {
    for (std::size_t i = 0; i < m_scales.size(); ++ i)
      for (std::size_t j = 0; j < m_scales[i]->features.size(); ++ j)
        if (m_scales[i]->features[j] == att)
          return i;
    return (std::size_t)(-1);    
  }

  std::string name_att (Feature_handle att,
                        std::map<Feature_handle, std::size_t>& map_scale)
  {
    std::ostringstream oss;
    oss << att->name() << "_" << map_scale[att];
    return oss.str();
  }
  /// \endcond

  /// \name Input/Output
  /// @{

  /*!
    \brief Saves the current configuration in the stream `output`.

    This allows to easily save and recover a specific classification
    configuration, that is to say:

    - The size of the smallest scale
    - The features and their respective weights
    - The labels and the effects of the features on them

    The output file is written in an XML format that is readable by
    the `load_configuration()` method.
  */
  void save_configuration (std::ostream& output)
  {
    boost::property_tree::ptree tree;

    // tree.put("classification.parameters.grid_resolution", m_grid_resolution);
    // tree.put("classification.parameters.radius_neighbors", m_radius_neighbors);
    tree.put("classification.parameters.voxel_size", m_scales[0]->voxel_size);

    std::map<Feature_handle, std::size_t> map_scale;
    get_map_scale (map_scale);

    for (std::size_t i = 0; i < this->number_of_features(); ++ i)
      {
        Feature_handle att = this->feature(i);
        if (att->weight() == 0)
          continue;
        boost::property_tree::ptree ptr;
        
        ptr.put("id", name_att (att, map_scale));
        ptr.put("weight", att->weight());
        tree.add_child("classification.features.feature", ptr);
      }


    for (std::size_t i = 0; i < this->number_of_labels(); ++ i)
      {
        Label_handle label = this->label(i);
        boost::property_tree::ptree ptr;
        ptr.put("id", label->name());
        for (std::size_t j = 0; j < this->number_of_features(); ++ j)
          {
            Feature_handle att = this->feature(j);
            if (att->weight() == 0)
              continue;
            boost::property_tree::ptree ptr2;
            ptr2.put("id", name_att (att, map_scale));
            Classification::Feature::Effect effect = label->feature_effect(att);
            if (effect == Classification::Feature::PENALIZING)
              ptr2.put("effect", "penalized");
            else if (effect == Classification::Feature::NEUTRAL)
              ptr2.put("effect", "neutral");
            else if (effect == Classification::Feature::FAVORING)
              ptr2.put("effect", "favored");
            ptr.add_child("feature", ptr2);
          }
        tree.add_child("classification.labels.label", ptr);
      }

    // Write property tree to XML file
    boost::property_tree::xml_writer_settings<std::string> settings(' ', 3);
    boost::property_tree::write_xml(output, tree, settings);
  }

  
  /*!
    \brief Loads a configuration from the stream `input`.

    All data structures, features and labels specified in the input
    stream `input` are instantiated if possible (in particular,
    property maps needed should be provided), similarly to what is
    done in `generate_features()`.

    The input file should be in the XML format written by the
    `save_configuration()` method.

    \tparam VectorMap model of `ReadablePropertyMap` whose key type is
    the value type of the iterator of `PointRange` and value type is
    `Geom_traits::Vector_3`.
    \tparam ColorMap model of `ReadablePropertyMap`  whose key type is
    the value type of the iterator of `PointRange` and value type is
    `CGAL::Classification::RGB_Color`.
    \tparam EchoMap model of `ReadablePropertyMap` whose key type is
    the value type of the iterator of `PointRange` and value type is
    `std::size_t`.
    \param input input stream.
    \param normal_map property map to access the normal vectors of the input points (if any).
    \param color_map property map to access the colors of the input points (if any).
    \param echo_map property map to access the echo values of the input points (if any).
  */
  template<typename VectorMap = Default,
           typename ColorMap = Default,
           typename EchoMap  = Default>
  bool load_configuration (std::istream& input, 
                           VectorMap normal_map = VectorMap(),
                           ColorMap color_map = ColorMap(),
                           EchoMap echo_map = EchoMap())
  {
    typedef typename Default::Get<VectorMap, typename Geom_traits::Vector_3 >::type
      Vmap;
    typedef typename Default::Get<ColorMap, RGB_Color >::type
      Cmap;
    typedef typename Default::Get<EchoMap, std::size_t >::type
      Emap;

    return load_configuration_impl (input,
                                    get_parameter<Vmap>(normal_map),
                                    get_parameter<Cmap>(color_map),
                                    get_parameter<Emap>(echo_map));
  }



  /*!

    \brief Writes the classified point set in a colored and labeled
    PLY format in the stream `output`.

    The input points are written in a PLY format with the addition of
    the following PLY properties:

    - a property `label` that indicates which label is
    assigned to the point. The labels are indexed from 0 to N (the
    correspondancy is given as comments in the PLY header).

    - 3 properties `red`, `green` and `blue` to associate each label
    to a color (this is useful to visualize the classification in a
    viewer that supports PLY colors). Colors are picked randomly.
  */
  void write_classification_to_ply (std::ostream& output)
  {
    output << "ply" << std::endl
           << "format ascii 1.0" << std::endl
           << "comment Generated by the CGAL library www.cgal.org" << std::endl
           << "element vertex " << m_input.size() << std::endl
           << "property double x" << std::endl
           << "property double y" << std::endl
           << "property double z" << std::endl
           << "property uchar red" << std::endl
           << "property uchar green" << std::endl
           << "property uchar blue" << std::endl
           << "property int label" << std::endl;

    std::vector<RGB_Color> colors;
    
    std::map<Label_handle, std::size_t> map_labels;
    output << "comment label -1 is (unclassified)" << std::endl;

    for (std::size_t i = 0; i < this->number_of_labels(); ++ i)
      {
        map_labels.insert (std::make_pair (this->label(i), i));
        output << "comment label " << i << " is " << this->label(i)->name() << std::endl;
        
        RGB_Color c = {{ (unsigned char)(64 + rand() % 128),
                         (unsigned char)(64 + rand() % 128),
                         (unsigned char)(64 + rand() % 128) }};

        if (this->label(i)->name() == "vegetation")
          c = {{ 0, 255, 27 }};
        else if (this->label(i)->name() == "ground")
          c = {{ 245, 180, 0 }};
        else if (this->label(i)->name() == "roof")
          c = {{ 255, 0, 170 }};
        else if (this->label(i)->name() == "facade")
          c = {{ 100, 0, 255 }};
        colors.push_back (c);
      }
    map_labels.insert (std::make_pair (Label_handle(), this->number_of_labels()));
    
    output << "end_header" << std::endl;

    std::size_t i = 0;
    for (Iterator it = m_input.begin(); it != m_input.end(); ++ it)
      {
        Label_handle t = this->label_of(i);
        std::size_t idx = map_labels[t];

        if (idx == this->number_of_labels())
          output << get(m_item_map, *it) << " 0 0 0 -1" << std::endl;
        else
          output << get(m_item_map, *it) << " "
                 << (int)(colors[idx][0]) << " "
                 << (int)(colors[idx][1]) << " "
                 << (int)(colors[idx][2]) << " "
                 << idx << std::endl;
        ++ i;
      }
  }

  /// @}
  
private:
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
  void generate_features_impl (std::size_t nb_scales,
                                 VectorMap normal_map,
                                 ColorMap color_map,
                                 EchoMap echo_map)
  {
    CGAL::Real_timer t; t.start();

    m_scales.reserve (nb_scales);
    double voxel_size = - 1.;

    m_scales.push_back (new Scale (m_input, m_item_map, m_bbox, voxel_size));
    voxel_size = m_scales[0]->grid_resolution();
    
    for (std::size_t i = 1; i < nb_scales; ++ i)
      {
        voxel_size *= 2;
        m_scales.push_back (new Scale (m_input, m_item_map, m_bbox, voxel_size));
      }
    
    generate_point_based_features ();
    generate_normal_based_features (normal_map);
    generate_color_based_features (color_map);
    generate_echo_based_features (echo_map);
  }

  template <typename Feature_type>
  void generate_multiscale_feature_variant_0 ()
  {
    for (std::size_t i = 0; i < m_scales.size(); ++ i)
      {
        this->template add_feature<Feature_type>(*(m_scales[i]->eigen));
        m_scales[i]->features.push_back (this->feature (this->number_of_features() - 1));
      }
  }

  template <typename Feature_type>
  void generate_multiscale_feature_variant_1 ()
  {
    for (std::size_t i = 0; i < m_scales.size(); ++ i)
      {
        this->template add_feature<Feature_type>(m_item_map, *(m_scales[i]->eigen));
        m_scales[i]->features.push_back (this->feature (this->number_of_features() - 1));
      }
  }

  template <typename Feature_type>
  void generate_multiscale_feature_variant_2 ()
  {
    for (std::size_t i = 0; i < m_scales.size(); ++ i)
      {
        this->template add_feature<Feature_type>(m_item_map,
                                            *(m_scales[i]->grid),
                                            m_scales[i]->grid_resolution(),
                                            m_scales[i]->radius_neighbors());
        m_scales[i]->features.push_back (this->feature (this->number_of_features() - 1));
      }
  }

  template <typename Feature_type>
  void generate_multiscale_feature_variant_3 ()
  {
    for (std::size_t i = 0; i < m_scales.size(); ++ i)
      {
        this->template add_feature<Feature_type>(m_item_map,
                                            *(m_scales[i]->grid),
                                            m_scales[i]->grid_resolution(),
                                            m_scales[i]->radius_dtm());
        m_scales[i]->features.push_back (this->feature (this->number_of_features() - 1));
      }
  }

  template<typename VectorMap,typename ColorMap, typename EchoMap>
  bool load_configuration_impl (std::istream& input, 
                                VectorMap normal_map,
                                ColorMap color_map,
                                EchoMap echo_map)
  {
    typedef Classification::Feature::Echo_scatter<Geom_traits, PointRange, PointMap, EchoMap> Echo_scatter;
    typedef Classification::Feature::Hsv<Geom_traits, PointRange, ColorMap> Hsv;
    
    clear();
    
    boost::property_tree::ptree tree;
    boost::property_tree::read_xml(input, tree);

    double voxel_size = tree.get<double>("classification.parameters.voxel_size");
    
    m_scales.push_back (new Scale (m_input, m_item_map, m_bbox, voxel_size));
    
    CGAL::Real_timer t;
    std::map<std::string, Feature_handle> att_map;
    BOOST_FOREACH(boost::property_tree::ptree::value_type &v, tree.get_child("classification.features"))
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
            m_scales.push_back (new Scale (m_input, m_item_map, m_bbox, voxel_size));
          }
        
        double weight = v.second.get<double>("weight");

        // Generate the right feature if possible
        if (id == "anisotropy")
          this->template add_feature<Anisotropy>(*(m_scales[scale]->eigen));
        else if (id == "distance_to_plane")
          this->template add_feature<Distance_to_plane>(m_item_map, *(m_scales[scale]->eigen));
        else if (id == "eigentropy")
          this->template add_feature<Eigentropy>(*(m_scales[scale]->eigen));
        else if (id == "elevation")
          {
            t.start();
            this->template add_feature<Elevation>(m_item_map,
                                           *(m_scales[scale]->grid),
                                           m_scales[scale]->grid_resolution(),
                                           m_scales[scale]->radius_dtm());
            t.stop();
          }
        else if (id == "linearity")
          this->template add_feature<Linearity>(*(m_scales[scale]->eigen));
        else if (id == "omnivariance")
          this->template add_feature<Omnivariance>(*(m_scales[scale]->eigen));
        else if (id == "planarity")
          this->template add_feature<Planarity>(*(m_scales[scale]->eigen));
        else if (id == "sphericity")
          this->template add_feature<Sphericity>(*(m_scales[scale]->eigen));
        else if (id == "sum_eigen")
          this->template add_feature<Sum_eigen>(*(m_scales[scale]->eigen));
        else if (id == "surface_variation")
          this->template add_feature<Surface_variation>(*(m_scales[scale]->eigen));
        else if (id == "vertical_dispersion")
          this->template add_feature<Dispersion>(m_item_map,
                                          *(m_scales[scale]->grid),
                                          m_scales[scale]->grid_resolution(),
                                          m_scales[scale]->radius_neighbors());
        else if (id == "verticality")
          {
            if (boost::is_convertible<VectorMap,
                typename CGAL::Default_property_map<Iterator, typename Geom_traits::Vector_3> >::value)
              this->template add_feature<Verticality>(*(m_scales[scale]->eigen));
            else
              this->template add_feature<Verticality>(normal_map);
          }
        else if (id == "echo_scatter")
          {
            if (boost::is_convertible<EchoMap,
                typename CGAL::Default_property_map<Iterator, std::size_t> >::value)
              {
                CGAL_CLASSIFICATION_CERR << "Warning: echo_scatter required but no echo map given." << std::endl;
                continue;
              }
            this->template add_feature<Echo_scatter>(echo_map, *(m_scales[scale]->grid),
                                              m_scales[scale]->grid_resolution(),
                                              m_scales[scale]->radius_neighbors());
          }
        else if (boost::starts_with(id.c_str(), "hue")
                 || boost::starts_with(id.c_str(), "saturation")
                 || boost::starts_with(id.c_str(), "value"))
          {
            if (boost::is_convertible<ColorMap,
                typename CGAL::Default_property_map<Iterator, RGB_Color> >::value)
              {
                CGAL_CLASSIFICATION_CERR << "Warning: color feature required but no color map given." << std::endl;
                continue;
              }
            if (boost::starts_with(id.c_str(), "hue"))
              {
                double value = boost::lexical_cast<int>(id.c_str() + 4);
                this->template add_feature<Hsv>(color_map, 0, value, 22.5);
              }
            else if (boost::starts_with(id.c_str(), "saturation"))
              {
                double value = boost::lexical_cast<int>(id.c_str() + 11);
                this->template add_feature<Hsv>(color_map, 1, value, 12.5);
              }
            else if (boost::starts_with(id.c_str(), "value"))
              {
                double value = boost::lexical_cast<int>(id.c_str() + 6);
                this->template add_feature<Hsv>(color_map, 2, value, 12.5);
              }
          }
        else
          {
            CGAL_CLASSIFICATION_CERR << "Warning: unknown feature \"" << id << "\"" << std::endl;
            continue;
          }

        Feature_handle att = this->feature (this->number_of_features() - 1);
        m_scales[scale]->features.push_back (att);
        att->set_weight(weight);
        att_map[full_id] = att;
      }

    BOOST_FOREACH(boost::property_tree::ptree::value_type &v, tree.get_child("classification.labels"))
      {
        std::string label_id = v.second.get<std::string>("id");

        Label_handle new_label = this->add_label (label_id.c_str());
        
        BOOST_FOREACH(boost::property_tree::ptree::value_type &v2, v.second)
          {
            if (v2.first == "id")
              continue;
            std::string att_id = v2.second.get<std::string>("id");
            std::map<std::string, Feature_handle>::iterator it = att_map.find(att_id);
            if (it == att_map.end())
              continue;
            Feature_handle att = it->second;
            std::string effect = v2.second.get<std::string>("effect");
            if (effect == "penalized")
              new_label->set_feature_effect (att, Classification::Feature::PENALIZING);
            else if (effect == "neutral")
              new_label->set_feature_effect (att, Classification::Feature::NEUTRAL);
            else
              new_label->set_feature_effect (att, Classification::Feature::FAVORING);
          }
      }
    
    return true;
  }
};

} // namespace CGAL


#endif // CGAL_POINT_SET_CLASSIFIER_H
