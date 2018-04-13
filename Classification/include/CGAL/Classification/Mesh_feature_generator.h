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

#ifndef CGAL_CLASSIFICATION_MESH_FEATURE_GENERATOR_H
#define CGAL_CLASSIFICATION_MESH_FEATURE_GENERATOR_H

#include <CGAL/license/Classification.h>

#include <CGAL/property_map.h>
#include <CGAL/Classification/Mesh_neighborhood.h>
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

// New features
#include <CGAL/Classification/Feature/Eigenvalue.h>
#include <CGAL/Classification/Feature/Color_channel.h>

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

*/
template <typename Geom_traits,
          typename FaceListGraph,
          typename PointMap,
#if defined(DOXYGEN_RUNNING)
          typename ConcurrencyTag,
#elif defined(CGAL_LINKED_WITH_TBB)
          typename ConcurrencyTag = CGAL::Parallel_tag,
#else
          typename ConcurrencyTag = CGAL::Sequential_tag,
#endif
          typename DiagonalizeTraits = CGAL::Default_diagonalize_traits<double,3> >
class Mesh_feature_generator
{
  
public:
  typedef typename Geom_traits::Iso_cuboid_3             Iso_cuboid_3;

  /// \cond SKIP_IN_MANUAL
  typedef typename boost::graph_traits<FaceListGraph>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<FaceListGraph>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<FaceListGraph>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<FaceListGraph>::face_iterator face_iterator;
  typedef typename CGAL::Iterator_range<face_iterator> Face_range;
  
  typedef typename PointMap::value_type       Point;
  typedef CGAL::Identity_property_map<face_descriptor> Face_map;
  /// \endcond


public:
  
  typedef Classification::Planimetric_grid
  <Geom_traits, Face_range, PointMap>                    Planimetric_grid;
  typedef Classification::Mesh_neighborhood
  <FaceListGraph>                                        Neighborhood;
  typedef Classification::Local_eigen_analysis           Local_eigen_analysis;

  /// \cond SKIP_IN_MANUAL
  typedef Classification::Feature_handle                 Feature_handle;
  
  typedef Classification::Feature::Anisotropy            Anisotropy;
  typedef Classification::Feature::Distance_to_plane
  <Face_range, PointMap>                                 Distance_to_plane;
  typedef Classification::Feature::Eigentropy            Eigentropy;
  typedef Classification::Feature::Elevation
  <Geom_traits, Face_range, PointMap>                    Elevation;
  typedef Classification::Feature::Linearity             Linearity;
  typedef Classification::Feature::Omnivariance          Omnivariance;
  typedef Classification::Feature::Planarity             Planarity;
  typedef Classification::Feature::Sphericity            Sphericity;
  typedef Classification::Feature::Sum_eigenvalues       Sum_eigen;
  typedef Classification::Feature::Surface_variation     Surface_variation;
  typedef Classification::Feature::Vertical_dispersion
  <Geom_traits, Face_range, PointMap>                    Dispersion;
  typedef Classification::Feature::Verticality
  <Geom_traits>                                          Verticality;
  typedef Classification::Feature::Eigenvalue           Eigenvalue;

  typedef typename Classification::RGB_Color RGB_Color;
  /// \endcond
    
private:

  struct Scale
  {
    Neighborhood* neighborhood;
    Planimetric_grid* grid;
    Local_eigen_analysis* eigen;
    float voxel_size;
    
    Scale (const FaceListGraph& input,
           const Face_range& range,
           PointMap point_map,
           const Iso_cuboid_3& bbox, float voxel_size,
           std::size_t nb_scale,
           Planimetric_grid* lower_grid = NULL)
      : voxel_size (voxel_size)
    {
      CGAL::Real_timer t;
      t.start();
      neighborhood = new Neighborhood (input);
      t.stop();
      
      CGAL_CLASSIFICATION_CERR << "Neighborhood computed in " << t.time() << " second(s)" << std::endl;

      t.reset();
      t.start();
      
      eigen = new Local_eigen_analysis (typename Local_eigen_analysis::Input_is_face_graph(),
                                        input, neighborhood->n_ring_neighbor_query(nb_scale + 1),
                                        ConcurrencyTag(), DiagonalizeTraits());
      float mrange = eigen->mean_range();
      if (this->voxel_size < 0)
        this->voxel_size = mrange;
      t.stop();
      CGAL_CLASSIFICATION_CERR << "Eigen values computed in " << t.time() << " second(s)" << std::endl;
      CGAL_CLASSIFICATION_CERR << "Range = " << mrange << std::endl;
      t.reset();
      t.start();

      if (lower_grid == NULL)
        grid = new Planimetric_grid (range, point_map, bbox, this->voxel_size);
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

  const FaceListGraph& m_input;
  Face_range m_range;
  PointMap m_point_map;
  Feature_set& m_features;
  std::vector<std::pair<Feature_handle, std::size_t> > m_features_to_rename;
  
public:

  
  /*!
    \brief Generates all possible features from an input range.
  */
  Mesh_feature_generator(Feature_set& features,
                         const FaceListGraph& input,
                         PointMap point_map,
                         std::size_t nb_scales)
    : m_input (input), m_range(faces(input)), m_point_map (point_map), m_features (features)
  {

    m_bbox = CGAL::bounding_box
      (boost::make_transform_iterator (m_range.begin(), CGAL::Property_map_to_unary_function<PointMap>(m_point_map)),
       boost::make_transform_iterator (m_range.end(), CGAL::Property_map_to_unary_function<PointMap>(m_point_map)));
    generate_features_impl (nb_scales);

  }

  
  /// \cond SKIP_IN_MANUAL
  virtual ~Mesh_feature_generator()
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

  void add (Feature_handle fh, std::size_t i)
  {
    m_features_to_rename.push_back (std::make_pair (fh, i));
  }

  void generate_point_based_features ()
  {
#ifdef DO_NOT_USE_EIGEN_FEATURES
    for (int j = 0; j < 3; ++ j)
    {
      for (std::size_t i = 0; i < m_scales.size(); ++ i)
        add(m_features.add<Eigenvalue> (m_range, eigen(i), std::size_t(j)), i);
    }
#else
    for (std::size_t i = 0; i < m_scales.size(); ++ i)
      add(m_features.add<Anisotropy> (m_range, eigen(i)), i);
    for (std::size_t i = 0; i < m_scales.size(); ++ i)
      add(m_features.add<Eigentropy> (m_range, eigen(i)), i);
    for (std::size_t i = 0; i < m_scales.size(); ++ i)
      add(m_features.add<Linearity> (m_range, eigen(i)), i);
    for (std::size_t i = 0; i < m_scales.size(); ++ i)
      add(m_features.add<Omnivariance> (m_range, eigen(i)), i);
    for (std::size_t i = 0; i < m_scales.size(); ++ i)
      add(m_features.add<Planarity> (m_range, eigen(i)), i);
    for (std::size_t i = 0; i < m_scales.size(); ++ i)
      add(m_features.add<Sphericity> (m_range, eigen(i)), i);
    for (std::size_t i = 0; i < m_scales.size(); ++ i)
      add(m_features.add<Sum_eigen> (m_range, eigen(i)), i);
    for (std::size_t i = 0; i < m_scales.size(); ++ i)
      add(m_features.add<Surface_variation> (m_range, eigen(i)), i);
#endif
    for (std::size_t i = 0; i < m_scales.size(); ++ i)
      add(m_features.add<Distance_to_plane> (m_range, m_point_map, eigen(i)), i);
    for (std::size_t i = 0; i < m_scales.size(); ++ i)
      add(m_features.add<Dispersion> (m_range, m_point_map, grid(i), radius_neighbors(i)), i);
    for (std::size_t i = 0; i < m_scales.size(); ++ i)
      add(m_features.add<Elevation> (m_range, m_point_map, grid(i), radius_dtm(i)), i);
  }

  void generate_normal_based_features(const CGAL::Default_property_map<face_iterator, typename Geom_traits::Vector_3>&)
  {
    for (std::size_t i = 0; i < m_scales.size(); ++ i)
      add(m_features.add<Verticality> (m_range, eigen(i)), i);
  }

  template <typename T>
  const T& get_parameter (const T& t)
  {
    return t;
  }

  template <typename T>
  Default_property_map<face_iterator, T>
  get_parameter (const Default&)
  {
    return Default_property_map<face_iterator, T>();
  }

  void generate_features_impl (std::size_t nb_scales)
  {
    CGAL::Real_timer t; t.start();
    
    m_scales.reserve (nb_scales);
    
    float voxel_size = -1;
    
    m_scales.push_back (new Scale (m_input, m_range, m_point_map, m_bbox, voxel_size, 0));
    voxel_size = m_scales[0]->grid_resolution();
    for (std::size_t i = 1; i < nb_scales; ++ i)
    {
      voxel_size *= 2;
      m_scales.push_back (new Scale (m_input, m_range, m_point_map, m_bbox, voxel_size, i, m_scales[i-1]->grid));
    }
    t.stop();
    CGAL_CLASSIFICATION_CERR << "Scales computed in " << t.time() << " second(s)" << std::endl;
    t.reset();

    t.start();
    
#ifdef CGAL_LINKED_WITH_TBB
    m_features.begin_parallel_additions();
#endif
    
    generate_point_based_features ();
    generate_normal_based_features (CGAL::Default_property_map<face_iterator, typename Geom_traits::Vector_3>());

#ifdef CGAL_LINKED_WITH_TBB
    m_features.end_parallel_additions();
#endif
    
    for (std::size_t i = 0; i < m_features_to_rename.size(); ++ i)
      m_features_to_rename[i].first->set_name
        (m_features_to_rename[i].first->name() + "_" + std::to_string(m_features_to_rename[i].second));
    
    t.stop();
    CGAL_CLASSIFICATION_CERR << "Features computed in " << t.time() << " second(s)" << std::endl;
  }

};


} // namespace Classification
  
} // namespace CGAL



#endif // CGAL_CLASSIFICATION_MESH_FEATURE_GENERATOR_H
