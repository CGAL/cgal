// Copyright (c) 2018 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_CLASSIFICATION_PROPERTY_MAPS_H
#define CGAL_CLASSIFICATION_PROPERTY_MAPS_H

#include <CGAL/license/Classification.h>

#include <CGAL/centroid.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/property_map.h>

#include <CGAL/boost/graph/properties.h>

namespace CGAL
{

namespace Classification
{

/*!
  \ingroup PkgClassificationMesh

  \brief Property map that constructs the center of mass of the face
  of a mesh on-the-fly.

  \cgalModels `ReadablePropertyMap`

  \tparam FaceGraph model of `FaceGraph`.

  \tparam VertexPointMap model of `ReadablePropertyMap` with with
  `boost::graph_traits<FaceGraph>::%vertex_descriptor` as key type
  and `CGAL::Point_3` as value type.
*/
template <typename FaceGraph,
          typename VertexPointMap = typename boost::property_map<FaceGraph,vertex_point_t>::type >
class Face_descriptor_to_center_of_mass_map
{
public:
  typedef typename boost::property_traits<VertexPointMap>::value_type Point_3;
  typedef typename boost::graph_traits<FaceGraph>::face_descriptor key_type;
  typedef Point_3 value_type;
  typedef Point_3 reference;
  typedef boost::readable_property_map_tag category;

private:
  typedef typename boost::graph_traits<FaceGraph>::vertex_descriptor vertex_descriptor;

  const FaceGraph* m_mesh;
  VertexPointMap m_vpm;

public:

  Face_descriptor_to_center_of_mass_map ()
    : m_mesh (nullptr) { }
  Face_descriptor_to_center_of_mass_map (const FaceGraph* mesh)
    : m_mesh (mesh), m_vpm (get (vertex_point, *m_mesh)) { }
  Face_descriptor_to_center_of_mass_map (const FaceGraph* mesh, VertexPointMap vpm)
    : m_mesh (mesh), m_vpm (vpm) { }

  /// \cond SKIP_IN_MANUAL
  inline friend reference get (const Face_descriptor_to_center_of_mass_map& map, key_type f)
  {
    std::vector<Point_3> points;

    for(vertex_descriptor v : vertices_around_face(halfedge(f, *(map.m_mesh)), *(map.m_mesh)))
      points.push_back (get (map.m_vpm, v));

    return CGAL::centroid (points.begin(), points.end());
  }
  /// \endcond
};

/*!
  \ingroup PkgClassificationMesh

  \brief Property map that constructs a face descriptor with a
  `bbox()` method from a face descriptor.

  \cgalModels `ReadablePropertyMap`

  \tparam FaceGraph model of `FaceGraph`.

  \tparam VertexPointMap model of `ReadablePropertyMap` with with
  `boost::graph_traits<FaceGraph>::%vertex_descriptor` as key type
  and `CGAL::Point_3` as value type.
*/
template <typename FaceGraph,
          typename VertexPointMap = typename boost::property_map<FaceGraph,vertex_point_t>::type >
class Face_descriptor_to_face_descriptor_with_bbox_map
{
public:
  typedef typename boost::graph_traits<FaceGraph>::face_descriptor face_descriptor;

  /*!
    \brief Face descriptor with a precomputed bounding box.
  */
  class face_descriptor_with_bbox
  {
    face_descriptor m_descriptor;
    CGAL::Bbox_3 m_bbox;

  public:
    face_descriptor_with_bbox (const face_descriptor& descriptor,
                               const CGAL::Bbox_3& bbox)
      : m_descriptor (descriptor), m_bbox (bbox)
    { }

    const CGAL::Bbox_3 bbox() const { return m_bbox; }
    operator face_descriptor() const { return m_descriptor; }
  };

  typedef face_descriptor key_type;
  typedef face_descriptor_with_bbox value_type;
  typedef face_descriptor_with_bbox reference;
  typedef boost::readable_property_map_tag category;

private:
  typedef typename boost::graph_traits<FaceGraph>::vertex_descriptor vertex_descriptor;

  const FaceGraph* m_mesh;
  VertexPointMap m_vpm;

public:

  Face_descriptor_to_face_descriptor_with_bbox_map ()
    : m_mesh (nullptr) { }
  Face_descriptor_to_face_descriptor_with_bbox_map (const FaceGraph* mesh)
    : m_mesh (mesh), m_vpm (get (vertex_point, *m_mesh)) { }
  Face_descriptor_to_face_descriptor_with_bbox_map (const FaceGraph* mesh, VertexPointMap vpm)
    : m_mesh (mesh), m_vpm (vpm) { }

  /// \cond SKIP_IN_MANUAL
  inline friend reference get (const Face_descriptor_to_face_descriptor_with_bbox_map& map, key_type f)
  {
    CGAL::Bbox_3 bbox;

    for(vertex_descriptor v : vertices_around_face(halfedge(f, *(map.m_mesh)), *(map.m_mesh)))
      bbox = bbox + get(map.m_vpm, v).bbox();

    return value_type (f, bbox);
  }
  /// \endcond
};


} // namespace Classification

} // namespace CGAL

#endif // CGAL_CLASSIFICATION_PROPERTY_MAPS_H
