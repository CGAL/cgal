// Copyright (c) 2026 Tel-Aviv University (Israel). All rights reserved.
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
// Author(s): Efi Fogel         <efifogel@gmail.com>

namespace CGAL {
namespace Arrangement_on_curve_1 {

/*! \ingroup PkgArrangementOnCurve1Concepts
 * \cgalConcept
 *
 * \cgalHasModelsBegin
 * \cgalHasModels{CGAL::Arrangement_on_curve_1::Arrangement_on_curve_1<GeometryTraits_1, TopologyTraits>}
 * \cgalHasModelsEnd
 *
 * A model of the concept `ArrangementOnCurve_1` can be used to represent a 1D
 * subdivision of a continuous geometric curve (a 1D "master curve" or line
 * space) into alternating vertices and edges.
 */
class ArrangementOnCurve_1 {
public:
  /// \name Types
  /// @{

  //! the size type (convertible to `size_t`).
  typedef unspecified_type Size;

  //! The 1D arrangement point representation.
  typedef unspecified_type Point_1;

  //! Descriptor targeting a mutable vertex element.
  typedef unspecified_type Vertex_descriptor;

  //! Descriptor targeting a mutable edge element.
  typedef unspecified_type Edge_descriptor;

  //! Descriptor targeting an immutable vertex element.
  typedef unspecified_type Vertex_const_descriptor;

  //! Descriptor targeting an immutable edge element.
  typedef unspecified_type Edge_const_descriptor;

  /*! Represents a point localization result within the 1D curve structure.
   * A coordinate query point must match a `Vertex_descriptor` if it rests on an arrangement node,
   * or an `Edge_descriptor` if it falls securely between two nodes or out inside an infinity frontier.
   */
  typedef std::variant<Vertex_descriptor, Edge_descriptor> Location_result;
  typedef std::variant<Vertex_const_descriptor, Edge_const_descriptor> Const_location_result;
  /// @}

  /// \name Creation
  /// @{

  /*! Default Constructor. Allocates a new default-constructed instance on the heap.
   */
  ArrangementOnCurve_1();

  /// @}

  /// \name Accessors
  /// @{

  /*! obtains `true` if the arrangement does not contains any vertices.
   */
  bool empty() const;

  /*! obtains the total count of vertices resting along the line track.
   */
  Size number_of_vertices() const;

  /*! obtains the total number of bounded and unbounded edge sections.
   */
  Size number_of_edges() const;

  /*! obtains an iterator range tracking all constant vertex descriptors.
   */
  auto vertices() const;

  /*! obtains an iterator range tracking all constant edge descriptors.
   */
  auto edges() const;

  /*! obtains an lvalue property map matching vertex descriptors to their coordinate positions.
   */
  auto vertex_point_map() const;

  /*! obtains an lvalue property map matching vertex descriptors to their extended user attributes.
   */
  auto vertex_data_map() const;

  /*! obtains an lvalue property map matching edge descriptors to their extended user attributes.
   */
  auto edge_data_map() const;

  /*! obtains a handle to the leftmost unbounded edge spanning \f$(-\infty, v_{first})\f$.
   */
  Edge_descriptor unbounded_left_edge();

  /*! obtains a constant handle to the leftmost unbounded edge.
   */
  Edge_const_descriptor unbounded_left_edge() const;

  /*! obtains a handle to the rightmost unbounded edge spanning \f$(v_{last}, +\infty)\f$.
   */
  Edge_descriptor unbounded_right_edge();

  /*! obtains a constant handle to the rightmost unbounded edge.
   */
  Edge_const_descriptor unbounded_right_edge() const;

  /// @}

  /// \name Modification Modifiers
  /// @{

  /*! safety resets function allowing empty arrangements to bind to a separate existing traits memory frame.
   * \pre `is_empty() == true`
   */
  void reset_shared_geometry_traits(Shared_geometry_traits new_shared_traits_traits);

  /// @}
};

} // namespace Arrangement_on_curve_1
} // namespace CGAL
