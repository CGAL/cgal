// Copyright (c) 1997-2013 INRIA Sophia-Antipolis (France).
// All rights reserved.

/*!
\ingroup PkgPeriodic2Triangulation2Concepts
\cgalConcept

At the base level (see Section \ref
Section_2D_Triangulations_Software_Design), a face stores handles to
its four vertices and to its four neighbor faces. The vertices and
neighbors are indexed 0, 1 and 2. Neighbor \f$ i\f$ lies opposite to
vertex \f$ i\f$.

\cgalRefines ::TriangulationFaceBase_2

\cgalHasModel CGAL::Periodic_2_triangulation_face_base_2

\sa `TriangulationDataStructure_2`
\sa `TriangulationFaceBase_2`
\sa `Periodic_2TriangulationVertexBase_2`

*/

class Periodic_2TriangulationFaceBase_2
{
public:

/// \name Access Functions
/// @{

  /*!
  Returns the offset of vertex `i`.
  \pre \f$ i \in\{0, 1, 2\}\f$.
  */
  int offset(int i) const;

  /*!
  Returns true if the offset of vertex `i` is zero for \f$ i \in\{0, 1, 2\}\f$.
  */
  bool has_zero_offsets() const;

/// @}

/// \name Setting
/// @{

  /*!
  Sets the vertex offsets according to `off0` to `off2`.
  */
  void set_offsets(int off0, int off1, int off2);

/// @}

}; /* end Periodic_2TriangulationFaceBase_2 */

