// Copyright (c) 1997-2013 INRIA Sophia-Antipolis (France).
// All rights reserved.

namespace CGAL
{

/*!
\ingroup PkgPeriodic2Triangulation2VertexFaceClasses

The class `Periodic_2_triangulation_face_base_2` is a model of
the concept `Periodic_2TriangulationFaceBase_2` to be used by
`Triangulation_data_structure_2` to represent faces of a periodic
triangulation.

The first one `Traits` is the geometric traits, it is to be
instantiated by a model of the concept
`Periodic_2TriangulationTraits_2`. The second argument is the base
class to which the additional information for the periodic vertex is
added and should be a model of `TriangulationDSFaceBase_2`

As faces cannot span more than one domain per direction of space in a
periodic Delaunay triangulation, it is enough to store offsets in the
range \f$ \{0,1\}^2\f$. For optimization purposes we encode all three
offsets in one integer. Each offset needs two bits (for the offset in
the \f$ x\f$- and \f$ y\f$-direction). If we number the bits from <I>least
significant to most significant</I> then bits \f$ 2*i\f$ and \f$ 2*i+1\f$ contain
the offset corresponding to vertex \f$ i\f$.

The implementation of `has_zero_offsets()` results in checking
whether all offsets are zero.

\cgalModels ::Periodic_2TriangulationFaceBase_2

\sa `CGAL::Triangulation_face_base_2`
\sa `CGAL::Triangulation_face_base_with_info_2`

*/
template< typename Gt, typename Fb = Triangulation_face_base_2<Gt> >
class Periodic_2_triangulation_face_base_2
{
public:

/// @}

}; /* end Periodic_2_triangulation_face_base_2 */
} /* end namespace CGAL */
