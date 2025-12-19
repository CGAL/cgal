// Copyright (c) 1997-2013 INRIA Sophia-Antipolis (France).
// All rights reserved.

namespace CGAL
{

/*!
\ingroup PkgPeriodic2Triangulation2VertexFaceClasses

The class `Periodic_2_triangulation_vertex_base_2` is a model
of the concept `Periodic_2TriangulationVertexBase_2` to be used by
`Triangulation_data_structure_2` to represent vertices of a
periodic triangulation.

The first one `Traits` is the geometric traits, it is to be
instantiated by a model of the concept
`Periodic_2TriangulationTraits_2`. The second argument is the base
class to which the additional information for the periodic vertex is
added and should be a model of `TriangulationDSVertexBase_2`


\cgalModels{Periodic_2TriangulationVertexBase_2}

\sa `CGAL::Periodic_2_triangulation_face_base_2`
\sa `CGAL::Triangulation_vertex_base_2`
\sa `CGAL::Triangulation_vertex_base_with_info_2`

*/
template< >
class Periodic_2_triangulation_vertex_base_2
{
public:

}; /* end Periodic_2_triangulation_vertex_base_2 */
} /* end namespace CGAL */
