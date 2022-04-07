// Copyright (c) 1999-2018   INRIA Nancy - Grand Est (France).
// All rights reserved.


namespace CGAL {

/*!
\ingroup PkgPeriodic4HyperbolicTriangulation2VertexFaceClasses

The class `Periodic_4_hyperbolic_triangulation_vertex_base_2` is the default model for the
concept `Periodic_4HyperbolicTriangulationVertexBase_2`. It accepts two template parameters:

\tparam GT         Geometric traits type. This should be the same model of the concept
                        `Periodic_4HyperbolicDelaunayTriangulationTraits_2` that is used in the class
                        `Periodic_4_hyperbolic_Delaunay_triangulation_2`. This template parameter has
                        no default value.
\tparam Vb         Vertex base type. Should be a model of the concept `TriangulationVertexBase_2`.
                           The default value for this template parameter is `CGAL::Triangulation_vertex_base_2<GT>`

`Periodic_4_hyperbolic_triangulation_vertex_base_2` can be simply plugged in the triangulation
data structure of a periodic hyperbolic triangulation, or used as a base class to derive other
base vertex classes tuned for specific applications.

\cgalModels `Periodic_4HyperbolicTriangulationVertexBase_2`

\sa `Periodic_4_hyperbolic_triangulation_face_base_2`
*/


template< class GT, class Vb >
class Periodic_4_hyperbolic_triangulation_vertex_base_2 : public Vb {

};

}  // namespace CGAL
