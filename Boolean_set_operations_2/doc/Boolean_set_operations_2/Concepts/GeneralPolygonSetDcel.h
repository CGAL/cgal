
/*!
\ingroup PkgBooleanSetOperations2Concepts
\cgalConcept

A doubly-connected edge-list (\em Dcel for short) data-structure. It consists 
of three containers of records: vertices \f$ V\f$, halfedges \f$ E\f$, and faces \f$ F\f$. 
It maintains the incidence relation among them. The halfedges are ordered 
in pairs sometimes referred to as twins, such that each halfedge pair 
represent an edge. 

A model of the `GeneralPolygonSetDcel` simply refines `ArrangementDcel`,
the `Halfedge` and `Face` types being models of the concepts
`GeneralPolygonSetDcelHalfedge` and `GeneralPolygonSetDcelFace`
respectively

\cgalRefines `ArrangementDcel`

\cgalHasModel `CGAL::Gps_default_dcel<Traits>` 

\sa `GeneralPolygonSetDcelFace` 
\sa `GeneralPolygonSetDcelHalfedge` 

*/

class GeneralPolygonSetDcel {};

/* end GeneralPolygonSetDcel */

