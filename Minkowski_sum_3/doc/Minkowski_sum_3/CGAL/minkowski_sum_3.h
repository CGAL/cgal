namespace CGAL {

/*!
\ingroup PkgMinkowskiSum3
\anchor refminkowski_sum_3 

The function `minkowski_sum_3` computes the Minkowski sum of two 
given 3D Nef polyhedra \f$ N0\f$ and \f$ N1\f$. Note that the function runs in 
\f$ O(n^3m^3)\f$ time in the worst case, where \f$ n\f$ and 
\f$ m\f$ are the complexities of the two input polyhedra (the complexity of 
a `Nef_polyhedron_3` is the sum of its `Vertices`, 
`Halfedges` and `SHalfedges`). 

An input polyhedron may consist of: 
<OL> 
<LI>singular vertices 
<LI>singular edges 
<LI>singular convex facets without holes 
<LI>surfaces with convex facets that have no holes. 
<LI>three-dimensional features, whose coplanar facets have 
common selection marks (this includes open and closed solids) 
</OL> 

Taking a different viewpoint, the implementation is restricted as 
follows: 
<OL> 
<LI>The input polyhedra must be bounded (selected outer volume is ignored). 
<LI>All sets of coplanar facets of a full-dimensional 
feature must have the same selection mark (in case of different 
selection marks, unselected is assumed). 
<LI>All facets of lower-dimensional features need to be convex and 
must not have holes (non-convex facets and holes are ignored). 
</OL> 

\post If either of the input polyhedra is non-convex, it is modified during the computation, i.e., it is decomposed into convex pieces. 

\sa `CGAL::Nef_polyhedron_3<Traits>` 
\sa `CGAL::convex_decomposition_3` 

*/
Nef_polyhedron_3 minkowski_sum_3(Nef_polyhedron_3 N0, Nef_polyhedron_3 N1);

} /* namespace CGAL */
