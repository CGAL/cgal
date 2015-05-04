namespace CGAL {

/*!
\ingroup PkgConvexHull3Functions

\deprecated This function relies on a package deprecated since \cgal 4.6 and thus is also deprecated.

computes the convex hull polyhedron 
of the three-dimensional points in the range [`first`,`beyond`)
and assigns it to `P`. If `test_correctness` is set to
`true`, the tests described in \cgalCite{mnssssu-cgpvg-96} are
used to determine the correctness of the resulting polyhedron.


\attention This function is provided for completeness and educational 
purposes. When an efficient incremental implementation is needed, 
using `Delaunay_triangulation_3` together with 
`convex_hull_3_to_polyhedron_3()` is highly recommended. 

\tparam InputIterator  must be an input iterator where the value type must be  
`Polyhedron::Traits::Point`.

\tparam Polyhedron must provide a type `Polyhedron::Traits` 
that defines the following types 
- `Polyhedron::Traits::R`, which is a model of 
the representation class `R` required by 
`Convex_hull_d_traits_3<R>` 
- `Polyhedron::Traits::Point`  


\sa `convex_hull_3()` 
\sa `Convex_hull_d` 

\cgalHeading{Implementation}

This function uses the `d`-dimensional convex hull incremental construction 
algorithm \cgalCite{cms-frric-93} 
with `d` fixed to 3. The algorithm requires \f$ O(n^2)\f$ time in the 
worst case and \f$ O(n \log n)\f$ expected time. 

*/
template <class InputIterator, class Polyhedron>
void
convex_hull_incremental_3(InputIterator first, 
InputIterator beyond, 
Polyhedron& P, 
bool test_correctness = false);

} /* namespace CGAL */
