


/*!
\ingroup PkgMesh2Concepts
\cgalConcept


The concept `MeshingCriteria_2` defines the meshing criteria to be used in the
algorithm. It provides a predicate `Is_bad` that tests a triangle
according to criteria. The return type of `Is_bad` is an enum
`Mesh_2::Face_badness`.

The possible values of `CGAL::Mesh_2::Face_badness` are `CGAL::Mesh_2::NOT_BAD`,
`CGAL::Mesh_2::BAD` and `CGAL::Mesh_2::IMPERATIVELY_BAD`. If the predicate returns `BAD`,
the triangle is marked as bad and the algorithm will try to destroy it. If
the predicates returns `CGAL::Mesh_2::IMPERATIVELY_BAD`, the algorithm will destroy
the triangle unconditionally during its execution.

The termination of the algorithm is guaranteed when criteria are shape
criteria corresponding to a bound on smallest angles not less than
\f$ 20.7\f$ degrees (this corresponds to a radius-edge ratio bound not less
than \f$ \sqrt{2}\f$). Any size criteria that are satisfied by small enough
triangle can be added to the set of criteria without compromising
the termination.

Note that, in the presence of input angles smaller than \f$ 60\f$ degrees,
some bad shaped triangles can appear in the final mesh in the
neighboring of small angles. To achieve termination and the respect of
size criteria everywhere, the `Is_bad` predicate has to return
`CGAL::Mesh_2::IMPERATIVELY_BAD` when size criteria are not satisfied, and
`CGAL::Mesh_2::BAD` when shape criteria are not satisfied.

`MeshingCriteria_2` also provides a type `Quality` designed to code a quality
measure for triangles. The type `Quality` must be <I>less-than
comparable</I> as the meshing algorithm will order bad triangles by quality,
to split those with smallest quality first. The predicate `Is_bad`
computes the quality of the triangle as a by-product.

\cgalHasModel `CGAL::Delaunay_mesh_criteria_2<CDT>`
\cgalHasModel `CGAL::Delaunay_mesh_size_criteria_2<CDT>`


*/

class MeshingCriteria_2 {
public:


/// \name Types
/// @{

/*!
Handle to a face of the triangulation.
*/
typedef unspecified_type Face_handle;


/*!
Default constructible, copy constructible,
assignable, and less-than comparable type.
*/
typedef unspecified_type Quality;


/*!
Predicate object. Must provide two operators.

- The first operator `Mesh_2::Face_badness operator()(Face_handle fh,
  Quality& q)` returns
   - `CGAL::Mesh_2::NOT_BAD` if it satisfies the desired
     criteria for mesh triangles,
   - `CGAL::Mesh_2::BAD` if it does not, and
   - `CGAL::Mesh_2::IMPERATIVELY_BAD` if it does not and should be refined
     unconditionally.

  In addition, this operator assigns to `q` a value measuring the quality
  of the triangle pointed by `fh`.

- The second operator `Mesh_2::Face_badness operator()(Quality q)` returns
  - `CGAL::Mesh_2::NOT_BAD` if `q` is the quality of a good triangle,
  - `CGAL::Mesh_2::BAD` if the `q` represents a poor quality, and
  - `CGAL::Mesh_2::IMPERATIVELY_BAD` if `q` represents the quality of a bad
  triangle that should be refined unconditionally.
*/
typedef unspecified_type Is_bad;

/// @}


/// \name Access to predicate and constructor objects
/// @{

/*!

*/
Is_bad is_bad_object();





/// @}

}; /* end MeshingCriteria_2 */

