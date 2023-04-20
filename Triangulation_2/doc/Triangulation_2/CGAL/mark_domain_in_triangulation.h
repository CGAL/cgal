namespace CGAL {

/*!
\ingroup PkgTriangulation2Miscellaneous
marks faces connected with non constrained edges
as inside of the domain based on the nesting level.

The function starts from facets incident to the infinite vertex, with a nesting
level of 0. Then recursively considers the non-explored facets incident
to constrained edges bounding the former set and increase the nesting level by 1.
Facets in the domain are those with an odd nesting level.

\tparam CT a constrained triangulation
\tparam InDomainPmap a class model of `ReadWritePropertyMap` with `CT::Face_handle` as key type and `bool` as value type.
*/
template <typename CT, typename InDomainPmap>
void
mark_domain_in_triangulation(CT& ct, InDomainPmap ipm);



 /*!
\ingroup PkgTriangulation2Miscellaneous
As the function above.

Requires that the face type of `CT` has the methods `bool is_in_domain()`  and `void set_in_domain(bool)`.

\tparam CT a constrained triangulation
*/
template <typename CT>
void
mark_domain_in_triangulation(CT& ct);

} /* namespace CGAL */
