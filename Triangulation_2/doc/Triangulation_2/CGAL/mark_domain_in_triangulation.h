namespace CGAL {

/*!
\ingroup PkgTriangulation2Miscellaneous

\brief marks faces connected with non constrained edges
as inside of the domain based on the nesting level.

The function starts from facets incident to the infinite vertex, with a nesting
level of 0. Then recursively considers the non-explored facets incident
to constrained edges bounding the former set and increase the nesting level by 1.
Facets in the domain are those with an odd nesting level.

For this overload, the "in domain" information is set in `ipm`.

\tparam CT a constrained triangulation
\tparam InDomainPmap a class model of `ReadWritePropertyMap` with `CT::Face_handle` as key type and `bool` as value type.

\pre `ct.dimension() == 2`
*/
template <typename CT, typename InDomainPmap>
void
mark_domain_in_triangulation(CT& ct, InDomainPmap ipm);



 /*!
\ingroup PkgTriangulation2Miscellaneous

\brief marks faces connected with non constrained edges
as inside of the domain based on the nesting level.

The function starts from facets incident to the infinite vertex, with a nesting
level of 0. Then recursively considers the non-explored facets incident
to constrained edges bounding the former set and increase the nesting level by 1.
Facets in the domain are those with an odd nesting level.

For this overload, the "in domain" marker is contained in the face type of `CT`, which must provide
the methods `bool is_in_domain()`  and `void set_in_domain(bool)`.

\tparam CT a constrained triangulation

\pre `ct.dimension() == 2`
*/
template <typename CT>
void
mark_domain_in_triangulation(CT& ct);

} /* namespace CGAL */
