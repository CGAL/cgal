namespace CGAL {

/*!
\ingroup PkgPeriodic_3_mesh_3Domains

The class `Implicit_periodic_3_mesh_domain_3` implements a periodic domain
whose bounding surface is described implicitly as the zero level set
of a function defined over the three dimensional flat torus. In practice,
the domain of definition of the function must contain the input canonical cube,
see Section \ref Periodic_3_mesh_3InputDomain.
The domain to be discretized is assumed to be the domain where
the function has negative values.

This class is a model of the concept `Periodic_3MeshDomain_3`.

\tparam Function provides the definition of the function.
        This parameter stands for a model of the concept `ImplicitFunction`
        described in the surface mesh generation package.
        The number types `Function::FT` and `BGT::FT` are required to match.

\tparam BGT is a geometric traits which provides the basic operations to implement
        intersection tests and computations through a bisection method.
        This parameter must be instantiated with a model of the concept
        `BisectionGeometricTraits_3`.

The constructor of `Implicit_periodic_3_mesh_domain_3` takes as argument
a cuboid (the canonical representative of the fundamental domain) in which
we construct the mesh (see \ref PkgPeriodic_3_mesh_3).
This domain constructs intersection points between the surface and
segments (duals of a facet) by bisection. It requires an `error_bound`
such that the bisection process is stopped when the query segment is smaller
than the error bound. The `error_bound`, passed as argument to the domain constructor,
is a relative error bound expressed as a ratio of the longest space diagonal
of the cuboid.

\cgalModels `Periodic_3MeshDomain_3`

\sa `CGAL::make_periodic_3_mesh_3()`.
\sa `Labeled_periodic_3_mesh_domain_3`

\sa `BisectionGeometricTraits_3`
\sa `Implicit_mesh_domain_3`

*/
template< typename Function, typename BGT >
class Implicit_periodic_3_mesh_domain_3 {
public:

/// \name Creation
/// @{

/*!
\param f is the object of type `Function` that represents the implicit surface
\param cuboid is the fundamental domain
\param error_bound  is the relative error bound used to compute intersection
       points between the implicit surface and query segments. The bisection
       is stopped when the length of the intersected segment is less than
       the product of `error_bound` by the diagonal of `cuboid`.
*/
Implicit_periodic_3_mesh_domain_3(Function f,
                                  BGT::Iso_cuboid_3 cuboid,
                                  BGT::FT error_bound = FT(1e-6));

/// @}

}; /* end Implicit_periodic_3_mesh_domain_3 */
} /* end namespace CGAL */
