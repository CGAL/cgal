namespace CGAL {

/*!
\ingroup PkgMesh3Domains

\deprecated The class template `Implicit_mesh_domain_3` is deprecated
since CGAL-4.13, in favor of the class template `Labeled_mesh_domain_3` and
its static function
`Labeled_mesh_domain_3::create_implicit_mesh_domain()`.

The class `Implicit_mesh_domain_3` implements a domain whose bounding surface is
described
implicitly as the zero level set of a function.
The domain to be discretized is assumed to be the domain where
the function has negative values.
This class is a model of the concept `MeshDomain_3`.


\tparam Function provides the definition of the function.
This parameter stands for a model of the concept
`ImplicitFunction` described in the
surface mesh generation package.
The number types `Function::FT`
and `BGT::FT` are required to match.

\tparam BGT is a geometric traits which provides the basic operations to implement
intersection tests and computations
through a bisection method. This parameter must be instantiated
with a model of the concept `BisectionGeometricTraits_3`.

The constructor of `Implicit_mesh_domain_3`
takes as argument a bounding sphere which is required
to circumscribe the surface and to have its center inside the
domain.
This domain constructs intersection points
between
the surface and segments/rays/lines by
bisection. It needs an
`error_bound` such that the bisection process is stopped
when the query segment is smaller than the error bound.
The `error_bound` passed as argument to the domain constructor
is a relative error bound expressed as a ratio to the bounding sphere radius.

\cgalModels `MeshDomain_3`

\sa `BisectionGeometricTraits_3`
\sa `CGAL::make_mesh_3()`.

*/
template< typename Function, typename BGT >
class Implicit_mesh_domain_3
{
public:

/// \name Creation
/// @{

/*!
`f` is the object of type `Function` that represents the implicit
surface.

`bounding_sphere` is a bounding sphere of the implicit surface. The
value of `f` at the sphere center `c` must be
negative: \f$ f(c)<0\f$.

`error_bound` is the relative error bound
used to compute intersection points between the implicit surface
and query segments. The
bisection is stopped when the length of the intersected
segment is less than the product of `error_bound` by the
radius of `bounding_sphere`.
*/
  Implicit_mesh_domain_3(Function f,
                         const BGT::Sphere_3& bounding_sphere,
                         const BGT::FT& error_bound = FT(1e-3));

/// @}

}; /* end Implicit_mesh_domain_3 */
} /* end namespace CGAL */
