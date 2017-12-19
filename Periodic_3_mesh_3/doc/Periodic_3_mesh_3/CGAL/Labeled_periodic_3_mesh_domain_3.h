namespace CGAL {

/*!
\ingroup PkgPeriodic_3_mesh_3Domains

\brief The class `Labeled_periodic_3_mesh_domain_3` implements indexed periodic domains.

This class is a model of concept `Periodic_3MeshDomain_3`.

Any boundary facet is labeled `<a,b>`, with `a<b`, where `a` and `b` are the
tags of its incident subdomains.
Thus, a boundary facet of the domain is labeled `<0,b>`, where `b!=0`.

This class includes a function that provides the index of the subdomain
in which any query point lies. An intersection between a segment and bounding
surfaces is detected when both segment endpoints are associated with different
values of subdomain indices. The intersection is then constructed by bisection.
The bisection stops when the query segment is shorter than an error bound `e`
given by the product of the length of the diagonal of the cuboid
(in world coordinates) and a relative error bound passed as argument
to the constructor of `Labeled_periodic_3_mesh_domain_3`.

`Implicit_periodic_3_mesh_domain_3` is derived from `Labeled_periodic_3_mesh_domain_3`.
It uses a satisfactory labeling function if there is only one component to mesh.

\tparam LabelingFunction is the type of the input function.<br />
        This parameter stands for a model of the concept `ImplicitFunction`
        described in the surface mesh generation package.<br />
        The function is defined over the canonical representation of the three dimensional flat torus,
        the input canonical cube (see Section \ref Periodic_3_mesh_3InputDomain).
        The labeling function `f` must return `0` if the point isn't located in any subdomain.
        Usually, the return types of labeling functions are integer.<br />
        Let `p` be a point.
        <ul>
        <li>`f(p)=0` means that `p` is outside domain.</li>
        <li>`f(p)=a`, `a!=0` means that `p` is inside the subdomain `a`.</li>
        </ul>
        `Implicit_multi_domain_to_labeling_function_wrapper` is a good candidate for this template parameter
        if there are several components to mesh.

\tparam BGT is a geometric traits class that provides
        the basic operations to implement
        intersection tests and intersection computations
        through a bisection method. This parameter must be instantiated
        with a model of the concept `BisectionGeometricTraits_3`.

\cgalModels Periodic_3MeshDomain_3

\sa `CGAL::make_periodic_3_mesh_3()`.
\sa `Implicit_periodic_3_mesh_domain_3`

\sa `Labeled_mesh_domain_3`
\sa `Implicit_mesh_domain_3`

*/
template<class LabelingFunction, class BGT>
class Labeled_periodic_3_mesh_domain_3
{
public:
  typedef typename BGT::Iso_cuboid_3          Iso_cuboid_3;

/// \name Creation
/// @{

/*!
\brief Construction from a function, a fundamental domain, and a relative error bound.
\param f is the function.
\param cuboid is the fundamental domain in which we construct the mesh.
\param relative_error_bound is the relative error bound used to compute intersection points between the implicit surface and query segments. The
       bisection is stopped when the length of the intersected segment is less
       than the product of `relative_error_bound` by the diagonal of `cuboid`.
*/
Labeled_periodic_3_mesh_domain_3(const LabelingFunction& f,
                                 const BGT::Iso_cuboid_3& cuboid,
                                 const BGT::FT& relative_error_bound = FT(1e-6));

/// @}

/*!
 * \brief Return the user-chosen cube that is the canonical instance of the flat torus.
*/
const Iso_cuboid_3& canonical_periodic_domain();

}; /* end Labeled_periodic_3_mesh_domain_3 */
} /* end namespace CGAL */
