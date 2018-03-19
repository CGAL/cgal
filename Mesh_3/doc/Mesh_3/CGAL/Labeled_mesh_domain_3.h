namespace CGAL {

/*!
\ingroup PkgMesh_3Domains

\brief The class `Labeled_mesh_domain_3` implements indexed domains.

This class is a model of concept `MeshDomain_3`.

Any boundary facet is labeled <a,b>, a<b, where a and b are the
tags of its incident subdomain.
Thus, a boundary facet of the domain is labeled <0,b>, where b!=0.

This class includes a <em>labeling function</em> that provides, the subdomain index of any
query point. An intersection between a segment and bounding
surfaces is detected when both segment endpoints are associated with different
values of subdomain indices. The intersection is then constructed by bisection.
The bisection stops when the query segment is shorter than an error bound
`e` given by the product of the
length of the diagonal of the bounding box (in world coordinates), or the radius of the bounding sphere, and
a relative error bound passed as argument to the constructor of `Labeled_mesh_domain_3`.

Implicit_mesh_domain_3 is a heir of Labeled_mesh_domain_3. It uses a satisfying labeling function if there is only one component to mesh.

\tparam BGT is a geometric traits class that provides
the basic operations to implement
intersection tests and intersection computations
through a bisection method. This parameter must be instantiated
with a model of the concept `BisectionGeometricTraits_3`.

\cgalHeading{Labeling function}

A labeling function `f` must return `0` if the point isn't located in any subdomain. The return type of labeling functions is an integer.

Let `p` be a Point.
<ul>
<li>`f(p)=0` means that `p` is outside domain.</li>
<li>`f(p)=a`, `a!=0` means that `p` is inside subdomain `a`.</li>
</ul>
`CGAL::Implicit_multi_domain_to_labeling_function_wrapper` is a good candidate for this template parameter
if there are several components to mesh.

The function type can be any model of the concept `Callable` compatible with the signature `Subdomain_index(const Point_3&)`: it can be a function, a function object, a lambda expression... that takes a `%Point_3` as argument, and returns a type convertible to `Subdomain_index`.

\cgalModels MeshDomain_3

\sa `Implicit_mesh_domain_3`
\sa `Implicit_multi_domain_to_labeling_function_wrapper`
\sa `CGAL::make_mesh_3()`.

*/
template<typename BGT>
class Labeled_mesh_domain_3
{
public:

/// \name Types
///@{

/// The subdomain index of this model of `MeshDomain_3`.
typedef int Subdomain_index;

/// The type of that stores the function internally
typedef CGAL::cpp11::function<Subdomain_index(const Point_3&)> Labeling_function;

///@}
/// \name Types imported from the geometric traits class
///@{

/// The point type of the geometric traits class
typedef typename Geom_traits::Point_3      Point_3;
/// The sphere type of the geometric traits class
typedef typename Geom_traits::Sphere_3     Sphere_3;
/// The iso-cuboid type of the geometric traits class
typedef typename Geom_traits::Iso_cuboid_3 Iso_cuboid_3;
/// The number type (a field type) of the geometric traits class
typedef typename Geom_traits::FT           FT;
///@}

/// \name Creation
/// @{
/*!  \brief Construction from a function, a bounding
object and a relative error bound.

This constructor uses named parameters (from the <em>Boost Parameter
Library</em>). They can be specified in any order.

\cgalHeading{Named Parameters}
- <b>`parameters::function` (mandatory)</b>  the labeling function, compatible with `Labeling_function`.
- <b>`parameters::bounding_object` (mandatory)</b>  the bounding object is either a bounding sphere (of type `Sphere_3`, a bounding box (type `Bbox_3`), or a bounding `Iso_cuboid_3`. It bounds the meshable space.
- <b>`parameters::relative_error_bound` (optional)</b>  the relative error bound used to compute intersection points between the implicit surface and query segments. The
bisection is stopped when the length of the intersected segment is less than the product of `relative_error_bound` by the diameter of the bounding object. Its default value is `FT(1e-3)`.

\cgalHeading{Example}
From the example (\ref Mesh_3/mesh_implicit_domains_2.cpp):
\snippet Mesh_3/mesh_implicit_domains_2.cpp Domain creation

 */
template <typename ... A_i>
Labeled_mesh_domain_3(const A_i&...);

///@}

/// \name Creation of domains from 3D images


/*!
\brief Construction from a 3D gray image

This static method is a <em>named constructor</em>. It constructs a domain
described by a 3D gray image. A 3D gray image is a grid of voxels,
where each voxel is associated with a gray level value.  Unless otherwise specified by the parameter `image_values_to_subdom_indices`, the domain to
be discretized is the union of voxels that lie inside a surface
described by an isolevel value, called \a isovalue. The voxels lying
inside the domain have gray level values that are larger than the
isovalue.

The value of voxels is interpolated to  a gray level value at any query point.

This constructor uses named parameters (from the <em>Boost Parameter
Library</em>). They can be specified in any order.

\cgalHeading{Named Parameters}
The parameters are optional unless otherwise specified.
<ul>

<li> <b>`parameters::image` (mandatory)</b> the input 3D image. Must
be a `CGAL::Image_3` object.

<li><b>`parameters::iso_value`</b> the isovalue, inside
  `image`, of the surface describing the boundary of the object to be
  meshed. Its default value is `0`.

<li><b>`parameters::image_values_to_subdom_indices`</b> a function or
  a function object, compatible with the signature
  `Subdomain_index(word_type)`, where the type `word_type` is the type
  of words of the 3D image. This function returns the subdomain index
  corresponding to a pixel value. If this parameter is used, then the
  parameter `iso_value` is ignored.

<li><b>`parameter::value_outside`</b> the value attached to voxels
 outside of the domain to be meshed. It should be lower than
 `iso_value`. Its default value is `0`.

<li><b>`parameter::relative_error_bound`</b> is the relative error
  bound, relative to the diameter of the box of the image. Its default
  value is `FT(1e-3)`.  </ul>

\cgalHeading{Examples}

From the example (\ref Mesh_3/mesh_3D_gray_image.cpp), where the name
of the parameters is not specified, as they are given is the same
order as the parameters definition:

\snippet Mesh_3/mesh_3D_gray_image.cpp Domain creation

From the example (\ref Mesh_3/mesh_3D_gray_vtk_image.cpp):

\snippet Mesh_3/mesh_3D_gray_vtk_image.cpp Domain creation

 */
template <typename ... A_i>
static
Labeled_mesh_domain_3
Labeled_mesh_domain_3::create_gray_image_mesh_domain(A_i&...);

/*!
\brief Construction from a 3D labeled image

This static method is a <em>named constructor</em>. It constructs a
domain described by a 3D labeled image. A 3D labeled image is a grid
of voxels, where each voxel is associated with an index (a subdomain
index) characterizing the subdomain in which the voxel lies. The
domain to be discretized is the union of voxels that have an non-zero
values.

This constructor uses named parameters (from the <em>Boost Parameter
Library</em>). They can be specified in any order.

\cgalHeading{Named Parameters}
The parameters are optional unless otherwise specified.
<ul>

<li> <b>`parameters::image` (mandatory)</b> the input 3D image. Must
be a `CGAL::Image_3` object.

<li><b>`parameter::value_outside`</b> the value attached to voxels
 outside of the domain to be meshed. Its default value is `0`.

<li><b>`parameter::relative_error_bound`</b> is the relative error
  bound, relative to the diameter of the box of the image. Its default
  value is `FT(1e-3)`.  </ul>

\cgalHeading{Example}

From the example (\ref Mesh_3/mesh_3D_image.cpp):

\snippet Mesh_3/mesh_3D_image.cpp Domain creation

 */
template <typename ... A_i>
static
Labeled_mesh_domain_3
Labeled_mesh_domain_3::create_labeled_image_mesh_domain(A_i&...);

/// \name Deprecated constructors
///
/// Those three constructors have been deprecated since CGAL-4.13, and
/// replaced by the constructor using the <em>Boost Parameter Library</em>.
///
/// @{

/*!
\brief Construction from a labeling function, a bounding Sphere and a relative error bound.
\param f the labeling function.
\param bounding_sphere the bounding sphere of the meshable space.
\param relative_error_bound is the relative error bound used to compute intersection points between the implicit surface and query segments. The
bisection is stopped when the length of the intersected segment is less than the product of `relative_error_bound` by the radius of
`bounding_sphere`.
\deprecated This constructor is deprecated since CGAL-4.13, and
replaced by the constructor using the <em>Boost Parameter Library</em>.
*/ 
Labeled_mesh_domain_3(Labeling_function f,
                       const Sphere_3& bounding_sphere,
                       const FT& relative_error_bound = FT(1e-3));

/*!
\brief Construction from a labeling function, a bounding box and a relative error bound.
\param f the labeling function.
\param bbox the bounding box of the meshable space.
\param relative_error_bound is the relative error bound used to compute intersection points between the implicit surface and query segments. The
bisection is stopped when the length of the intersected segment is less than the product of `relative_error_bound` by the diagonal of
`bounding_box`.
\deprecated This constructor is deprecated since CGAL-4.13, and
replaced by the constructor using the <em>Boost Parameter Library</em>.
*/
Labeled_mesh_domain_3(Labeling_function f,
                       const Bbox_3& bbox,
                       const FT& relative_error_bound = FT(1e-3));

/*!
\brief Construction from a function, a bounding Iso_cuboid_3 and a relative error bound.
\param f the function.
\param bbox the bounding box of the meshable space.
\param relative_error_bound is the relative error bound used to compute intersection points between the implicit surface and query segments. The
bisection is stopped when the length of the intersected segment is less than the product of `relative_error_bound` by the diagonal of
`bounding_box`.
\deprecated This constructor is deprecated since CGAL-4.13, and
replaced by the constructor using the <em>Boost Parameter Library</em>.
*/
Labeled_mesh_domain_3(Labeling_function f,
                       const Iso_cuboid_3& bbox,
                       const FT& relative_error_bound = FT(1e-3));

/// @}

}; /* end Labeled_mesh_domain_3 */
} /* end namespace CGAL */
