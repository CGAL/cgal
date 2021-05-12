namespace CGAL {

/*!
\ingroup PkgMesh3Domains

\deprecated The class template `Labeled_image_mesh_domain_3` is deprecated
since CGAL-4.13, in favor of the class template `Labeled_mesh_domain_3` and
its static function
`Labeled_mesh_domain_3::create_labeled_image_mesh_domain()`.

The class `Labeled_image_mesh_domain_3` implements a domain described by a 3D labeled image. A 3D
labeled image is a grid of voxels, where each voxel is associated with an index
(a subdomain index) characterizing the subdomain in which the voxel lies. This
class is a model of the concept `MeshDomain_3`. The domain to be discretized
is the union of voxels that have an non-default index (different from the
default constructed value of the type `Image::Type`).

This class includes a member function that provides, by interpolation, the index
of the subdomain in which any query point lies. An intersection between a segment and bounding
surfaces is detected when both segment endpoints are associated with different
values of subdomain indices. The intersection is then constructed by bisection.
The bisection stops when the query segment is shorter than a given error bound
`e`. This error bound is given by `e=d`\f$ \times\f$`bound` where `d` is the
length of the diagonal of the bounding box (in world coordinates) and
`bound` is the argument passed to the constructor of `Labeled_image_mesh_domain_3`.


\tparam Image is the type of the input image.
This parameter must be `CGAL::Image_3`.

\tparam BGT is a geometric traits class which provides
the basic operations to implement
intersection tests and intersection computations
through a bisection method. This parameter must be instantiated
with a model of the concept `BisectionGeometricTraits_3`.

\cgalModels `MeshDomain_3`

An executable that uses `Labeled_image_mesh_domain_3` must be linked with
the <I>CGAL_ImageIO</I> library.

\sa `BisectionGeometricTraits_3`
\sa `CGAL::make_mesh_3()`.

*/
template< typename Image, typename BGT >
class Labeled_image_mesh_domain_3 {
public:

/// \name Creation
/// @{

/*!
Construction from an image.
The parameter `error_bound` is relative to the size of the image.
*/
  Labeled_Image_mesh_domain_3(const Image& image,
                              const BGT::FT& error_bound = FT(1e-3));

/// @}

}; /* end Labeled_image_mesh_domain_3 */
} /* end namespace CGAL */
