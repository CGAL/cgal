namespace CGAL {

/*!
\ingroup PkgMesh_3Domains

The class `Gray_image_mesh_domain_3` implements a domain described by a 3D
gray image. A 3D gray image is a grid of voxels, 
where each voxel is associated with a gray level value.
This class is a model of the concept `MeshDomain_3`.
The domain to be discretized is the union of voxels that lie inside a surface
described by an isolevel value, called \a isovalue. The voxels lying inside the
domain have gray level values that are larger than the isovalue.

This class includes a member function that provides, by interpolation,
a gray level value at any query point.
An intersection between a segment and bounding 
surfaces is detected when both segment endpoints are associated with gray level
values which are on both sides of the isovalue.
The intersection is then constructed by bisection. 
The bisection stops when the query segment is shorter than a given error bound 
`e`. This error bound is given by `e=d`\f$ \times\f$`bound` where `d` is the 
length of the diagonal of the bounding box (in world coordinates) and 
`bound` is the argument passed to the constructor of `Labeled_image_mesh_domain_3`.


\tparam Image is the type of the input image. 
This parameter must be a model of the concept 
`LabeledImage_3`. 

\tparam BGT is a geometric traits class which provides 
the basic operations to implement 
intersection tests and intersection computations 
through a bisection method. This parameter must be instantiated 
with a model of the concept `BisectionGeometricTraits_3`. 

\tparam Image_word_type is the data type encoded in the `Image`
input file


\cgalModels `MeshDomain_3`

\sa `BisectionGeometricTraits_3` 
\sa `CGAL::make_mesh_3()`. 

*/
template<typename Image,
         typename BGT,
         typename Image_word_type>
class Gray_image_mesh_domain_3
{
public:

/// \name Creation 
/// @{

/*!
Construction from an image.
The object to be meshed is described by the voxels that have a gray-level
value higher than the input isovalue.
@param image the input image
@param iso_value the isovalue, inside `image`,
       of the surface describing the boundary of the object to be meshed.
@param value_outside the value attached to voxels outside of the domain
       to be meshed. It should be lower than `iso_value`
@param error_bound is relative to the size of the image. 
*/ 
  Gray_image_mesh_domain_3(
      const Image& image,
      const Image_word_type iso_value,
      const Image_word_type value_outside = 0.,
      const BGT::FT& error_bound = BGT::FT(1e-3));

/// @}

}; /* end Labeled_image_mesh_domain_3 */
} /* end namespace CGAL */
