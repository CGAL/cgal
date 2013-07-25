/*!
\ingroup PkgMesh_3SecondaryConcepts
\cgalConcept

The concept `LabeledImage_3` describes the requirements for the second template 
parameter of the class `CGAL::Labeled_image_mesh_domain_3<Image,BGT>` 
which represents mesh domains defined by 3D labeled images. A 3D labeled 
image is a 3D array of elements of an integral 
type `Type`. `Type` can be `bool`, `char`, `short`, 
`int`, or `long` (signed or not). Such an array is 
associated to a 3D axis-aligned regular grid, in 
\f$ \mathbb{R}^3\f$. A cell of this grid is denoted by <I>voxel</I>. A voxel is 
an iso-cuboid of size `vx()`, `vy()`, and `vz()`. 

\cgalHasModel `CGAL::Image_3<Kernel, T>`, for any \cgal kernel `K` and any integral type `T`

*/

class LabeledImage_3 {
public:

/// \name Types 
/// @{

/*!
Type of voxel data. Must be an 
integral type. 
*/ 
typedef unspecified_type Type; 

/*!
Ring number type. 
*/ 
typedef unspecified_type RT; 

/// @} 

/// \name Operations 
/// @{

/*!
First dimension of the 3D array, 
i.e., the number of voxels along the x coordinate axis. 
*/ 
int xdim(); 

/*!
Second dimension of the 3D array. 
*/ 
int ydim(); 

/*!
Third dimension of the 3D array. 
*/ 
int zdim(); 

/*!
Size of each voxel along x coordinate axis. 
*/ 
RT vx(); 

/*!
Size of each voxel along y coordinate axis. 
*/ 
RT vy(); 

/*!
Size of each voxel along z coordinate axis. 
*/ 
RT vz(); 

/*!
Pointer to the first element of the 3D 
image. The size of the array must be `xdim()`\f$ \times\f$`ydim()`\f$ \times\f$`zdim()`. 
*/ 
const Type* data(); 

/// @}

}; /* end LabeledImage_3 */
