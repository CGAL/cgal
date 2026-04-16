/*!
\ingroup PkgFeatureGraphConcepts
\cgalConcept

The concept 'FeatureImage_3' describes a 3D image that contains a label for each voxel.
It also stores the dimension, voxel size as well as the position in 3D space of the voxel at index (0, 0, 0).

\cgalHasModelsBegin
\cgalHasModels{CGAL::Image_3}
\cgalHasModelsEnd

*/

class FeatureImage_3 {
public:

/// \name Access Functions
/// @{

/*!
Returns the value of the voxel at the index (i, j, k)
*/
float value(const std::size_t i, const std::size_t j, const std::size_t k) const;

/*!
Returns the X dimension of the image.
*/
std::size_t xdim() const;
/*!
Returns the Y dimension of the image.
*/
std::size_t ydim() const;
/*!
Returns the Z dimension of the image.
*/
std::size_t zdim() const;

/*!
Returns the voxel size in the X dimension
*/
double vx() const;
/*!
Returns the voxel size in the Y dimension
*/
double vy() const;
/*!
Returns the voxel size in the Z dimension
*/
double vz() const;

/*!
Returns the position in 3D space along the X axis of the voxel at index (0, 0, 0).
*/
float tx() const;
/*!
Returns the position in 3D space along the Y axis of the voxel at index (0, 0, 0).
*/
float ty() const;
/*!
Returns the position in 3D space along the Z axis of the voxel at index (0, 0, 0).
*/
float tz() const;

/// @}

}; /* end FeatureImage_3 */