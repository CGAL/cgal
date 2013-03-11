
namespace CGAL {

/*!
\ingroup PkgStlExtension

The class `Spatial_grid_lock_data_structure_3` allows to lock 
point coordinates (x, y, z) in a 3D grid. For example, it can be used by 
concurrent algorithms to lock simplices.

\tparam Grid_lock_tag allows to choose the locking strategy used by the
structure. The following tags are available:
- `Tag_non_blocking_with_atomics` is non-blocking (i.e. trye_lock will always
return immediately) and uses atomics to perform lock operations.
- `Tag_non_blocking_with_mutexes` is non-blocking (i.e. trye_lock will always
return immediately) and uses mutexes to perform lock operations.
- `Tag_priority_blocking_with_atomics` provides both non-blocking and
partially-locking `try_lock` versions. The partially-locking version will 
block/spin if the thread owning the lock has a "lower ID" (each threads is 
given an arbitrary ID) or return immediately otherwise. It uses atomics to 
perform lock operations.

\cgalModels `SpatialLockDataStructure_3`
*/

template <typename Grid_lock_tag>
class Spatial_grid_lock_data_structure_3
{
public:

/// \name Creation
/// @{

/// Construct the lock grid of size `bbox`, with `num_grid_cells_per_axis` 
/// cells per axis.
Spatial_grid_lock_data_structure_3(const Bbox_3 &bbox,
                                   int num_grid_cells_per_axis);

/// @}

}; /* end Spatial_grid_lock_data_structure_3 */
} /* end namespace CGAL */
