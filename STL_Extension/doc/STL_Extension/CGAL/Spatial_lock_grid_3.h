
namespace CGAL {

/*!
\ingroup PkgStlExtension

The class `Spatial_lock_grid_3` allows to lock 
points with coordinates (x, y, z) in a 3D grid.
The point type is called `P3` here. `P3` must provide x(), y() and z() functions,
returning the respective point coordinates as numbers whose type is a 
model of the concept of `CGAL::RealEmbeddable`.

It is a model of `SurjectiveLockDataStructure`, with `T` being
`P3` and `S` being the function that maps a point to 
the cell of the 3D grid containing this point.

For example, it can be used by concurrent algorithms to lock simplices.

\tparam Grid_lock_tag allows to choose the locking strategy used by the
structure. The following tags are available:
- `Tag_non_blocking` is non-blocking (i.e.\ try_lock will always
return immediately) and uses atomics to perform lock operations.
- `Tag_priority_blocking` provides both non-blocking and
partially-blocking `try_lock` versions. The partially-blocking version will 
block (spin) if the thread owning the lock has a lower \"ID\" (each thread is 
given an arbitrary ID) or return immediately otherwise. It uses atomics to 
perform lock operations.

\cgalModels `SurjectiveLockDataStructure`
*/

template <typename Grid_lock_tag>
class Spatial_lock_grid_3
{
public:

/// \name Creation
/// @{

  /// Constructs the lock grid of size `bbox`, with `num_grid_cells_per_axis` 
  /// cells per axis.
  Spatial_lock_grid_3(const Bbox_3 &bbox,
                      int num_grid_cells_per_axis);

/// @}

/// \name Operations 
/// @{ 

  /// Sets the bounding box of the domain.
  void set_bbox(const CGAL::Bbox_3 &bbox);

/// @}

}; /* end Spatial_lock_grid_3 */
} /* end namespace CGAL */
