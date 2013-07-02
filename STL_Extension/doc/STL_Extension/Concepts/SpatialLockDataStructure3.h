
/*!
\ingroup PkgStlExtensionConcepts
\cgalConcept

The concept `SpatialLockDataStructure_3` is intended to be used by concurrent
algorithms. It allows to lock points (x, y, z) in a 3D bounding box.
Each template type called `P3` below must provide x(), y() and z() functions,
returning the respective point coordinates as numbers whose type is a 
model of the concept of `CGAL::RealEmbeddable`.

Note that it is allowed to "lock too much". E.g., the embedded data structure 
may be a voxel grid and locking a point may be locking the voxel containing this point.
The only requirement is that when a thread owns the lock of a point, no other
thread can lock a point with the exact same coordinates.

\cgalHasModel `CGAL::Spatial_grid_lock_data_structure_3`

*/
class SpatialLockDataStructure_3 {
public:
  
/// \name Operations 
/// @{ 
  /// Set the bounding box of the domain
  /// Depending on the way one manages the 3D domain, this may be an
  /// empty function.
  void set_bbox(const CGAL::Bbox_3 &bbox);
  
  /// Test if `point` is locked (by this thread or by any other thread).
  template <typename P3>
  bool is_locked(const P3 &point);

  /// Test if `point` is locked by this thread.
  template <typename P3>
  bool is_locked_by_this_thread(const P3 &point);

  /// Try to lock a point in the domain. Returns true if the point is already locked by this thread or if the point could be locked.
  template <typename P3>
  bool try_lock(const P3 &point);

  /// Try to lock the point `p`. Returns true if the point is already locked by this thread or if the point could be locked.
  /// \tparam no_spin If true, force non-blocking operation (never blocks). If false, use the default behavior (same as previous function).
  template <bool no_spin, typename P3>
  bool try_lock(const P3 &p);

  /// Unlock all the locations locked by this thread
  void unlock_all_points_locked_by_this_thread();
  
  /// Unlock all the locations locked by this thread except point `p`
  template <typename P3>
  void unlock_all_tls_locked_locations_but_one_point(const P3 &p);
/// @} 
};