
/*!
\ingroup PkgSTLExtensionConcepts
\cgalConcept

The concept `SurjectiveLockDataStructure` describes functionalities
of a data structure whose goal is to lock
objects in a multi-threaded environment.
Such data structures are intended to be used by concurrent
algorithms.

Note that it is allowed to \"lock too much\". E.g., the data structure
might be a voxel grid and locking a point might be locking the
voxel containing this point.
Thus, a point may also be locked because
another point in the same voxel has been locked.
The only requirement is that when a thread owns the lock of an object, no other
thread can lock the same object.

We call `S` the surjective function such that `S(object)`
is the \"thing\" that is locked when one tries to lock `object`.
In the previous example, `S(point)` is the voxel containing `point`.

\cgalHasModelsBegin
\cgalHasModels{CGAL::Spatial_lock_grid_3}
\cgalHasModelsEnd

*/
class SurjectiveLockDataStructure {
public:

/// \name Operations
/// @{
  /// Test if `object` is locked (by this thread or by any other thread).
  template <typename T>
  bool is_locked(const T &object);

  /// Test if `object` is locked by this thread.
  template <typename T>
  bool is_locked_by_this_thread(const T &object);

  /// Try to lock `object`. Returns `true` if the object is already locked by this thread or if the object could be locked.
  template <typename T>
  bool try_lock(const T &object);

  /// Try to lock `object`. Returns `true` if the object is already locked by this thread or if the object could be locked.
  /// \tparam no_spin If `true`, force non-blocking operation (in any case, the
  ///                 function will return immediately, i.e.\ it will not
  ///                 wait for the resource to be free).
  ///                 If `false`, use the default behavior (same as previous function).
  template <bool no_spin, typename T>
  bool try_lock(const T &object);

  /// Unlock everything that is locked by this thread.
  void unlock_everything_locked_by_this_thread();

  /// Unlock everything that is locked by this thread except `S(object)`.
  template <typename T>
  void unlock_everything_locked_by_this_thread_but_one(const T &object);
/// @}
};
