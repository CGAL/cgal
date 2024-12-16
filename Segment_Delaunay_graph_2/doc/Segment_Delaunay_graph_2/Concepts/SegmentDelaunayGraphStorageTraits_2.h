
/*!
\ingroup PkgSegmentDelaunayGraph2Concepts
\cgalConcept

The concept `SegmentDelaunayGraphStorageTraits_2` provides the
requirements for the storage traits of a segment Delaunay graph.

To avoid redundancy in the storage of points, the input points of a segment Delaunay graph
are stored in a container, and the various types of sites (input points and segments,
points of intersection, subsegments with one or two points of intersection as endpoints)
only store handles to the points in the container.

See section \ref Segment_Delaunay_graph_2StronglyIntersecting for more information.

\cgalRefines{CopyConstructible,Assignable,DefaultConstructible}

\cgalHasModelsBegin
\cgalHasModels{CGAL::Segment_Delaunay_graph_storage_traits_2<K>}
\cgalHasModelsEnd

\sa `SegmentDelaunayGraphTraits_2`
*/

class SegmentDelaunayGraphStorageTraits_2 {
public:

/// \name Types
/// @{

  /*!
  The geometric traits type. It must be a model of `SegmentDelaunayGraphTraits_2`.
  */
  typedef unspecified_type Geom_traits;

  /*!
  A container of unique points, used to associate a unique handle to each
  unique input geometric position.
  */
  typedef std::set<typename Geom_traits::Point_2> Point_container;

  /*!
  */
  typedef Point_container::iterator Point_handle;

  /*!
  */
  typedef Point_container::const_iterator const_Point_handle;

  /*!
  Type for the storage site. It must be a model of `SegmentDelaunayGraphStorageSite_2`.
  */
  typedef unspecified_type Storage_site_2;

  /*!
  Type of the storage site construction functor.
  */
  typedef CGAL::SegmentDelaunayGraph_2::Construct_storage_site_2<Self> Construct_storage_site_2;

/// @}

/// \name Creation
/// @{

  /*!
  Constructor.
  */
  SegmentDelaunayGraphStorageTraits_2(const Geom_traits& gt = Geom_traits());

/// @}

/// \name Access Functions
/// @{

  /*!
  returns the geometric traits.
  */
  const Geom_traits& geom_traits() const;

/// @}

/// \name Constructions
/// @{

  /*!
  returns a functor of type `Construct_storage_site_2`.
  */
  Construct_storage_site_2 construct_storage_site_2_object() const
/// @}

}; /* end SegmentDelaunayGraphStorageTraits_2 */

