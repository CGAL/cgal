
namespace CGAL {

/*!
\ingroup PkgSegmentDelaunayGraph2Ref

The class `Segment_Delaunay_graph_storage_site_with_info_2` is a model for the concept
`SegmentDelaunayGraphStorageSite_2`.

\tparam STraits must be a model of the `SegmentDelaunayGraphStorageTraits_2` concept.
\tparam Base must be a model of `SegmentDelaunayGraphStorageSite_2`

\cgalModels{SegmentDelaunayGraphStorageSite_2}

\sa `CGAL::Segment_Delaunay_graph_site_2<K>`
*/
template< typename STraits, typename Info_, typename Base >
class Segment_Delaunay_graph_storage_site_with_info_2 : public Base {
public:

  /// \name Types
  /// The class
  /// @{

  typedef STraits Storage_traits;
  /*!
   */
  typedef Info_ Info;

  typedef typename Storage_traits::Geom_traits    Geom_traits;
  typedef typename Geom_traits::Site_2            Site_2;
  typedef typename Storage_traits::Point_handle   Point_handle;

  struct Has_info_tag {};
  /// @}

  /*!
   */
  const Info& info() const;

  /*!
   */
  void  set_info(const Info&) const;

  /*!
   */
  bool info_has_been_set() const;

}; /* end Segment_Delaunay_graph_storage_site_with_info_2 */
} /* end namespace CGAL */
