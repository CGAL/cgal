namespace CGAL {

/*!
\ingroup PkgSegmentDelaunayGraph2Ref

The class `Segment_Delaunay_graph_storage_traits_with_info_2` provides a model for the
`SegmentDelaunayGraphStorageTraits_2` concept.


\tparam Gt must be a model of the concept `SegmentDelaunayGraphTraits_2`.

\cgalModels{SegmentDelaunayGraphStorageTraits_2}

\sa `CGAL::Segment_Delaunay_graph_storage_site_with_info_2`
*/
  template<class Gt, class Info, class Converter, class Merger>
class Segment_Delaunay_graph_storage_traits_with_info_2
    : public Segment_Delaunay_graph_storage_traits_2<Gt>
{
public:
  /*!
  */
  typedef Info_                                        Info;
  typedef Converter                                    Convert_info;
  typedef Merger                                       Merge_info;
  typedef typename Base::Geom_traits                   Geom_traits;

  typedef
  Segment_Delaunay_graph_storage_site_with_info_2<Self,
                                                  Info,
                                                  Base_storage_site_2> Storage_site_2;

};

} //namespace CGAL
