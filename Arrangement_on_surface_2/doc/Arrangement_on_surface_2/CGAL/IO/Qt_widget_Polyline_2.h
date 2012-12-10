namespace CGAL {

/*!
\ingroup PkgArrangement2op_left_shift
Inserts a polyline into a given `Qt_widget` stream.
Only the basic geometric and topological features of the polylines are
written. Auxiliary data that might be attached is lost.
*/
template<typename SegmentTraits>
Qt_widget& operator<< (std::ostream& os,
                       const Polyline_2<SegmentTraits>& polyline);
}
