/*!
\ingroup PkgArrangement2op_left_shift
Inserts a conic arc into a given `Qt_widget` stream.
Only the basic geometric and topological features of the conic arcs are
written. Auxiliary data that might be attached is lost.
*/
template<typename ConicArc>
Qt_widget& operator<< (std::ostream& os,
                       const _Conic_x_monotone_arc_2<ConicArc>& cv);
