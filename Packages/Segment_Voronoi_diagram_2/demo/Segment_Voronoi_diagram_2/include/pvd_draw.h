#ifndef SVD_DRAW_H
#define SVD_DRAW_H



template<class T, class Widget>
void draw_diagram(Widget& widget, const T& svd)
{
  widget << CGAL::BLUE;
#if !defined (__POWERPC__)
  widget << CGAL::PointSize(3);
  widget << CGAL::LineWidth(3);
#endif

  typename T::Finite_edges_iterator eit = svd.finite_edges_begin();
  for (; eit != svd.finite_edges_end(); ++eit) {
    if ( eit->first->vertex( svd.cw(eit->second) )->info() !=
	 eit->first->vertex( svd.ccw(eit->second) )->info() ) {
      svd.draw_dual_edge(*eit, widget);	
    }
#if 0
    Site_2 p = eit->first->vertex(  cw(eit->second) )->site();
    Site_2 q = eit->first->vertex( ccw(eit->second) )->site();

    bool is_endpoint_of_seg =
      ( p.is_segment() && q.is_point() &&
	is_endpoint_of_segment(q, p) ) ||
      ( p.is_point() && q.is_segment() &&
	is_endpoint_of_segment(p, q) );

    if ( !is_endpoint_of_seg ) {
      svd.draw_dual_edge(*eit, widget);	
    }
#endif
  }
}




#endif // SVD_DRAW_H
