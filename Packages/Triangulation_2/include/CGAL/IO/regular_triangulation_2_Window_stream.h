#ifdef CGAL_REGULAR_TRIANGULATION_2_H
#ifndef CGAL_WINDOW_STREAM_REGULAR_TRIANGULATION_2_H
#define CGAL_WINDOW_STREAM_REGULAR_TRIANGULATION_2_H



template < class Gt, class Tds >
CGAL_Window_stream&
operator<<(CGAL_Window_stream& os,
           /*const*/ CGAL_Regular_triangulation_2<Gt,Tds> &T)
{
	CGAL_Regular_triangulation_2<Gt,Tds>::Edge_iterator 
	  it = T.edges_begin();

	while(it != T.edges_end())
	  {	os << CGAL_GREEN<<T.segment(it);
		++it;
	}

	if (T.number_of_vertices()>1 )
	{	
	  CGAL_Regular_triangulation_2<Gt,Tds>::Vertex_iterator 
	    v_i=T.vertices_begin(),
	    done=T.vertices_end();

		while(v_i != done)
		{
			os<< CGAL_RED<<(v_i->point());
			++v_i;
		}
	}else if (T.number_of_vertices()==1)
	{
		os<<CGAL_RED<<T.finite_vertex()->point();
	}

	typedef typename Gt::Point    Point;
	typedef list<Point> Point_list;
	
	if( T.number_of_vertices() < 2) {return os;}

	CGAL_Regular_triangulation_2<Gt,Tds>::Face_iterator fit;
	for(fit=T.faces_begin(); fit != T.faces_end(); ++fit){
	  Point_list::iterator plit;
	 for(plit=fit->point_list().begin(); 
	     plit!=fit->point_list().end(); ++plit) {
	    os<<CGAL_BLUE<< *plit;}
	}

	
	CGAL_Regular_triangulation_2<Gt,Tds>::Face_circulator 
	  fc = T.infinite_vertex()->incident_faces(),done(fc);
	do {
	  Point_list::iterator plit;
	  for(plit=fc->point_list().begin(); 
	     plit!=fc->point_list().end(); ++plit) {
	    os<<CGAL_BLUE<< *plit;}
	}while(++fc != done); 
 
    return os;
}
#endif // CGAL_WINDOW_STREAM_REGULAR_TRIANGULATION_2_H
#endif // CGAL_REGULAR_TRIANGULATION_2_H

