#ifdef CGAL_TRIANGULATION_2_H


CGAL_BEGIN_NAMESPACE

template <class Gt,class Tds>
PS_Stream& operator << (PS_Stream& ps, const Triangulation_2<Gt,Tds> &t)

{
 t.draw_triangulation(ps);
  return ps;
}
CGAL_END_NAMESPACE

#endif // CGAL_TRIANGULATION_2_H


#ifdef CGAL_DELAUNAY_TRIANGULATION_2_H


CGAL_BEGIN_NAMESPACE

template < class Gt, class Tds >
PS_Stream& operator << (PS_Stream& ps, const Delaunay_triangulation_2<Gt,Tds> &t)
{
 t.draw_triangulation(ps);
 return ps; 
}

template < class Gt, class Tds >
void draw_psdual(PS_Stream & ps,Delaunay_triangulation_2<Gt,Tds> &t)
{
  ps <<border_color(CGAL::Color(255,0,0));
  Delaunay_triangulation_2<Gt,Tds>::Edge_iterator eit, ebegin=t.edges_begin(), eend=t.edges_end();

 for (eit=ebegin; eit != eend; ++eit)
    {
        CGAL::Object o = t.dual(eit);
        typename Gt::Ray r;
        typename Gt::Segment s;
        if (CGAL::assign(s,o)) ps << s;
        if (CGAL::assign(r,o)) ps << r;
    }
  
} 

CGAL_END_NAMESPACE

#endif // CGAL_DELAUNAY_TRIANGULATION_2_H

#ifdef CGAL_CONSTRAINED_TRIANGULATION_2_H


CGAL_BEGIN_NAMESPACE

template < class Gt, class Tds>
PS_Stream& operator<<(PS_Stream& ps,const Constrained_triangulation_2<Gt,Tds> &t)
{

 t.draw_triangulation(ps);
 return ps;
  

}


CGAL_END_NAMESPACE

#endif // CGAL_CONSTRAINED_TRIANGULATION_2_H


#ifdef CGAL_REGULAR_TRIANGULATION_2_H


CGAL_BEGIN_NAMESPACE

template < class Gt, class Tds >
 PS_Stream& operator << (PS_Stream& ps,/*const*/ Regular_triangulation_2<Gt,Tds> &t)
 {
    ps << border_color(CGAL::Color(0,255,0));
    t.draw_triangulation(ps);
   return ps;
	

	if (t.number_of_vertices()>1 )
	{	
	  Regular_triangulation_2<Gt,Tds>::Vertex_iterator 
	    v_i=t.vertices_begin(),
	    done=t.vertices_end();

	while(v_i != done)
		{
			ps <<border_color(CGAL::Color(255,0,0)) <<(v_i->point());
			++v_i;
		}
	}else if (t.number_of_vertices()==1)
	{
		ps<<border_color(CGAL::Color(255,0,0)) <<t.finite_vertex()->point();
	}

        typedef typename Gt::Point  Point;
	typedef list<Point> Point_list;
	
	if( t.number_of_vertices() < 2) {return ps;}

	Regular_triangulation_2<Gt,Tds>::Face_iterator fit;
	for(fit=t.faces_begin(); fit != t.faces_end(); ++fit){
	  Point_list::iterator plit;
	 for(plit=fit->point_list().begin(); 
	     plit!=fit->point_list().end(); ++plit) {
	    ps <<border_color(CGAL::Color(0,0,255))<< *plit;}
	}
	Regular_triangulation_2<Gt,Tds>::Face_circulator 
	  fc = t.infinite_vertex()->incident_faces(),done(fc);
	do {
	  Point_list::iterator plit;
	  for(plit=fc->point_list().begin(); 
	     plit!=fc->point_list().end(); ++plit) {
	    ps <<border_color(CGAL::Color(0,0,255))<< *plit;}
       	   }
        while(++fc != done); 
 
        return ps;
}


CGAL_END_NAMESPACE

#endif // CGAL_REGULAR_TRIANGULATION_2_H





	




