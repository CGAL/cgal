#ifdef CGAL_TRIANGULATION_2_H


CGAL_BEGIN_NAMESPACE

template <class Gt,class Tds>
PS_Stream& operator << (PS_Stream& ps, const Triangulation_2<Gt,Tds> &t)

{
  //t.draw_triangulation(ps);
typedef  Triangulation_2<Gt, Tds>::Edge_iterator EI ;
EI it = t.edges_begin();

      while(it != t.edges_end())
	{
        ps << CGAL::line_width(4) <<border_color(CGAL::Color(0,0,255)) 
	   << t.segment(it);
	ps << border_color(CGAL::Color(0,255,0)) << CGAL::point_size(5) <<
	  point_style(CGAL::PS_Stream::FDOT) <<
	  t.segment(it).source();
	ps << t.segment(it).target();
        ++it;
	}
       
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

#ifdef CGAL_WEIGHTED_POINT_H
#ifndef CGAL_PS_STREAM_WEIGHTED_POINT_H
#define CGAL_PS_STREAM_WEIGHTED_POINT_H
CGAL_BEGIN_NAMESPACE
template < class Point, class We >
PS_Stream&
operator<<(PS_Stream& ps, const Weighted_point< Point, We > &p)
{
typedef Circle_2< CGAL::Cartesian<double> > Circle;
 
double cx = to_double(p.point().x()),
         cy = to_double(p.point().y()),
         r = to_double(p.weight());
Point o(cx,cy);
Circle ci(o,r);
 ps <<point_style(CGAL::PS_Stream::FDOT) <<o;
 ps <<ci;
 
  return ps;
}

CGAL_END_NAMESPACE
#endif // CGAL_PS_STREAM_WEIGHTED_POINT_H
#endif // CGAL_WEIGHTED_POINT_H

#ifdef CGAL_REGULAR_TRIANGULATION_2_H


CGAL_BEGIN_NAMESPACE

template < class Gt, class Tds >
 PS_Stream& operator << (PS_Stream& ps,/*const*/ Regular_triangulation_2<Gt,Tds> &t)
 {
   //  ps << border_color(CGAL::Color(0,255,0));
    t.draw_triangulation(ps);
   
	

	// if (t.number_of_vertices()>1 )
// 	{	
// 	  Regular_triangulation_2<Gt,Tds>::Vertex_iterator 
// 	    v_i=t.vertices_begin(),
// 	    done=t.vertices_end();

// 	while(v_i != done)
// 		{
// 			ps <<border_color(CGAL::Color(255,0,0)) <<(v_i->point());
// 			++v_i;
// 		}
// 	}else if (t.number_of_vertices()==1)
// 	{
// 		ps<<border_color(CGAL::Color(255,0,0)) <<t.finite_vertex()->point();
// 	}

        typedef typename Gt::Point  Point;
	typedef list<Point> Point_list;
	
	if( t.number_of_vertices() <= 2) {return ps;}

	typedef Regular_triangulation_2<Gt,Tds>::Face_iterator FI;
	FI fit;
	for(fit=t.faces_begin(); fit != t.faces_end(); ++fit) {
	  Point_list::iterator plit;
	 for(plit=fit->point_list().begin(); 
	     plit!=fit->point_list().end(); ++plit) {
	   ps << point_style(CGAL::PS_Stream::FDOT) <<border_color(CGAL::Color(0,255,0)) << (*plit).point();}
	}

	 typedef Regular_triangulation_2<Gt,Tds>::Face_circulator FC;
	 FC fc;
	 FC fc_end;
	  fc = t.infinite_vertex()->incident_faces();
	  fc_end= fc;
	do {
	  Point_list::iterator plit;
	  for(plit=fc->point_list().begin(); 
	     plit!=fc->point_list().end(); ++plit) {
	    ps<< point_style(CGAL::PS_Stream::FDOT) <<border_color(CGAL::Color(0,255,0))<< (*plit).point();}
	}while(++fc != fc_end); 
 
    return ps;
}


CGAL_END_NAMESPACE

#endif // CGAL_REGULAR_TRIANGULATION_2_H





	




