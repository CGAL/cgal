// ============================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/IO/triangulation_Window_stream.h
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Olivier Devillers
//                 Andreas Fabri
//                 Monique Teillaud
//                 Mariette Yvinec
//
// coordinator   : Mariette Yvinec  <Mariette Yvinec@sophia.inria.fr>
//
// ============================================================================


#ifdef CGAL_TRIANGULATION_2_H
#ifndef CGAL_WINDOW_STREAM_TRIANGULATION_2_H
#define CGAL_WINDOW_STREAM_TRIANGULATION_2_H

CGAL_BEGIN_NAMESPACE

template < class Gt, class Tds>
Window_stream&
operator<<(Window_stream& os,
           const Triangulation_2<Gt, Tds> &T)
{
    Triangulation_2<Gt, Tds>::Edge_iterator it = T.edges_begin();

    while(it != T.edges_end()){
        os << T.segment(it);
        ++it;
    }

    return os;
}

CGAL_END_NAMESPACE

#endif // CGAL_WINDOW_STREAM_TRIANGULATION_2_H
#endif // CGAL_TRIANGULATION_2_H

#ifdef CGAL_DELAUNAY_TRIANGULATION_2_H
#ifndef CGAL_WINDOW_STREAM_DELAUNAY_TRIANGULATION_2_H
#define CGAL_WINDOW_STREAM_DELAUNAY_TRIANGULATION_2_H

CGAL_BEGIN_NAMESPACE

template < class Gt, class Tds >
Window_stream&
operator<<(Window_stream& os,
           const Delaunay_triangulation_2<Gt,Tds> &T)
{
  return (os << (Triangulation_2<Gt, Tds>) T);
  // c'est pas une bonne idee parceque le cast
  // utilise le createur 
  //Triangulation_2(const Triangulation_2<Gt,Tds> &tr)
  //qui recopie toute les faces 

//   Delaunay_triangulation_2<Gt,Tds>::Edge_iterator it = T.edges_begin();
// 
//     while(it != T.edges_end()){
//         os << T.segment(it);
//         ++it;
//     }
//   
//     return os;
}

CGAL_END_NAMESPACE

#endif // CGAL_WINDOW_STREAM_DELAUNAY_TRIANGULATION_2_H
#endif // CGAL_DELAUNAY_TRIANGULATION_2_H

#ifdef CGAL_CONSTRAINED_TRIANGULATION_2_H
#ifndef CGAL_WINDOW_STREAM_CONSTRAINED_TRIANGULATION_2_H
#define CGAL_WINDOW_STREAM_CONSTRAINED_TRIANGULATION_2_H

CGAL_BEGIN_NAMESPACE

template < class Gt, class Tds>
Window_stream&
operator<<(Window_stream& os,
           const Constrained_triangulation_2<Gt,Tds> &T)
{

  //return (os << (Triangulation_2<Gt, Tds>) T);
  
  Constrained_triangulation_2<Gt,Tds>::Edge_iterator it = T.edges_begin();

    while(it != T.edges_end()){
        os << T.segment(it);
        ++it;
    }
   return os;
}

CGAL_END_NAMESPACE

#endif // CGAL_WINDOW_STREAM_CONSTRAINED_TRIANGULATION_2_H
#endif // CGAL_CONSTRAINED_TRIANGULATION_2_H


#ifdef CGAL_WEIGHTED_POINT_H
#ifndef CGAL_WINDOW_STREAM_WEIGHTED_POINT_H
#define CGAL_WINDOW_STREAM_WEIGHTED_POINT_H

CGAL_BEGIN_NAMESPACE

template < class Point, class We >
Window_stream&
operator<<(Window_stream& os, const Weighted_point< Point, We > &p)
{
  double cx = to_double(p.point().x()),
         cy = to_double(p.point().y()),
         r = to_double(p.weight());

  os<<p.point();
  os.draw_circle(cx, cy , /*sqrt*/(r));
  return os;
}

template < class Point, class We >
Window_stream& operator>>(Window_stream &os, Weighted_point< Point, We > &wp)
{
  double cx, cy, x1, y1;
  os.read_mouse(cx,cy);
  os.read_mouse_circle(cx,cy, x1, y1);
  Point center(cx, cy);

  We sr = We (sqrt( square(cx-x1)+square(cy-y1) ) );

  os.draw_circle(cx, cy , sr );
  wp = Weighted_point< Point, We >(center, sr);
  return os;
}

CGAL_END_NAMESPACE

#endif // CGAL_WINDOW_STREAM_WEIGHTED_POINT_H
#endif // CGAL_WEIGHTED_POINT_H


#ifdef CGAL_REGULAR_TRIANGULATION_2_H
#ifndef CGAL_WINDOW_STREAM_REGULAR_TRIANGULATION_2_H
#define CGAL_WINDOW_STREAM_REGULAR_TRIANGULATION_2_H

CGAL_BEGIN_NAMESPACE

template < class Gt, class Tds >
Window_stream&
operator<<(Window_stream& os,
           /*const*/ Regular_triangulation_2<Gt,Tds> &T)
{
	Regular_triangulation_2<Gt,Tds>::Edge_iterator 
	  it = T.edges_begin();

	while(it != T.edges_end())
	  {	os << GREEN<<T.segment(it);
		++it;
	}

	if (T.number_of_vertices()>1 )
	{	
	  Regular_triangulation_2<Gt,Tds>::Vertex_iterator 
	    v_i=T.vertices_begin(),
	    done=T.vertices_end();

		while(v_i != done)
		{
			os<< RED<<(v_i->point());
			++v_i;
		}
	}else if (T.number_of_vertices()==1)
	{
		os<<RED<<T.finite_vertex()->point();
	}

	typedef typename Gt::Point    Point;
	typedef list<Point> Point_list;
	
	if( T.number_of_vertices() < 2) {return os;}

	Regular_triangulation_2<Gt,Tds>::Face_iterator fit;
	for(fit=T.faces_begin(); fit != T.faces_end(); ++fit){
	  Point_list::iterator plit;
	 for(plit=fit->point_list().begin(); 
	     plit!=fit->point_list().end(); ++plit) {
	    os<<BLUE<< *plit;}
	}
	
	Regular_triangulation_2<Gt,Tds>::Face_circulator 
	  fc = T.infinite_vertex()->incident_faces(),done(fc);
	do {
	  Point_list::iterator plit;
	  for(plit=fc->point_list().begin(); 
	     plit!=fc->point_list().end(); ++plit) {
	    os<<BLUE<< *plit;}
	}while(++fc != done); 
 
    return os;
}

CGAL_END_NAMESPACE

#endif // CGAL_WINDOW_STREAM_REGULAR_TRIANGULATION_2_H
#endif // CGAL_REGULAR_TRIANGULATION_2_H
