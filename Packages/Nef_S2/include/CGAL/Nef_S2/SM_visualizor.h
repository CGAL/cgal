#ifndef SM_VISUALIZOR_H
#define SM_VISUALIZOR_H

#include <CGAL/basic.h>
#include <CGAL/Nef_S2/SM_decorator.h>
#include <CGAL/Nef_S2/SM_triangulator.h>
#include <CGAL/Nef_S2/Sphere_geometry_OGL.h>

#define USING(t) typedef typename Sphere_map_::t t
//#define LGREY CGAL::Color(170,170,170)
//#define DGREY CGAL::Color(30,30,30)
#define LGREY CGAL::Color(170,170,200)
#define DGREY CGAL::Color(30,30,50)
CGAL_BEGIN_NAMESPACE

template <typename Sphere_map_>
class SM_BooleColor 
{
  typedef typename Sphere_map_::Vertex_const_handle   Vertex_const_handle;   
  typedef typename Sphere_map_::Halfedge_const_handle Halfedge_const_handle;   
  typedef typename Sphere_map_::Halfloop_const_handle Halfloop_const_handle;   
  typedef typename Sphere_map_::Face_const_handle     Face_const_handle;   
  typedef typename Sphere_map_::Mark Mark;   
public:
  Color color(Vertex_const_handle, Mark m) const
  { return ( m ? CGAL::BLACK : CGAL::WHITE ); }
  Color color(Halfedge_const_handle, Mark m) const
  { return ( m ? CGAL::BLACK : CGAL::WHITE ); }
  Color color(Halfloop_const_handle, Mark m) const
  { return ( m ? CGAL::BLACK : CGAL::WHITE ); }
  Color color(Face_const_handle, Mark m) const
  { return ( m ? DGREY : LGREY ); }
};


/*{\Moptions outfile=SM_visualizor.man }*/
/*{\Manpage {SM_visualizor}{Sphere_map_,Sphere_kernel_}
{Drawing plane maps}{V}}*/

template <typename Sphere_map_, typename Sphere_kernel_,
	  typename Color_ = SM_BooleColor<Sphere_map_> >
class SM_visualizor : public 
  SM_triangulator< SM_decorator<Sphere_map_,Sphere_kernel_> >
{
/*{\Mdefinition An instance |\Mvar| of the data type |\Mname| is a
decorator to draw the structure of a sphere map into the surface
of a OpenGL sphere. It is generic with respect to the template 
concept.}*/

/*{\Mgeneralization SM_decorator}*/
/*{\Mtypes 3}*/
public:
  typedef Sphere_map_    Sphere_map;
  typedef Sphere_kernel_ Sphere_kernel;
  typedef SM_visualizor<Sphere_map_,Sphere_kernel_,Color_> Self;
  typedef SM_const_decorator<Sphere_map_,Sphere_kernel_>   Explorer;
  typedef SM_decorator<Sphere_map_,Sphere_kernel_>         Decorator;
  typedef SM_triangulator<Decorator>                       Base;

  USING(Vertex_const_handle);   
  USING(Halfedge_const_handle); 
  USING(Face_const_handle);     
  USING(Vertex_const_iterator);
  USING(Halfedge_const_iterator);
  USING(Face_const_iterator);
  USING(Mark);

  typedef typename Sphere_kernel::Sphere_point    Sphere_point;
  typedef typename Sphere_kernel::Sphere_segment  Sphere_segment;
  typedef typename Sphere_kernel::Sphere_circle   Sphere_circle;
  typedef typename Sphere_kernel::Sphere_triangle Sphere_triangle;

  typedef Color_ Color_objects;
  /*{\Mtypemember The color data accessor.}*/  

protected:
  Explorer                E_;
  const Color_objects&    CO_;
  Sphere_map              MT_;
  CGAL::OGL::Unit_sphere& S_;
public:

/*{\Mcreation 4}*/
SM_visualizor(const Sphere_map& M, CGAL::OGL::Unit_sphere& S,
	      const Color_objects& C = Color_objects())
/*{\Mcreate creates an instance |\Mvar| of type |\Mname| to visualize
    the vertices, edges, and faces of |D| in an open GL window.}*/
  : Base(M,MT_), E_(M), CO_(C), MT_(), S_(S)
  { triangulate(); }

/*{\Moperations 2 1}*/

/* |draw_map| draws all object of the referenced sphere map:
   1) edges, loops, and vertices are taken from E_
   2) faces are drawn via the calculated triangulation in MT_  */

void draw_map() const
/*{\Mop draw the whole plane map.}*/
{
  // draw sphere segments underlying edges of E_:
  Halfedge_const_iterator e;
  CGAL_forall_edges(e,E_) {
    if ( source(e) == target(e) ) {
      S_.push_back(E_.circle(e), CO_.color(e,E_.mark(e))); 
    } else {
      S_.push_back(Sphere_segment(E_.point(E_.source(e)),
				  E_.point(E_.target(e)),
				  E_.circle(e)),CO_.color(e,E_.mark(e)));
    }
  }
  // draw sphere circles underlying loops of E_:

  if ( E_.has_loop() )
    S_.push_back(
      Sphere_circle(E_.circle(E_.halfloop())),
      CO_.color(E_.halfloop(),E_.mark(E_.halfloop())));

  // draw points underlying vertices of E_:
  Vertex_const_iterator v;
  CGAL_forall_vertices(v,E_)
    S_.push_back(E_.point(v),CO_.color(v,E_.mark(v)));

  Unique_hash_map<Halfedge_const_iterator,bool> Done(false);
  CGAL_forall_halfedges(e,*this) {
    if ( Done[e] ) continue;
    Halfedge_const_handle en(next(e)),enn(next(en));
    TRACEV(Base::incident_triangle(e));
    TRACEN(incident_mark(e)<<incident_mark(en)<<incident_mark(enn));
    CGAL_nef_assertion(Base::incident_mark(e)==Base::incident_mark(en) &&
		   Base::incident_mark(en)==Base::incident_mark(enn));
    Mark m = Base::incident_mark(e);
    Sphere_triangle t = Base::incident_triangle(e);
    S_.push_back(t, (m ? DGREY : LGREY) );
    Done[e]=Done[en]=Done[enn]=true;
  }

  Done.clear(false);
  CGAL_forall_halfedges(e,*this) {
    if ( Done[e] ) continue;
    S_.push_back_triangle_edge(Sphere_segment(E_.point(E_.source(e)),
					      E_.point(E_.target(e)),
					      E_.circle(e)));
    Done[e]=Done[twin(e)]=true;
  }

}

/* |draw_triangulation| draws all object of the underlying triangulation:
   1) edges, loops, and vertices are taken from E_
   2) faces are drawn via the calculated triangulation in MT_  */

void draw_triangulation() const
{ 
  // draw sphere segments underlying edges of triangulation:
  Halfedge_const_iterator e;
  CGAL_forall_edges(e,*this) {
    S_.push_back(Sphere_segment(point(source(e)),point(target(e)),
				circle(e)),CO_.color(e,mark(e)));
  }

  // draw points underlying vertices of triangulation:
  Vertex_const_iterator v;
  CGAL_forall_vertices(v,*this)
    S_.push_back(point(v),CO_.color(v,mark(v)));

  Done.clear(false);
  CGAL_forall_halfedges(e,*this) {
    if ( Done[e] ) continue;
    Halfedge_const_handle en(next(e)),enn(next(en));
    CGAL_nef_assertion(incident_mark(e)==incident_mark(en)&&
		   incident_mark(en)==incident_mark(enn));
    Mark m = incident_mark(e);
    Sphere_triangle t = incident_triangle(e);
    S_.push_back(t, (m ? DGREY : LGREY) );
    Done[e]=Done[en]=Done[enn]=true;
  }


}

}; // end of SM_visualizor 




CGAL_END_NAMESPACE
#undef USING
//#undef LGREY
//#undef DGREY
#endif // SM_VISUALIZOR_H

