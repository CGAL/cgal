#ifndef DRAW_PREFERENCES_H
#define DRAW_PREFERENCES_H

#include <CGAL/IO/Window_stream.h>

#ifndef CGAL_IO_FILE_DRAWER_H
#include <CGAL/IO/Pm_drawer.h>
#endif

#ifndef CGAL_IO_DRAW_PM_H
#include <CGAL/IO/draw_pm.h>
#endif

CGAL_BEGIN_NAMESPACE

template< class Arr_2, 
          class Ccb_halfedge_circulator, 
          class Holes_iterator >
class My_Arr_drawer : public Pm_drawer< Arr_2, Window_stream >{
private:
  typedef Pm_drawer<Arr_2,Window_stream>  Base;
public:
  My_Arr_drawer( Window_stream& W ): Pm_drawer<Arr_2,Window_stream>( W ){}
  
  void draw_face(typename Base::Face_handle f) {
    if (f->does_outer_ccb_exist()) {
      Ccb_halfedge_circulator cc=f->outer_ccb();
      do {
	W << cc->curve();
      } while (++cc != f->outer_ccb());  
    }

    Holes_iterator hit=f->holes_begin(),eit=f->holes_end();
    for (;hit!=eit; ++hit) {
      Ccb_halfedge_circulator cc=*hit; 
      do {
	W << cc->curve();
	} while (++cc != *hit);  
    }      
  }

  void draw_vertices(typename Base::Vertex_const_iterator Vertices_begin, 
		     typename Base::Vertex_const_iterator Vertices_end) {
    W << GREEN;
    Base::draw_vertices(Vertices_begin, Vertices_end);
  }
  
  void draw_halfedges(typename Base::Halfedge_const_iterator Halfedges_begin, 
		      typename Base::Halfedge_const_iterator Halfedges_end) {
    W << BLUE;
    Base base(window());
    base.draw_halfedges(Halfedges_begin, Halfedges_end);
  }
  
};
 
CGAL_END_NAMESPACE
#endif
