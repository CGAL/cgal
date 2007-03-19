#ifndef CGAL_OVERLAY_FUNCTOR_H
#define CGAL_OVERLAY_FUNCTOR_H


template < class Arrangement_ >
class Overlay_functor
{
public:

  typedef Arrangement_                                  Arrangement;

  typedef typename Arrangement::Face_const_handle       Face_const_handle;
  typedef typename Arrangement::Vertex_const_handle     Vertex_const_handle;
  typedef typename Arrangement::Halfedge_const_handle   Halfedge_const_handle;

  typedef typename Arrangement::Face_handle             Face_handle;
  typedef typename Arrangement::Vertex_handle           Vertex_handle;
  typedef typename Arrangement::Halfedge_handle         Halfedge_handle;


  void create_face (Face_const_handle f1,
                    Face_const_handle f2,
                    Face_handle res_f)
  {
    const QColor& c1 = f1->color();
    const QColor& c2 = f2->color();
    QColor res_c ( (c1.red() + c2.red())/2,
                   (c1.green() + c2.green())/2,
                   (c1.blue() + c2.blue())/2 );

    res_f->set_color(res_c);
  }

  void create_vertex(Halfedge_const_handle ,
                     Halfedge_const_handle ,
                     Vertex_handle )
  {}

  void create_vertex(Vertex_const_handle ,
                     Vertex_const_handle ,
                     Vertex_handle )
  {}

  void create_vertex(Vertex_const_handle ,
                     Halfedge_const_handle ,
                     Vertex_handle )
  {}

  void create_vertex(Halfedge_const_handle ,
                     Vertex_const_handle ,
                     Vertex_handle )
  {}

  void create_vertex(Face_const_handle ,
                     Vertex_const_handle ,
                     Vertex_handle )
  {}

  void create_vertex(Vertex_const_handle ,
                     Face_const_handle ,
                     Vertex_handle )
  {}

  void create_edge(Halfedge_const_handle ,
                   Halfedge_const_handle ,
                   Halfedge_handle )
  {}

  void create_edge(Halfedge_const_handle ,
                   Face_const_handle ,
                   Halfedge_handle )
  {}

  void create_edge(Face_const_handle ,
                   Halfedge_const_handle ,
                   Halfedge_handle )
  {}


};

#endif
