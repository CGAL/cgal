#ifndef OVERLAY_COLOR_TRAITS_H
#define OVERLAY_COLOR_TRAITS_H

template <typename Arrangement> struct Overlay_color_traits {
  typedef unsigned int                                  Color;
  typedef typename Arrangement::Vertex_const_handle     V_const_handle;
  typedef typename Arrangement::Halfedge_const_handle   H_const_handle;
  typedef typename Arrangement::Face_const_handle       F_const_handle;
  typedef typename Arrangement::Vertex_handle           V_handle;
  typedef typename Arrangement::Halfedge_handle         H_handle;
  typedef typename Arrangement::Face_handle             F_handle;

  // Compute the average of the red, green, and blue components separately.
  Color blend(Color color1, Color color2) const
  {
    return
      ((((color1 & 0x000000ff) + (color2 & 0x000000ff)) / 2) & 0x000000ff) |
      ((((color1 & 0x0000ff00) + (color2 & 0x0000ff00)) / 2) & 0x0000ff00) |
      ((((color1 & 0x00ff0000) + (color2 & 0x00ff0000)) / 2) & 0x00ff0000);
  }

  void create_face(F_const_handle f1, F_const_handle f2, F_handle f) const
  { f->set_data(blend(f1->data(), f2->data())); }
  void create_vertex(H_const_handle h1, H_const_handle h2, V_handle v) const
  { v->set_data(blend(h1->data(), h2->data())); }
  void create_vertex(V_const_handle v1, V_const_handle v2, V_handle v) const
  { v->set_data(blend(v1->data(), v2->data())); }
  void create_vertex(V_const_handle v1, H_const_handle h2, V_handle v) const
  { v->set_data(blend(v1->data(), h2->data())); }
  void create_vertex(H_const_handle h1, V_const_handle v2, V_handle v) const
  { v->set_data(blend(h1->data(), v2->data())); }
  void create_vertex(F_const_handle f1, V_const_handle v2, V_handle v) const
  { v->set_data(blend(f1->data(), v2->data())); }
  void create_vertex(V_const_handle v1, F_const_handle f2, V_handle v) const
  { v->set_data(blend(v1->data(), f2->data())); }
  void create_edge(H_const_handle h1, H_const_handle h2, H_handle h) const
  {
    h->set_data(blend(h1->data(), h2->data()));
    h->twin()->set_data(blend(h1->data(), h2->data()));
  }
  void create_edge(H_const_handle h1, F_const_handle f2, H_handle h) const
  {
    h->set_data(blend(h1->data(), f2->data()));
    h->twin()->set_data(blend(h1->data(), f2->data()));
  }
  void create_edge(F_const_handle f1, H_const_handle h2, H_handle h) const
  {
    h->set_data(blend(f1->data(), h2->data()));
    h->twin()->set_data(blend(f1->data(), h2->data()));
  }
};

#endif
