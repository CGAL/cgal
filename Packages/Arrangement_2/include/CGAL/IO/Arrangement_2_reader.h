#ifndef CGAL_IO_ARRANGEMENT_2_READER_H
#define CGAL_IO_ARRANGEMENT_2_READER_H

#include <CGAL/IO/Arrangement_2_ascii_formatter.h>
#include <CGAL/Inverse_index.h>
#include <CGAL/circulator.h>
#include <algorithm>
#include <iostream>

CGAL_BEGIN_NAMESPACE

template <class Arrangement_2>
class Arrangement_2_reader
{
public: 
  typedef typename Arrangement_2::Size                    Size;
  typedef typename Arrangement_2::Vertex_iterator         Vertex_iterator;
  typedef typename Arrangement_2::Halfedge_iterator       Halfedge_iterator;
  typedef typename Arrangement_2::Face_iterator           Face_iterator;

  typedef typename Arrangement_2::Vertex_handle           Vertex_handle;
  typedef typename Arrangement_2::Halfedge_handle         Halfedge_handle;
  typedef typename Arrangement_2::Face_handle             Face_handle;

  typedef typename Arrangement_2::Vertex_const_iterator   Vertex_const_iterator;
  typedef typename Arrangement_2::Halfedge_const_iterator Halfedge_const_iterator;
  typedef typename Arrangement_2::Face_const_iterator     Face_const_iterator;

  typedef typename Arrangement_2::Vertex_const_handle     Vertex_const_handle;
  typedef typename Arrangement_2::Halfedge_const_handle   Halfedge_const_handle;
  typedef typename Arrangement_2::Face_const_handle       Face_const_handle;

  typedef typename Arrangement_2::Holes_const_iterator    Holes_const_iterator;
  typedef typename Arrangement_2::Ccb_halfedge_const_circulator
                                              Ccb_halfedge_const_circulator;
  typedef typename Arrangement_2::Isolated_vertices_const_iterator
                                              Isolated_vertices_const_iterator;

  typedef typename Arrangement_2::Dcel                    Dcel;  

  typedef typename Arrangement_2::Traits_2                Traits_2;
  typedef typename Traits_2::X_monotone_curve_2           X_monotone_curve_2;
  typedef typename Traits_2::Point_2                      Point_2;

protected:
  typedef typename Dcel::Vertex                           D_vertex;
  typedef typename Dcel::Halfedge                         D_halfedge;
  typedef typename Dcel::Face                             D_face;

  typedef typename Arrangement_2::Stored_point_2          Stored_point_2;
  typedef typename Arrangement_2::Stored_curve_2          Stored_curve_2;

public:

  Arrangement_2_reader(Arrangement_2& arr) :
    m_arr(arr)
  { }

  virtual ~Arrangement_2_reader()
  {
  }
 
  template <class Formatter_>
  void operator()(Formatter_& formatter)
  {
    formatter.read_arr_begin();

    Size number_of_vertices = 0,
         number_of_halfedges = 0,
         number_of_faces = 0;
    number_of_vertices = formatter.read_size("number_of_vertices");
    number_of_halfedges = formatter.read_size("number_of_halfedges");
    number_of_faces = formatter.read_size("number_of_faces");

    Size i;
    // read the vertices, and save in vertices vector
    formatter.read_vertices_begin();
    for(i=0; i<number_of_vertices; ++i)
    {
      D_vertex* v = read_vertex(formatter);
      m_vertices.push_back(v);
    }
    formatter.read_vertices_end();

    // read the haledges and save in the halfedges vector
    formatter.read_halfedges_begin();
    for(i=0; i<number_of_halfedges; ++i, ++i)
    {
      D_halfedge *h = read_edge(formatter);
      
      m_halfedges.push_back(h);
      m_halfedges.push_back(h->opposite());
    }
    formatter.read_halfedges_end();

    // read faces
    formatter.read_faces_begin();
    for(i=0; i<number_of_faces; ++i)
      read_face(formatter, i);
    formatter.read_faces_end();

    formatter.read_arr_end();
  }

protected:

  template <class Formatter_>
  D_vertex* read_vertex(Formatter_& formatter)
  {
    formatter.read_vertex_begin();
    D_vertex* new_vertex = m_arr.dcel.new_vertex();
    Point_2 p;
    formatter.read_vertex_point(p);
    Stored_point_2 *p_p =m_arr._new_point (p);
    new_vertex->set_point(p_p);
    formatter.read_vertex(new_vertex);
    formatter.read_vertex_end();
    return new_vertex;
  }
 
  template <class Formatter_>
  D_halfedge *read_edge(Formatter_& formatter)
  {
    formatter.read_edge_begin();
    D_halfedge *new_h = m_arr.dcel.new_edge();
    std::size_t src_idx = formatter.read_halfedge_endpoint_index("endpoint1");
    std::size_t trg_idx = formatter.read_halfedge_endpoint_index("endpoint2");
    X_monotone_curve_2 cv;
    formatter.read_halfedge_curve(cv);
    Stored_curve_2  *p_cv = m_arr._new_curve (cv);
    new_h->set_curve(p_cv);

    // connect vertex & halfedge
    m_vertices[trg_idx]->set_halfedge(new_h);
    new_h->set_vertex(m_vertices[trg_idx]);
    // take care for twin
    m_vertices[src_idx]->set_halfedge(new_h->opposite());
    new_h->opposite()->set_vertex(m_vertices[src_idx]);
   
    formatter.read_halfedge(new_h);
    formatter.read_halfedge(new_h->opposite());
    formatter.read_edge_end();
    return new_h;
  }

  template <class Formatter_>
  void read_face(Formatter_& formatter, int idx)
  {
    formatter.read_face_begin();
    // read the outer ccb:
    formatter.read_outer_ccb_begin();
    int outer_size = formatter.read_number();
    D_face* new_f;
    if (outer_size == 0)
      new_f = m_arr.un_face;  // unbounded face
    else
    {
      new_f = m_arr.dcel.new_face();
      D_halfedge* first_halfedge = read_ccb_halfedges(formatter, new_f,
                                                      outer_size);
      new_f->set_halfedge(first_halfedge);
    }
    formatter.read_outer_ccb_end();

    // read the holes of the face
    formatter.read_holes_begin();
    std::size_t k, n_holes = formatter.read_number();
    for (k=0; k<n_holes; ++k)
    {
      formatter.read_inner_ccb_begin();
      std::size_t inner_size = formatter.read_number();

      D_halfedge* first_halfedge = read_ccb_halfedges(formatter, new_f,
                                                      inner_size);
      new_f->add_hole(first_halfedge);
   
      formatter.read_inner_ccb_end();
    }

    formatter.read_holes_end();

    // read the isolated vertices of the face
    formatter.read_isolated_vertices_begin();
    std::size_t n_isolated = formatter.read_number();
    std::size_t index;
    for(k=0; k<n_isolated; ++k)
    {
      formatter.read_vertex_index(index);
      D_vertex* v = m_vertices[index];
      v->set_halfedge(new_f->halfedge());
      new_f->add_isolated_vertex(v);
    }
    formatter.read_isolated_vertices_end();

    // this is used for writing face additional data
    formatter.read_face(new_f);
    formatter.read_face_end();
  }

  // read a circular boundary of a face and return the first halfedge of it
  template <class Formatter_>
  D_halfedge* read_ccb_halfedges(Formatter_& formatter, D_face* face,
                                 std::size_t boundary_size)
  {
    CGAL_assertion(boundary_size > 0);
    formatter.read_ccb_halfedges_begin();
 
    // find the first halfedge, and set its face
    std::size_t first_index = 0;
    formatter.read_halfedge_index(first_index);
    D_halfedge* first_halfedge = m_halfedges[first_index];
    first_halfedge->set_face(face);

    std::size_t index, prev_index = first_index;
    for (unsigned int j = 1; j < boundary_size; ++j)
    {
      formatter.read_halfedge_index(index);
      D_halfedge* nh = m_halfedges[index];

      // connect the previous halfedge to this edge
      D_halfedge* prev_nh = m_halfedges[prev_index];
      prev_nh->set_next(nh);

      nh->set_face(face);
      prev_index = index;
    }

    // making the last halfedge point to the first one (cyclic order).
    D_halfedge* prev_halfedge = m_halfedges[prev_index];
    prev_halfedge->set_next(first_halfedge);

    formatter.read_ccb_halfedges_end();

    // return the first halfedge
    return first_halfedge;
  }
   
protected:
  Arrangement_2&            m_arr;
  std::vector<D_vertex*>    m_vertices;
  std::vector<D_halfedge*>  m_halfedges;
};

CGAL_END_NAMESPACE

#endif // CGAL_IO_ARRANGEMENT_2_READER_H 
