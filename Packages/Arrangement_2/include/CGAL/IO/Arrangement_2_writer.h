#ifndef CGAL_IO_ARRANGEMENT_2_WRITER_H
#define CGAL_IO_ARRANGEMENT_2_WRITER_H

#include <CGAL/IO/Arrangement_2_ascii_formatter.h>
#include <CGAL/Inverse_index.h>
#include <CGAL/circulator.h>
#include <algorithm>

CGAL_BEGIN_NAMESPACE

template <class Arrangement_2>
class Arrangement_2_writer
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
  
  typedef Inverse_index<Halfedge_const_iterator>          Halfedges_index;
  typedef Inverse_index<Vertex_const_iterator>            Vertices_index;

  Arrangement_2_writer(const Arrangement_2& arr) :
    m_arr(arr)
  {
    m_h_index = new Halfedges_index(arr.halfedges_begin(),
                                    arr.halfedges_end());
    m_v_index = new Vertices_index(arr.vertices_begin(),
                                   arr.vertices_end());
  }

  virtual ~Arrangement_2_writer()
  {
    if (m_h_index) delete m_h_index;
    if (m_v_index) delete m_v_index;
  }

  template <class Formatter_>
  void operator()(Formatter_& formatter) const
  {
    formatter.write_arr_begin(m_arr);
    formatter.write_size("number_of_vertices", m_arr.number_of_vertices());
    formatter.write_size("number_of_halfedges", m_arr.number_of_halfedges());
    formatter.write_size("number_of_faces", m_arr.number_of_faces());

    // write the vertices
    formatter.write_vertices_begin();
    Vertex_const_iterator v_iter = m_arr.vertices_begin();
    for (; v_iter != m_arr.vertices_end(); ++v_iter)
      write_vertex(formatter, v_iter);
    formatter.write_vertices_end();

    // write the edges & halfedges
    formatter.write_halfedges_begin();
    Halfedge_const_iterator h_iter = m_arr.halfedges_begin();
    for (; h_iter != m_arr.halfedges_end(); ++h_iter, ++h_iter)
      write_edge(formatter, h_iter);
    formatter.write_halfedges_end();

    // write faces
    formatter.write_faces_begin();
    Face_const_iterator f_iter = m_arr.faces_begin();
    for (; f_iter != m_arr.faces_end(); ++f_iter)
      write_face(formatter, f_iter);
    formatter.write_faces_end();

    formatter.write_arr_end();
  }

protected:

  template <class Formatter_>
  void write_vertex(Formatter_& formatter, Vertex_const_iterator v) const
  {
    formatter.write_vertex_begin(v, get_index(v));
    formatter.write_vertex_point(v);
    // here users that expand the vertex can write additional data
    formatter.write_vertex(v); 
    formatter.write_vertex_end(v);
  }

  template <class Formatter_>
  void write_edge(Formatter_& formatter, Halfedge_const_iterator hi) const
  {
    Halfedge_const_handle h = hi;
    formatter.write_edge_begin(h/*, get_index(hi)*/);
    formatter.write_halfedge_endpoint_index(get_index(h->source()), "endp1");
    formatter.write_halfedge_endpoint_index(get_index(h->target()), "endp2");
    formatter.write_edge_curve(h);
    formatter.write_halfedge(h);
    formatter.write_halfedge(h->twin());
    formatter.write_edge_end(h);
  }

  template <class Formatter_>
  void write_face(Formatter_& formatter, Face_const_iterator fi) const
  {
    Face_const_handle f = fi;
    formatter.write_face_begin(f);
    // write the outer ccb:
    // first write the number of halfedges on it, then write the halfedges
    // indices
    formatter.write_outer_ccb_begin(f);                                                  
    if (f->is_unbounded())
    {
      formatter.write_number("halfedges_on_outer_boundary", 0);
    }
    else
    {
      Ccb_halfedge_const_circulator hec = f->outer_ccb();

      std::size_t n = circulator_size(hec);
      formatter.write_number("halfedges_on_outer_boundary", n);

      write_ccb_halfedges(formatter, hec);
    }
    formatter.write_outer_ccb_end(f);

    // write the holes of the face
    formatter.write_holes_begin(f);
    std::size_t n_holes = std::distance(f->holes_begin(), f->holes_end());
    formatter.write_number("number_of_holes", n_holes);

    Holes_const_iterator hole_iter = (Holes_const_iterator)f->holes_begin();
    for (; hole_iter != (Holes_const_iterator)f->holes_end(); ++hole_iter)
    {
      formatter.write_inner_ccb_begin(f);

      Ccb_halfedge_const_circulator hec = (*hole_iter);      
      std::size_t n = circulator_size(hec);
      formatter.write_number("halfedges_on_inner_boundary", n);
      write_ccb_halfedges(formatter, hec);
      
      formatter.write_inner_ccb_end(f);
    }

    formatter.write_holes_end(f);

    // write the isolated vertices of the face
    formatter.write_isolated_vertices_begin(f);
    std::size_t n_isolated = std::distance(f->isolated_vertices_begin(),
                                           f->isolated_vertices_end());
    formatter.write_number("number_of_isolated_vertices", n_isolated);
    Isolated_vertices_const_iterator iso_vi = f->isolated_vertices_begin();
    for(; iso_vi != f->isolated_vertices_end(); ++iso_vi)
      formatter.write_vertex_index(get_index(iso_vi));
    
    formatter.write_isolated_vertices_end(f);
    
    // this is used for writing face additional data
    formatter.write_face(f);
    formatter.write_face_end(f);
  }

  template <class Formatter_>   
  void write_ccb_halfedges(Formatter_& formatter,
                           Ccb_halfedge_const_circulator hec_begin) const      
  {
    formatter.write_ccb_halfedges_begin();
    Ccb_halfedge_const_circulator hec = hec_begin;
    do {
      formatter.write_halfedge_index(get_index(hec));
      ++hec;
    } while (hec != hec_begin);
    formatter.write_ccb_halfedges_end();
  }
  
  std::size_t get_index(Vertex_const_iterator v) const
  {
    CGAL_assertion(m_v_index);
    return (*m_v_index)[v];
  }

  std::size_t get_index(Halfedge_const_iterator h) const
  {
    CGAL_assertion(m_h_index);
    return (*m_h_index)[h];
  }

  std::size_t my_circulator_size(Ccb_halfedge_const_circulator c) const
  {
    // Simply count.
    if ( c == CGAL_CIRC_NULL)
        return 0;

    std::size_t     n = 0;
    Ccb_halfedge_const_circulator   d = c;
    do {
        ++n;
        ++d;
    } while( c != d);
    return n;  
  }
  
protected:
  const Arrangement_2& m_arr;
  Halfedges_index*     m_h_index;
  Vertices_index*      m_v_index;
};

CGAL_END_NAMESPACE

#endif // CGAL_IO_ARRANGEMENT_2_WRITER_H 
