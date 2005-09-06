#ifndef CGAL_IO_ARRANGEMENT_2_ASCII_FORMATTER_H
#define CGAL_IO_ARRANGEMENT_2_ASCII_FORMATTER_H

#include <CGAL/basic.h>
#include <iostream>

CGAL_BEGIN_NAMESPACE

template <class Arrangement_2>
class Arrangement_2_ascii_formatter
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

public:  
  Arrangement_2_ascii_formatter(std::ostream& o) : m_out(&o), m_in(NULL)
  {}

  Arrangement_2_ascii_formatter(std::istream& i, std::ostream& o) : m_out(&o), m_in(&i)
  {}

  Arrangement_2_ascii_formatter(std::istream& i) : m_out(NULL), m_in(&i)
  {}

  virtual ~Arrangement_2_ascii_formatter()
  {}

  std::ostream& out() const { return *m_out; }
  std::istream& in() { return *m_in; }

  // write functions
  void write_arr_begin(const Arrangement_2& m_arr)
  {
	  CGAL_assertion(m_out);
    m_old_mode = get_mode(*m_out);
    set_ascii_mode(*m_out);
    write_comment("************************************************");
    write_comment("Begin planar map");
    write_comment("************************************************");
  }

  void write_arr_end() const
  {
    write_comment("************************************************");
    write_comment("End planar map");
    write_comment("************************************************");
    set_mode(*m_out, m_old_mode);
  }

  void write_value(unsigned int val, char delimiter = '\n') const
  {
    out() << val << delimiter;
  }

  void write_comment(const char *str) const
  {
    out() <<"# "<< str << std::endl;
  }

  void write_size(const char *title, Size size)
  { 
	  out() << size << std::endl; 
  }

  void write_number(const char *title, int num)
  {
    write_comment(title);
    write_value(num);  
  }

  void write_vertices_begin() const
  {
    write_comment("Vertices:");
    write_comment("------------------------------------------");    
  }
  void write_vertices_end() const {}

  void write_halfedges_begin() const
  {
    write_comment("Halfedges:");
    write_comment("------------------------------------------");    
  }
  void write_halfedges_end() const {}

  void write_faces_begin() const
  {
    write_comment("Faces:");
    write_comment("------------------------------------------");
  }
  void write_faces_end() const {}

  // vertex
  void write_vertex_begin(Vertex_const_handle v, std::size_t idx) const
  {}
  void write_vertex_end(Vertex_const_handle v) const
  {}
  
  void write_vertex_point(Vertex_const_handle v) const
  {
    out() << v->point() << std::endl;
  }

  void write_vertex(Vertex_const_handle v) const {}

  // edge & halfedge
  void write_edge_begin(Halfedge_const_handle h) const
  {}
  void write_edge_end(Halfedge_const_handle h) const
  {}

  void write_halfedge_endpoint_index(std::size_t idx, const char *title=0)
  {
    write_value(idx, ' ');
  }

  void write_edge_curve(Halfedge_const_handle h) const
  {
    out() << h->curve() << std::endl;
  }
  void write_halfedge(Halfedge_const_handle h) const {}

  // face
  void write_face_begin(Face_const_handle f) const
  {
    if (f->is_unbounded())
      write_comment("start UNBOUNDED face");
    else
      write_comment("start BOUNDED face");
    write_comment("------------------------------------------");      
  }
  void write_face_end(Face_const_handle f) const
  {
    write_comment("end face");
    write_comment("------------------------------------------");    
  }

  void write_outer_ccb_begin(Face_const_handle f) const
  {
    write_comment("outer ccb");
  }
  void write_outer_ccb_end(Face_const_handle f) const
  {}

  void write_holes_begin(Face_const_handle f) const
  {}
  void write_holes_end(Face_const_handle f) const
  {}

  void write_inner_ccb_begin(Face_const_handle f) const
  {
    write_comment("inner ccb");
  }
  void write_inner_ccb_end(Face_const_handle f) const
  {}

  void write_face(Face_const_handle f) const {}

  void write_ccb_halfedges_begin() const {}
  void write_ccb_halfedges_end() const
  {
    out() << std::endl;    
  }

  void write_halfedge_index(std::size_t idx) const
  {
    write_value(idx, ' ');
  }

  void write_isolated_vertices_begin(Face_const_handle f) const
  {
    write_comment("isolated vertices");
  }
  void write_isolated_vertices_end(Face_const_handle f) const
  {
    out() << std::endl;
  }

  void write_vertex_index(std::size_t idx) const
  {
    write_value(idx, ' ');
  }  

  // read functions
  void read_arr_begin() 
  {
   	CGAL_assertion(m_in);
    m_old_mode = get_mode(*m_in);
    set_ascii_mode(*m_in);
	  skip_comments(in());
  }
  void read_arr_end() 
  {
    skip_comments(in());
    set_mode(*m_in, m_old_mode);
  }

  Size read_size(const char *title=0)
  {
    Size size;
    in() >> size;
	  skip_until_EOL(in());
	  return size;
  }

  int read_number(const char *title=0)
  {
    skip_comments(in());
  	int tmp;
  	in() >> tmp;
  	return tmp;
  }

  void read_vertices_begin() { skip_comments(in()); }
  void read_vertices_end() {}

  void read_halfedges_begin() { skip_comments(in()); }
  void read_halfedges_end() {}

  void read_faces_begin() { skip_comments(in()); }
  void read_faces_end() {}

  // vertex
  void read_vertex_begin() {}
  void read_vertex_end() {}
  void read_vertex_point(Point_2& p) 
  {
    in() >> p;
	skip_until_EOL(in());
  }
  void read_vertex(D_vertex* v) {}
  
  // edge & halfedge
  void read_edge_begin() {}
  void read_edge_end() {}
  std::size_t read_halfedge_endpoint_index(const char *title=0) 
  {
    std::size_t result;  
	  in() >> result;
	  return result;
  }
  void read_halfedge_curve(X_monotone_curve_2& cv) 
  {
    in() >> cv;
	skip_until_EOL(in());
  }

  void read_halfedge(D_halfedge* h) {}
 
  // face
  void read_face_begin() { skip_comments(in()); }
  void read_face_end() { skip_comments(in()); }

  void read_outer_ccb_begin() { skip_comments(in()); }
  void read_outer_ccb_end() {}

  void read_halfedge_index(std::size_t& index) { in() >> index; }

  void read_holes_begin() {}
  void read_holes_end() {}

  void read_inner_ccb_begin() { skip_comments(in()); }
  void read_inner_ccb_end() {}

  void read_ccb_halfedges_begin() {}
  void read_ccb_halfedges_end() 
  {
    skip_until_EOL(in());
  }

  void read_isolated_vertices_begin() { skip_comments(in()); }
  void read_isolated_vertices_end() 
  {
    skip_until_EOL(in());
  }

  void read_vertex_index(std::size_t& index) { in() >> index; }

  void read_face(D_face* f) {}


  // istream modifier skips chars until end of line.
  std::istream& skip_until_EOL( std::istream& in) 
  {
    char c;
    while ( in.get(c) && c != '\n')
      ;
    return in;
  }
  
  // istream modifier that checks for OFF comments and removes them.
  std::istream& skip_comments( std::istream& in) 
  {
    char c;
    while( (in >> c) && c == '#')
      skip_until_EOL(in);
    in.putback(c);
    return in;
  }

protected:
  std::ostream * m_out;
  std::istream * m_in;
  IO::Mode       m_old_mode;
};

CGAL_END_NAMESPACE

#endif // CGAL_IO_ARRANGEMENT_2_ASCII_FORMATTER_H 
