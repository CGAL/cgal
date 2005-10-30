// Copyright (c) 2005  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Michal Meyerovitch <gorgymic@post.tau.ac.il>
//                 Ron Wein           <wein@post.tau.ac.il>
//                 (based on old version by Ester Ezra)
#ifndef CGAL_IO_ARRANGEMENT_2_ASCII_FORMATTER_H
#define CGAL_IO_ARRANGEMENT_2_ASCII_FORMATTER_H

/*! \file
 * The header file for the Arrangement_2_ascii_formatter<Arrangement> class.
 */

#include <CGAL/basic.h>
#include <iostream>

CGAL_BEGIN_NAMESPACE

/*! \class
 * A class defining a textual (ASCII) input/output format for arrangements
 * and supports reading and writing an arrangement from or to input/output
 * streams.
 */
template <class Arrangement_>
class Arrangement_2_ascii_formatter
{

public:

  typedef Arrangement_                                    Arrangement_2;
  typedef typename Arrangement_2::Size                    Size;
  typedef typename Arrangement_2::Dcel                    Dcel;
  typedef typename Arrangement_2::Traits_2                Traits_2;
  typedef typename Traits_2::X_monotone_curve_2           X_monotone_curve_2;
  typedef typename Traits_2::Point_2                      Point_2;

  typedef typename Arrangement_2::Vertex_iterator         Vertex_iterator;
  typedef typename Arrangement_2::Halfedge_iterator       Halfedge_iterator;
  typedef typename Arrangement_2::Face_iterator           Face_iterator;

  typedef typename Arrangement_2::Vertex_handle           Vertex_handle;
  typedef typename Arrangement_2::Halfedge_handle         Halfedge_handle;
  typedef typename Arrangement_2::Face_handle             Face_handle;

  typedef typename Arrangement_2::Vertex_const_iterator
                                                      Vertex_const_iterator;
  typedef typename Arrangement_2::Halfedge_const_iterator
                                                      Halfedge_const_iterator;
  typedef typename Arrangement_2::Face_const_iterator 
                                                      Face_const_iterator;

  typedef typename Arrangement_2::Vertex_const_handle 
                                                      Vertex_const_handle;
  typedef typename Arrangement_2::Halfedge_const_handle
                                                      Halfedge_const_handle;
  typedef typename Arrangement_2::Face_const_handle   
                                                      Face_const_handle;
 
  typedef typename Arrangement_2::Holes_const_iterator
                                             Holes_const_iterator;
  typedef typename Arrangement_2::Ccb_halfedge_const_circulator
                                             Ccb_halfedge_const_circulator;
  typedef typename Arrangement_2::Isolated_vertices_const_iterator
                                             Isolated_vertices_const_iterator;


protected:

  typedef typename Dcel::Vertex                           DVertex;
  typedef typename Dcel::Halfedge                         DHalfedge;
  typedef typename Dcel::Face                             DFace;
  
  // Data members:
  std::ostream  *m_out;
  std::istream  *m_in;
  IO::Mode       m_old_mode;

public:  

  /*! Construct an output formatter. */
  Arrangement_2_ascii_formatter (std::ostream& os) :
    m_out (&os),
    m_in(NULL)
  {}

  /*! Construct an input formatter. */
  Arrangement_2_ascii_formatter (std::istream& is) :
    m_out(NULL),
    m_in(&is)
  {}

  /*! Destructor. */
  virtual ~Arrangement_2_ascii_formatter()
  {}

  /*! Get the output stream. */
  inline std::ostream& out ()
  {
    CGAL_assertion (m_out != NULL);
    return (*m_out);
  }

  /*! Get the input stream. */
  inline std::istream& in ()
  {
    CGAL_assertion (m_in != NULL);
    return (*m_in);
  }

  /// \name Global write functions.
  //@{

  /*! Write a begin-arrangement comment. */
  void write_arr_begin (const Arrangement_2& arr)
  {
    CGAL_assertion (m_out != NULL);
    m_old_mode = get_mode(*m_out);
    set_ascii_mode (*m_out);
    write_comment("BEGIN ARRANGEMENT");

    return;
  }

  /*! Write an end-arrangement comment. */
  void write_arr_end()
  {
    write_comment("END ARRANGEMENT");
    set_mode (*m_out, m_old_mode);

    return;
  }

  /*! Write an unsigned integer value. */
  void write_value (unsigned int val, char delimiter = '\n')
  {
    out() << val << delimiter;
    return;
  }

  /*! Write a comment line. */
  void write_comment (const char *str)
  {
    out() << "# " << str << std::endl;
    return;
  }

  /*! Write a labeled size value. */
  void write_size (const char *label, Size size)
  { 
    write_comment (label);
    write_value (size);  
    return;
  }

  /*! Write a begin-vertices comment. */
  void write_vertices_begin()
  {
    write_comment("BEGIN VERTICES");
    return;
  }

  /*! Write an end-vertices comment. */
  void write_vertices_end()
  {
    write_comment("END VERTICES");
    return;
  }

  /*! Write a begin-edges comment. */
  void write_edges_begin()
  {
    write_comment("BEGIN EDGES");
    return;
  }

  /*! Write an end-edges comment. */
  void write_edges_end()
  {
    write_comment("END EDGES");
    return;
  }

  /*! Write a begin-faces comment. */
  void write_faces_begin()
  {
    write_comment("BEGIN FACES");
    return;
  }

  /*! Write an end-faces comment. */
  void write_faces_end()
  {
    write_comment("END FACES");
    return;
  }
  //@}

  /// \name Write a vertex.
  //@{
  void write_vertex_begin (Vertex_const_handle , std::size_t )
  {}

  void write_vertex_end (Vertex_const_handle )
  {}
  
  void write_point (const Point_2& p)
  {
    out() << p << std::endl;
    return;
  }

  virtual void write_vertex_data (Vertex_const_handle v)
  {}
  //@}

  /// \name Write an edge.
  //@{
  void write_edge_begin (Halfedge_const_handle )
  {}

  void write_edge_end (Halfedge_const_handle )
  {}

  void write_vertex_index (std::size_t idx)
  {
    write_value (idx, ' ');
  }

  void write_curve (const X_monotone_curve_2& cv)
  {
    out() << cv << std::endl;
    return;
  }

  virtual void write_halfedge_data (Halfedge_const_handle he)
  {}
  //@}

  /// \name Write a face.
  //@{
  void write_face_begin (Face_const_handle f)
  {
    if (f->is_unbounded())
      write_comment("BEGIN FACE (unbounded)");
    else
      write_comment("BEGIN FACE (bounded)");
    return;
  }

  void write_face_end (Face_const_handle )
  {
    write_comment("END FACE");
  }

  void write_outer_ccb_begin (Face_const_handle )
  {
    write_comment("Outer CCB:");
  }

  void write_outer_ccb_end(Face_const_handle )
  {}

  void write_holes_begin(Face_const_handle )
  {}

  void write_holes_end(Face_const_handle )
  {}

  void write_inner_ccb_begin(Face_const_handle )
  {
    write_comment("Inner CCB:");
  }

  void write_inner_ccb_end(Face_const_handle )
  {}

  virtual void write_face_data (Face_const_handle )
  {}

  void write_ccb_halfedges_begin()
  {}
  
  void write_ccb_halfedges_end()
  {
    out() << std::endl;    
  }

  void write_halfedge_index (std::size_t idx)
  {
    write_value (idx, ' ');
  }

  void write_isolated_vertices_begin (Face_const_handle )
  {
    write_comment("Isolated vertices:");
  }

  void write_isolated_vertices_end(Face_const_handle )
  {
    out() << std::endl;
  }
  //@}

  /// \name Global read functions.
  //@{

  /*! Start reading an arrangement. */
  void read_arr_begin () 
  {
    CGAL_assertion (m_in != NULL);
    m_old_mode = get_mode(*m_in);
    set_ascii_mode(*m_in);
    _skip_comments();
  }

  /*! Read the arrangement edge. */
  void read_arr_end() 
  {
    _skip_comments();
    set_mode(*m_in, m_old_mode);
  }

  /*! Read a size value (with a label comment line before it). */
  std::size_t read_size (const char *title = NULL)
  {
    std::size_t   val;

    _skip_comments();
    in() >> val;
    _skip_until_EOL();

    return (val);
  }

  /*! Reading the arrangement vertices. */
  void read_vertices_begin()
  {
    _skip_comments();
  }

  void read_vertices_end()
  {
    _skip_comments();
  }

  /*! Reading the arrangement edges. */
  void read_edges_begin()
  {
    _skip_comments();
  }

  void read_edges_end()
  {
    _skip_comments();
  }

  /*! Reading the arrangement faces. */
  void read_faces_begin()
  {
    _skip_comments();
  }

  void read_faces_end()
  {
    _skip_comments();
  }
  //@}

  /// \name Reading a vertex.
  //@{
  void read_vertex_begin ()
  {}
  
  void read_vertex_end ()
  {}

  void read_point (Point_2& p) 
  {
    in() >> p;
    _skip_until_EOL();

    return;
  }

  virtual void read_vertex_data (Vertex_handle v)
  {}
  //@}

  /// \name Reading an edge.
  //@{
  void read_edge_begin ()
  {}
  
  void read_edge_end ()
  {}
  
  std::size_t read_vertex_index () 
  {
    std::size_t  val;

    in() >> val;
    return (val);
  }

  void read_curve (X_monotone_curve_2& cv) 
  {
    in() >> cv;
    _skip_until_EOL();

    return;
  }

  virtual void read_halfedge_data (Halfedge_handle )
  {}
 
  /// \name Reading a face.
  //@{
  void read_face_begin ()
  {
    _skip_comments();
  }
  
  void read_face_end ()
  {
    _skip_comments();
  }

  void read_outer_ccb_begin ()
  {
    _skip_comments();
  }
  
  void read_outer_ccb_end ()
  {}

  std::size_t read_halfedge_index ()
  { 
    std::size_t  val;

    in() >> val;
    return (val);
  }

  void read_holes_begin ()
  {}
  
  void read_holes_end ()
  {}

  void read_inner_ccb_begin ()
  {
    _skip_comments();
  }
  
  void read_inner_ccb_end ()
  {}

  void read_ccb_halfedges_begin()
  {}
  
  void read_ccb_halfedges_end() 
  {
    _skip_until_EOL ();
  }

  void read_isolated_vertices_begin ()
  {
    _skip_comments();
  }
  
  void read_isolated_vertices_end () 
  {
    _skip_until_EOL();
  }

  virtual void read_face_data (Face_handle )
  {}
  //@}

protected:

  /*! Skip until end of line. */
  void _skip_until_EOL () 
  {
    CGAL_assertion (m_in != NULL);

    int     c;
    while ((c = m_in->get()) != EOF && c != '\n');
    return;
  }
  
  /*! Skip comment lines. */
  void _skip_comments () 
  {
    CGAL_assertion (m_in != NULL);

    int     c;
    while ((c = m_in->get()) != EOF && c == '#')
      _skip_until_EOL();
    m_in->putback (c);

    return;
  }

};

CGAL_END_NAMESPACE

#endif // CGAL_IO_ARRANGEMENT_2_ASCII_FORMATTER_H 
