// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-2.3-I-44 $
// release_date  : $CGAL_Date: 2001/03/09 $
//
// file          : include/CGAL/IO/Pm_file_writer.h
// package       : pm (5.45)
// maintainer    : Eyal Flato <flato@math.tau.ac.il>
// source        : 
// revision      : 
// revision_date : 
// author(s)     : Eti Ezra <estere@post.tau.ac.il>
//
//
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// Chapter       : 
// ======================================================================

#ifndef CGAL_IO_PM_FILE_WRITER_H
#define CGAL_IO_PM_FILE_WRITER_H

#include <CGAL/IO/binary_file_io.h>
#include <CGAL/Inverse_index.h>
#include <CGAL/IO/File_header.h>
#include <iostream>
#include <cstddef>

CGAL_BEGIN_NAMESPACE

template <class PM>
class Pm_file_writer : public File_header{
public:
  typedef typename PM::Vertex_iterator           Vertex_iterator;
  typedef typename PM::Halfedge_iterator         Halfedge_iterator;
  typedef typename PM::Face_iterator             Face_iterator;
  typedef typename PM::Vertex_const_iterator     Vertex_const_iterator;
  typedef typename PM::Halfedge_const_iterator   Halfedge_const_iterator;
  typedef typename PM::Face_const_iterator       Face_const_iterator;

  typedef typename PM::Vertex_handle             Vertex_handle;
  typedef typename PM::Halfedge_handle           Halfedge_handle;
  typedef typename PM::Face_handle               Face_handle;
  typedef typename PM::Vertex_const_handle       Vertex_const_handle;
  typedef typename PM::Halfedge_const_handle     Halfedge_const_handle;  
  typedef typename PM::Face_const_handle         Face_const_handle;
 
  typedef typename PM::Holes_const_iterator      Holes_const_iterator;
  typedef typename PM::Ccb_halfedge_const_circulator
                                                 Ccb_halfedge_const_circulator;
  
  typedef Inverse_index<Halfedge_const_iterator> H_index;
  typedef Inverse_index<Vertex_const_iterator>   V_index;

  Pm_file_writer(std::ostream & o, const PM & pm, bool verbose = false) : 
    File_header(verbose), m_out(&o), m_pm(pm), m_h_index(0), m_v_index(0)
    { }

  Pm_file_writer(std::ostream & o, const PM & pm, const File_header & h) : 
    File_header(h),  m_out(&o), m_pm(pm), m_h_index(0), m_v_index(0)
    { }
  
  std::ostream & out() const { return *m_out; }

  void write_title(const char * feature_name) const
  { 
    if (comments()){
      out() <<"# ------------------------------------- "<< feature_name
            << std::endl;
      out() <<"# --------------------------------------------------------"
            << std::endl;
    }
  }

  void write_comment(const char * str) const
  {
    if (comments())
      out() <<"# "<< str << std::endl;
  }

  void write_comment(const char * str, int i) const
  {
    if (comments()){
      out() <<"# "<< i <<" "<<str << std::endl;
      out() << "# ------------------------------------------"<< std::endl;
    }
  }

  void write_value(unsigned int val, char delimiter = '\n') const
  {
    out() << val << delimiter;
  }
  
  /*void  write_header(const char   *title_name) {
    //set_vertices(vertices);
    //set_halfedges(halfedges);
    //set_faces(faces);
    
    // Print header.
    write_title(title_name);
    } */
  
  void write_pm_vhf_sizes(int num_vertices, int num_halfedges,
                          int num_faces) const
  {
    out() << num_vertices << ' ' << num_halfedges << ' ' << num_faces
          << std::endl;
  }
  
  /*void write_footer() {
    if (comments())
      write_title("End of Planar map");
    //      out() << "#------------------- End of Planar map #";
    //out() << std::endl;
    }*/

  void write_vertices_header() const
  {
    if (no_comments())
      out() << std::endl;
    else {
      out() << "# " << number_of_vertices() << " vertices " << std::endl;
      out() << "# ------------------------------------------" << std::endl;
    }
  }

protected:
  virtual void write_vertex(Vertex_const_handle v) const
  {
    out() << v->point() << std::endl;
  }

public:
  void write_halfedges_header() const
  {
    if (no_comments())
      out() << std::endl;
    
    else {
      out() << "# " << number_of_halfedges() << " halfedges "<< std::endl;
      out() << "# ------------------------------------------"<< std::endl;
    }
  }

  /*  void write_halfedge(Halfedge_handle h) {
    if (verbose()) 
      write_vertex(h->target());
    else
      write_index(v_index[ Vertex_const_iterator(h->target())]);
    
    out() << h->curve() << std::endl;
  }

  void write_halfedge(Halfedge_const_handle h) {
    if (verbose()) 
      write_vertex(h->target());
    else
      write_index(v_index[ Vertex_const_iterator(h->target())]);
  
    out() << h->curve() << std::endl;
    }*/

  void write_halfedge_verbose(Halfedge_const_handle h) const
  {
      out() << h->curve();
      out() << "  towards  ";
      write_vertex(h->target());
  }

protected:
  virtual void write_halfedge(Halfedge_const_handle h) const
  {
    if (verbose()) write_halfedge_verbose(h);
    else {
      write_index(get_v_index(Vertex_const_iterator(h->target())));
      out() << h->curve() << std::endl;
    }
  }

public:
  void write_faces_header() const
  {
    if (no_comments())
      out() << std::endl;
    
    else {
      out() << "# " << number_of_faces() << " faces "<< std::endl;
      out() << "# ------------------------------------------"<< std::endl;
    }
  }

protected:
  virtual void write_face(Face_const_handle f)
  {
    int ck;
    Holes_const_iterator iccbit;

    write_comment("writing face");
    out() << "# ------------------------------------------"<< std::endl;

    if (f->is_unbounded()){
      write_comment("UNBOUNDED");
      write_comment("number halfedges on outer boundary");
      write_value(0);
    }
    else {
      write_comment("outer ccb");
      
      Ccb_halfedge_const_circulator first = f->outer_ccb(), iter = first;

      std::size_t n = 0;
      do {
        n++;
        iter++;
      } while (iter != first);

      write_comment("number halfedges on outer boundary");
      write_value(n);
      
      do{
        if (verbose())
          write_halfedge_verbose(iter);
        else
          write_index(get_h_index(Halfedge_const_iterator(iter)));
        ++iter;
      } while (iter != first);
      
      out() << std::endl;
    }
    
    for (ck=0, iccbit = (Holes_const_iterator)f->holes_begin();
         iccbit != (Holes_const_iterator)f->holes_end(); ++iccbit, ck++);
    write_comment("number of holes");
    write_value(ck);
    
    for (iccbit = (Holes_const_iterator)f->holes_begin();
         iccbit != (Holes_const_iterator)f->holes_end(); ++iccbit)
    {
      write_comment("inner ccb");
      
      Ccb_halfedge_const_circulator ccb_circ = (*iccbit);
      Ccb_halfedge_const_circulator iter = ccb_circ;
      
      std::size_t n = 0;
      do {
        n++;
        iter++;
      } while (iter != ccb_circ);
      
      write_comment("number halfedges on inner boundary");
      write_value(n);
      
      do {
        if (verbose())
          write_halfedge_verbose(iter);
        else
          write_index(get_h_index(Halfedge_const_iterator(iter)));
        ++iter;
      } while (iter != ccb_circ);
      
      out() <<std::endl;
    }

    write_comment("finish writing face");
    out() << "# ------------------------------------------"<< std::endl;
  }

public:
  void write_vertices(Vertex_const_iterator Vertices_begin,
                      Vertex_const_iterator Vertices_end) const
  {
    Vertex_const_iterator v_iter;
    for (v_iter = Vertices_begin; v_iter != Vertices_end; v_iter++)
      write_vertex(v_iter);
  }

  //
  // The handling of the index arrays v_index and h_index is bad.
  // One solution is to replace write_halfedge() and write_face() with
  // corresponding classes with the index arrays as data members and
  // functors to perform the respective operations. Their constructor will
  // take the pm as a parameter to compure the index arrays.
  // 
  // void write_halfedges(Halfedge_const_iterator Halfedges_begin,
  //                      Halfedge_const_iterator Halfedges_end)
  // {
  //   Write_halfedge write_halfedge(m_pm);
  //   Halfedge_const_iterator h_iter;
  //   for (h_iter = Halfedges_begin; h_iter != Halfedges_end; h_iter++)
  //     write_halfedge(h_iter);
  // }
  //
  // The class Pm_file_writer could be templated with Write_vertex,
  // Write_halfedge, and Write_face.
  //
  // Also, I don't see any reason why write_pm should be a global function
  // and not a member of Pm_file_writer
  //
  // Efi

  void write_halfedges(Halfedge_const_iterator Halfedges_begin,
                       Halfedge_const_iterator Halfedges_end)
  {
    V_index v_index(m_pm.vertices_begin(), m_pm.vertices_end());
    m_v_index = &v_index;
    Halfedge_const_iterator h_iter;
    for (h_iter = Halfedges_begin; h_iter != Halfedges_end; h_iter++)
      write_halfedge(h_iter);
    m_v_index = 0;
  }

  void write_faces(Face_const_iterator Faces_begin,
                   Face_const_iterator Faces_end)
  {
    H_index h_index(m_pm.halfedges_begin(), m_pm.halfedges_end());
    m_h_index = &h_index;
    Face_const_iterator f_iter;
    for (f_iter = Faces_begin; f_iter != Faces_end; ++f_iter)
      write_face(f_iter);
    m_h_index = 0;
  }
  
  //void skip_line(){
  //  out() << std::endl;
  // }

protected:
  void write_index(std::size_t index) const
  {
    out() << index << ' ';
  }
  
  std::size_t get_v_index(Vertex_const_iterator vic) const
  {
    CGAL_assertion(m_v_index);
    return (*m_v_index)[vic];
  }

  std::size_t get_h_index(Halfedge_const_iterator hic) const
  {
    CGAL_assertion(m_h_index);
    return (*m_h_index)[hic];
  }

  //void write_string(const char *str) {
  //  out() << str << std::endl;
  // }

  // Data Members
  // ------------

  std::ostream * m_out;
  const PM & m_pm;
  H_index * m_h_index;
  V_index * m_v_index;
};

CGAL_END_NAMESPACE
#endif // CGAL_IO_FILE_WRITER_PM_H //
