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
#define CGAL_IO_PM_FILE_WRITER_H 1

#ifndef CGAL_IO_BINARY_FILE_IO_H
#include <CGAL/IO/binary_file_io.h>
#endif // CGAL_IO_BINARY_FILE_IO_H

#ifndef CGAL_INVERSE_INDEX_H
#include <CGAL/Inverse_index.h>
#endif

#ifndef CGAL_IO_FILE_HEADER_H
#include <CGAL/IO/File_header.h>
#endif // CGAL_IO_FILE_HEADER_H

#ifndef CGAL_PROTECT_IOSTREAM
#include <iostream>
#define CGAL_PROTECT_IOSTREAM
#endif

#ifndef CGAL_PROTECT_CSTDDEF
#include <cstddef>
#define CGAL_PROTECT_CSTDDEF
#endif

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
 
  typedef typename PM::Holes_const_iterator           Holes_const_iterator;
  typedef typename PM::Ccb_halfedge_const_circulator  Ccb_halfedge_const_circulator;
  
  typedef Inverse_index<Halfedge_const_iterator>       H_index;
  typedef Inverse_index<Vertex_const_iterator>         V_index;

  Pm_file_writer(std::ostream& o, const PM& pm, bool verbose = false) : 
    File_header(pm.number_of_vertices(), pm.number_of_halfedges(), pm.number_of_faces(), verbose), 
    m_out(&o), pm_(pm), 
    v_index(pm.vertices_begin(), pm.vertices_end()),
    h_index(pm.halfedges_begin(), pm.halfedges_end()){}

  Pm_file_writer(std::ostream& o, const PM& pm, const File_header& h) : 
    File_header(h),  m_out(&o), pm_(pm), 
    v_index(pm.vertices_begin(), pm.vertices_end()),
    h_index(pm.halfedges_begin(), pm.halfedges_end()) { 
    
    h.set_number_of_vertices(pm.number_of_vertices());
    h.set_number_of_halfedges(pm.number_of_vertices());
    h.set_number_of_faces(pm.number_of_vertices());
  }
  
  std::ostream&          out()          { return *m_out;   }

  void write_title(const char *feature_name) { 
    if (comments()){
      out() <<"# ------------------------------------- "<< feature_name << std::endl;
      out() <<"# --------------------------------------------------------"<<std::endl;
    }
  }

  void write_comment(const char *str){
    if (comments())
      out() <<"# "<< str << std::endl;
  }

  void write_comment(const char *str, int i){
    if (comments()){
      out() <<"# "<< i <<" "<<str << std::endl;
      out() << "# ------------------------------------------"<< std::endl;
    }
  }

  void write_value(unsigned int val, char delimiter = '\n') {
    out() << val <<delimiter;
  }
  
  /*void  write_header(const char   *title_name) {
    //set_vertices(vertices);
    //set_halfedges(halfedges);
    //set_faces(faces);
    
    // Print header.
    write_title(title_name);
    } */
  
  void write_pm_vhf_sizes(){
    out() << number_of_vertices() << ' ' << number_of_halfedges() << ' ' << number_of_faces();
    out() << std::endl;
  }
  
  /*void write_footer() {
    if (comments())
      write_title("End of Planar map");
    //      out() << "#------------------- End of Planar map #";
    //out() << std::endl;
    }*/

  void write_vertices_header() {
    if (no_comments())
      out() << std::endl;
    else {
      out() << "# " << number_of_vertices() << " vertices "<< std::endl;
      out() << "# ------------------------------------------"<< std::endl;
    }
  }

  void write_vertex(Vertex_handle v) {
    out() << v->point() << std::endl;
  }

  void write_vertex(Vertex_const_handle v) {
    out() << v->point() << std::endl;
  }
  
  void write_halfedges_header() {
    if ( no_comments())
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
      write_index( v_index[ Vertex_const_iterator(h->target())] );
    
    out() << h->curve() << std::endl;
  }

  void write_halfedge(Halfedge_const_handle h) {
    if (verbose()) 
      write_vertex(h->target());
    else
      write_index( v_index[ Vertex_const_iterator(h->target())] );
  
    out() << h->curve() << std::endl;
    }*/

  void write_halfedge(Halfedge_handle h) {
    if (verbose()) {
      out() << h->curve();
      out() << "  towards  ";
      write_vertex(h->target());
    }
    else {
      write_index( v_index[ Vertex_const_iterator(h->target())] );
      out() << h->curve() << std::endl;
    }
  }

  void write_halfedge(Halfedge_const_handle h) {
    if (verbose()) {
      out() << h->curve();
      out() << "  towards  ";
      write_vertex(h->target());
    }
    else {
      write_index( v_index[ Vertex_const_iterator(h->target())] );
      out() << h->curve() << std::endl;
    }
  }

  void write_faces_header() {
    if ( no_comments())
      out() << std::endl;
    
    else {
      out() << "# " << number_of_faces() << " faces "<< std::endl;
      out() << "# ------------------------------------------"<< std::endl;
    }
  }
  
  void write_face(Face_handle f) {
    int ck;
    Holes_const_iterator iccbit;

    write_comment("writing face");
    out() << "# ------------------------------------------"<< std::endl;

    if (f->is_unbounded()){
      write_comment("UNBOUNDED");
    }
    else {
      write_comment("outer ccb");
      
      Ccb_halfedge_const_circulator first = f->outer_ccb(), iter = first;
      
      //std::size_t n = circulator_size(first);
      std::size_t n = 0;
      do {
        n++;
        iter++;
      } while (iter != first);
      write_comment("number halfedges on outer boundary");
      write_value(n);
      
      do{
        if (verbose())
          write_halfedge(iter);
        else
          write_index( h_index[Halfedge_const_iterator(iter)] );
        ++iter;
      } while (iter != first);
      
      out() << std::endl;
    }
    
    for (ck=0, iccbit = (Holes_const_iterator)f->holes_begin(); iccbit != (Holes_const_iterator)f->holes_end(); ++iccbit, ck++);
    write_comment("number of holes");
    write_value(ck);
    
    for (iccbit = (Holes_const_iterator)f->holes_begin(); iccbit != (Holes_const_iterator)f->holes_end(); ++iccbit){
      write_comment("inner ccb");
      
      Ccb_halfedge_const_circulator ccb_circ = (*iccbit);
      Ccb_halfedge_const_circulator iter = ccb_circ;
      
      //std::size_t n = circulator_size(iter);
      std::size_t n = 0;
      do {
        n++;
        iter++;
      } while (iter != ccb_circ);
      
      write_comment("number halfedges on inner boundary");
      write_value(n);
      
      do {
        if (verbose())
          write_halfedge(iter);
        else
          write_index( h_index[ Halfedge_const_iterator(iter)] );
        ++iter;
      } while (iter != ccb_circ);
      
      out() <<std::endl;
    }

    write_comment("finish writing face");
    out() << "# ------------------------------------------"<< std::endl;
  }
  
  void write_face(Face_const_handle f) {
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
          write_halfedge(iter);
        else
          write_index( h_index[Halfedge_const_iterator(iter)] );
        ++iter;
      } while (iter != first);
      
      out() << std::endl;
    }
    
    for (ck=0, iccbit = (Holes_const_iterator)f->holes_begin(); iccbit != (Holes_const_iterator)f->holes_end(); ++iccbit, ck++);
    write_comment("number of holes");
    write_value(ck);
    
    for (iccbit = (Holes_const_iterator)f->holes_begin(); iccbit != (Holes_const_iterator)f->holes_end(); ++iccbit){
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
          write_halfedge(iter);
        else
          write_index( h_index[ Halfedge_const_iterator(iter)] );
        ++iter;
      } while (iter != ccb_circ);
      
      out() <<std::endl;
    }

    write_comment("finish writing face");
    out() << "# ------------------------------------------"<< std::endl;
  }

  void write_vertices(Vertex_iterator Vertices_begin, Vertex_iterator Vertices_end) {
    for (Vertex_iterator v_iter = Vertices_begin; v_iter != Vertices_end; v_iter++)
      write_vertex(v_iter);
  }
  
  void write_vertices(Vertex_const_iterator Vertices_begin, Vertex_const_iterator Vertices_end) {
    for (Vertex_const_iterator v_iter = Vertices_begin; v_iter != Vertices_end; v_iter++)
      write_vertex(v_iter);
  }
   
  void write_halfedges(Halfedge_iterator Halfedges_begin, Halfedge_iterator Halfedges_end) {
    for (Halfedge_iterator h_iter = Halfedges_begin; h_iter != Halfedges_end; h_iter++)
      write_halfedge(h_iter);
  }

  void write_halfedges(Halfedge_const_iterator Halfedges_begin, Halfedge_const_iterator Halfedges_end) {
    for (Halfedge_const_iterator h_iter = Halfedges_begin; h_iter != Halfedges_end; h_iter++)
      write_halfedge(h_iter);
  }
  
  void write_faces(Face_iterator Faces_begin, Face_iterator Faces_end) {
    for (Face_iterator f_iter = Faces_begin; f_iter != Faces_end; f_iter++)
      write_face(f_iter);
  }

  void write_faces(Face_const_iterator Faces_begin, Face_const_iterator Faces_end) {
    for (Face_iterator f_iter = Faces_begin; f_iter != Faces_end; f_iter++)
      write_face(f_iter);
  }
  
  //void skip_line(){
  //  out() << std::endl;
  // }

protected:
  void write_index( std::size_t index) {
    out() << index <<' ';
  }
  
  //void write_string(const char *str) {
  //  out() << str << std::endl;
  // }

  std::ostream* m_out;
  const PM&     pm_;
  V_index       v_index;
  H_index       h_index;
};

CGAL_END_NAMESPACE
#endif // CGAL_IO_FILE_WRITER_PM_H //
// EOF //
