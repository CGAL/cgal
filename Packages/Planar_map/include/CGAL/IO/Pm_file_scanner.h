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
// file          : include/CGAL/IO/Pm_file_scanner.h
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

#ifndef CGAL_IO_PM_FILE_SCANNER_H
#define CGAL_IO_PM_FILE_SCANNER_H 1

#ifndef CGAL_BASIC_H
#include <CGAL/basic.h>
#endif
#ifndef CGAL_KNOWN_BIT_SIZE_INTEGERS_H
#include <CGAL/known_bit_size_integers.h>
#endif
#ifndef CGAL_PROTECT_CSTDDEF
#include <cstddef>
#define CGAL_PROTECT_CSTDDEF
#endif
#ifndef CGAL_IO_BINARY_FILE_IO_H
#include <CGAL/IO/binary_file_io.h>
#endif // CGAL_IO_BINARY_FILE_IO_H

#ifndef CGAL_IO_FILE_HEADER_H
#include <CGAL/IO/File_header.h>
#endif // CGAL_IO_FILE_HEADER_H

#ifndef CGAL_PROTECT_IOSTREAM
#include <iostream>
#define CGAL_PROTECT_IOSTREAM
#endif

CGAL_BEGIN_NAMESPACE

template <class PM_>
class Pm_file_scanner : public File_header {
public:

  typedef PM_                                    Planar_map;

  typedef typename Planar_map::Traits            Traits;
  typedef typename Traits::Point                 Point;
  typedef typename Traits::X_curve               X_curve;

  typedef typename Planar_map::Dcel::Vertex	         D_vertex;
  typedef typename Planar_map::Dcel::Halfedge            D_halfedge;
  typedef typename Planar_map::Dcel::Face	         D_face;

  Pm_file_scanner( std::istream& in) : File_header(), m_in(in) { skip_comment(); }
  
  Pm_file_scanner( std::istream& in, const File_header& header) : File_header(header), m_in(in) { skip_comment(); }
  
  std::istream& in() { return m_in; }
 
  void scan_pm_vhf_sizes (){
    std::size_t v, h, f;

    in() >> v >> h >> f;

    set_number_of_vertices(v);
    set_number_of_halfedges(h);
    set_number_of_faces(f);
  }
 
  void scan_vertex (D_vertex* v){

    Point p;
    skip_comment();

    // providing default reading function.
    in() >> p;
    if ( !in()){
      std::cerr << "can't read face number"<<std::endl;
      in().clear( std::ios::badbit);
      return;
    }
    
    skip_to_next_vertex();
    if ( !in()){
      std::cerr << "can't skip to next vertex"<<std::endl;
      in().clear( std::ios::badbit);
      return;
    }
    
    v->set_point(p);
  }

  X_curve scan_halfedge (D_halfedge* h){
    // providing default reading function.
    X_curve cv;

    skip_comment();

    in() >> cv;
    if ( !in()){
      std::cerr << "can't read face number"<<std::endl;
      in().clear( std::ios::badbit);
      return cv;
    }
    
    skip_to_next_halfedge();
    if ( ! in()){
      std::cerr << "can't skip to next halfedge"<<std::endl;
      in().clear( std::ios::badbit);
      return cv;
    } 
    
    h->set_curve(cv);

    halfedges_vec.push_back(h);

    return cv;
  }
  
  //void  update_halfedges_vec(D_halfedge* h)
  // {
  // halfedges_vec.push_back(h);
  // }

  void scan_face(D_face* f) {
    
    std::size_t  num_of_holes, num_halfedges_on_outer_ccb, i = 0;
    
    scan_face_number(num_halfedges_on_outer_ccb, i);
    if ( !in()){
      std::cerr << "can't read face number"<<std::endl;
      in().clear( std::ios::badbit);
      return;
    }
    
    //  not an unbounded face. Scanning the outer ccb.
    if (num_halfedges_on_outer_ccb > 0){
      std::size_t  index, prev_index, first_index;
      
      for (unsigned int j = 0; j < num_halfedges_on_outer_ccb; j++) {
        
        scan_index(index);
        if ( !in()){
          std::cerr << "can't read halfedge's index on face"<<std::endl;
          in().clear( std::ios::badbit);
          return;
        }
        
        D_halfedge* nh = halfedges_vec[index];
        
        // for debugging.
        //std::cout<<"source of haledge : "<<nh->vertex()->point()<<std::endl;
        
        if (j > 0) {
          D_halfedge* prev_nh = halfedges_vec[prev_index];
          prev_nh->set_next(nh);
        }
        else {
          f->set_halfedge(nh);
          first_index = index;
        }
        
        nh->set_face(f); 
        
        prev_index = index;
      }
      
      // making the last halfedge point to the first one (cyclic order).
      D_halfedge* nh = halfedges_vec[first_index];
      D_halfedge* prev_nh = halfedges_vec[prev_index];
      prev_nh->set_next(nh);
    }
    
    scan_face_number(num_of_holes, i);
    if ( !in()){
      std::cerr << "can't read number holes in face"<<std::endl;
      in().clear( std::ios::badbit);
      return;
    }
    
    // take care the hols.
    for (unsigned int k = 0; k < num_of_holes; k++){
      std::size_t  num_halfedges_on_inner_ccb;
      
      scan_face_number(num_halfedges_on_inner_ccb, i);
      if ( !in()){
        std::cerr << "can't read number of halfedges in hole"<<std::endl;
        in().clear( std::ios::badbit);
        return;
      }
      
      std::size_t  index, prev_index, first_index;
      for (unsigned int j = 0; j < num_halfedges_on_inner_ccb; j++) {
        scan_index(index);
        if ( !in()){
          std::cerr << "can't read halfedge's index on hole"<<std::endl;
          in().clear( std::ios::badbit);
          return;
        }
        
        D_halfedge* nh = halfedges_vec[index];
        
        // for debugging.
        //std::cout<<"source of haledge : "<<nh->vertex()->point()<<std::endl;
        
        if (j > 0) {
          D_halfedge* prev_nh = halfedges_vec[prev_index];
          prev_nh->set_next(nh);
        }
        else {
          f->add_hole(nh);
          first_index = index;
        }
        
        nh->set_face(f); 
        
        prev_index = index;
      }
      
      // making the last halfedge point to the first one (cyclic order).
      D_halfedge* nh = halfedges_vec[first_index];
      D_halfedge* prev_nh = halfedges_vec[prev_index];
      prev_nh->set_next(nh);
    }
    
    skip_to_next_face();
    if ( !in()){
      std::cerr << "can't skip to next face"<<std::endl;
      in().clear( std::ios::badbit);
      return;
    }
  }

  void scan_index(std::size_t& index){ 
    
    skip_comment();
    skip_comment();

    in() >> index;
    
    if( ! in()) {
      return;
    } 

    index -= index_offset();
  }
  
protected:
  std::istream&               m_in;
  std::vector<D_halfedge* >   halfedges_vec;
  
  void scan_face_vertex_index(std::size_t& index , std::size_t current_face){ 
    
    in() >> index;
  
    if( ! in()) {
      return;
    }
    
    index -= index_offset();
    
    /*  if( index < 0 || index >= size_of_halfedges()) {
        in().clear( std::ios::badbit);
        
        std::cerr << " " << std::endl;
        std::cerr << "File_scanner_pm::" << std::endl;
        std::cerr << "scan_facet_vertex_index(): input error: "
        "facet " << current_face << ": vertex index "
        << index + index_offset() << ": is out of range."
        << std::endl;
        
        //set_off_header( false);
        return;
        }*/
    
    //std::cout<< index <<std::endl; // for debuging.
  }

  void scan_face_number( std::size_t& size,  std::size_t current_face) {
    
    CGAL_assertion( current_face < number_of_faces());
    
    skip_comment();
    
    in() >> size;

    //std::cout<< size <<std::endl; // for debuging.
  }

  void skip_to_next_vertex(){
    //in() >> skip_until_EOL
    skip_until_EOL(in());
  }
  
  void skip_to_next_halfedge(){
    skip_until_EOL(in());
  }
  
  void skip_to_next_face(){
    skip_until_EOL(in());
  }

  void skip_comment() { 

    char c;
    while( (in() >> c) && c == '#'){
      in().putback(c);
      skip_comment_OFF(in());
      //in() >> skip_comment_OFF; 
    }
    in().putback(c);
  }
};


CGAL_END_NAMESPACE
#endif // CGAL_IO_FILE_SCANNER_PM_H //
// EOF //



