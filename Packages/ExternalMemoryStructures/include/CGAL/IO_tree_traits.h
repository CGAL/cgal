// ======================================================================
//
// Copyright (c) 1998 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 2000, August 11
//
// file          : include/CGAL/ExternalMemoryStructures/IO_tree_traits.h
// package       : ExternalMemoryStructures (0.631)
// maintainer    : Philipp Kramer <kramer@inf.ethz.ch>
// chapter       : $CGAL_Chapter: Basic / External Data Structures $
// source        : 
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Gabriele Neyer<neyer@inf.ethz.ch>
//
// coordinator   : ETH Zurich (Peter Widmayer <widmayer@inf.ethz.ch>)
//
// Implementation of the IO_tree_traits.
// ======================================================================
#ifndef __IO_tree_traits_H__
#define __IO_tree_traits_H__

#include <CGAL/basic.h>
#include <cstdio>
#include <fstream>
#include <unistd.h>

CGAL_BEGIN_NAMESPACE

//template <class InOut, int Mincap,  int Maxcap>
template <class InOut>
class IO_tree_traits {
public:
  typedef ptrdiff_t difference_type;
  InOut *IO;
  bool get(long num, char **x){
    return IO->get(num,x);
  }
  void open( int pagesize, char* name, 
	     int m = (std::fstream::binary|std::fstream::in|std::fstream::out)){
    IO=new InOut(pagesize,name,m);
  }
  void close(){
    IO->close();
  }
  difference_type get_pos() {
    return  IO->get_pos();
  }
  bool insert(difference_type pos, char **x) {
    return  IO->insert(pos, x);
  }

  difference_type number_of_elements(){
    return IO->number_of_elements();
  }
  bool update(difference_type pos, char** x){
    return IO->update(pos, x);
  }
  bool erase(difference_type pos){
    return IO->erase(pos);
  }
  bool deleted(difference_type pos){
    return IO->deleted(pos);
  }
};


//template <class InOut, int Mincap,  int Maxcap>
template <class InOut>
class IO_internal_traits {
public:
  typedef ptrdiff_t difference_type;
  InOut *IO;
  bool get(long num, char **x){
    return IO->get(num,x);
  }
  void open(int pagesize, char* name,
            int m = (std::fstream::in | std::fstream::out|std::fstream::binary))
  {
    IO=new InOut(pagesize);
  }
  void open( int pagesize){
    IO=new InOut(pagesize);
  }
  void close(){
    IO->close();
  }
  bool insert(difference_type pos, char **x) {
    return  IO->insert(pos, x);
  }

  difference_type get_pos() {
    return  IO->get_pos();
  }
  difference_type number_of_elements(){
    return IO->number_of_elements();
  }
  bool update(difference_type pos, char** x){
    return IO->update(pos, x);
  }
  bool erase(difference_type pos){
    return IO->erase(pos);
  }
  bool deleted(difference_type pos){
    return IO->deleted(pos);
  }
};


CGAL_END_NAMESPACE
#endif

