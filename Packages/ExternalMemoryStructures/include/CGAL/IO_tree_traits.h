// Copyright (c) 1998  ETH Zurich (Switzerland).
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
// Author(s)     : Gabriele Neyer<neyer@inf.ethz.ch>
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

