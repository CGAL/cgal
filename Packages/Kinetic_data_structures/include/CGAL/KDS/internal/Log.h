// Copyright (c) 2005  Stanford University (USA).
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
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_KDS_LOG_H_
#define CGAL_KDS_LOG_H_
#include <CGAL/KDS/basic.h>
#include <iostream>
#include <fstream>
#include <ios>

CGAL_KDS_BEGIN_INTERNAL_NAMESPACE


class Logs{
public:
  // The different types of logs supported
  /*  MAPLE is a log which should be able to be fed directly in to
    maple and preferably will produce obviously good or bad outwhen
    when evaluated.
  */
  typedef CGAL::KDS::Log_level Level;
  typedef enum Target {COUT, FILE, DEVNULL} Target;

  Level level() const {
    return level_;
  }

  void set_level(Level l){
    level_= l;
  }

  std::ostream &stream(Level l) {
    if (is_output(l)){
      if (target_== COUT){
	return std::cout;
      } else {
	return fstream_;
      }
    } else {
      return null_;
    }
  }
  bool is_output(Level l){
    if (level_== CGAL::KDS::LOG_NONE && l != level_) return false;
    else if (level_==CGAL::KDS::LOG_SOME && l== CGAL::KDS::LOG_SOME) return true;
    else return true; 
  }
  Target target() const {
    return target_;
  }
  void set_target(Target t){
    target_=t;
  }
  
  void set_filename(const char *name){
    fstream_.open(name);
  }

  static Logs& get() {
    static Logs lg;
    return lg;
  }
  
  bool output_maple() const {
    return output_maple_;
  }
  void set_output_maple(bool b) {
    output_maple_=b;
  }
  std::ofstream &maple_stream() {
    if (!maple_is_open_){
      maple_is_open_=true;
      maple_.open("maple.log");
    }
    return maple_;
  }
  
protected:
  Logs(){
    level_= LOG_NONE; 
    target_= COUT;
    null_.open("/dev/null");
    //maple_.open("maple.log");
    maple_is_open_=false;
    output_maple_=true;
  }


  Target target_;
  Level level_;
  std::ofstream fstream_;
  std::ofstream null_;
  std::ofstream maple_;
  bool maple_is_open_;
  bool output_maple_;
};



CGAL_KDS_END_INTERNAL_NAMESPACE

#endif
