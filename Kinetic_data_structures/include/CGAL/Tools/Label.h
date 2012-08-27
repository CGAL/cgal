// Copyright (c) 2005  Stanford University (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_LABEL_H_
#define CGAL_LABEL_H_

#include <CGAL/basic.h>
#include <iostream>
#include <sstream>

namespace CGAL {

//! A type which provides opaque, typed identifiers.
/*!
  Basically this is just an int with a type associated with it
  and one value (-1) picked out as the null value. It implements
  comparisons and things.
*/
template <class A_Type>
class Label
{
protected:
  typedef Label<A_Type> This;
  int id_;

public:
  //! Construct it from an int
  explicit Label(int i):id_(i){}
  //! Construct a default (null) label
  Label():id_(-1) {
    if(0) print();                        // make sure it is compiled
  }
  //! Make the next label
  
  CGAL_GETNR(This, successor, return This(id_+1));

  CGAL_IS(valid, return id_ >=0);

  CGAL_GETNR(unsigned int, index, CGAL_assertion(is_valid()); return id_);

  CGAL_COMPARISONS1(id_);

  std::ostream& write(std::ostream &out) const
  {
    if (id_ == -1) out << "(N)";
    else out << "(" << id_ << ")";
    return out;
  }
  void print() const
  {
    write(std::cout);
  }
  //! Convert to a string
  std::string string() const
  {
    std::ostringstream os;
    write(os);
    return os.str();
  }
};

CGAL_OUTPUT1(Label)

} //namespace CGAL

#endif
