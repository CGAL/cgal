// Copyright (c) 1997  Tel-Aviv University (Israel).
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
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>

#ifndef SWEEP_LINE_TRAITS_H
#define SWEEP_LINE_TRAITS_H


CGAL_BEGIN_NAMESPACE
template<class Traits>
struct Sweep_line_traits
{
  static Traits* static_traits;

  Sweep_line_traits(Traits *traits)
  {
    static_traits = traits;
  }

 
  static Traits* get_traits()
  {
    return static_traits;
  }
   
  ~Sweep_line_traits()
  {
  }

};


//TODO : move the below two line to Sweep_line_traits.C file and write instead
 //#include<CGAL/Sweep_line_2/Sweep_line_traits.C>
template<class Traits>
Traits* Sweep_line_traits<Traits>::static_traits = NULL;

CGAL_END_NAMESPACE

#endif
