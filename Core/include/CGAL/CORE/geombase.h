/****************************************************************************
 * Core Library Version 1.7, August 2004                                     
 * Copyright (c) 1995-2004 Exact Computation Project                         
 * All rights reserved.                                                      
 *                                                                           
 * This file is part of CORE (http://cs.nyu.edu/exact/core/).                
 * You can redistribute it and/or modify it under the terms of the GNU       
 * General Public License as published by the Free Software Foundation,      
 * either version 3 of the License, or (at your option) any later version.   
 *                                                                           
 * Licensees holding a valid commercial license may use this file in         
 * accordance with the commercial license agreement provided with the        
 * software.                                                                 
 *                                                                           
 * This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE   
 * WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. 
 *                                                                           
 *                                                                           
 * $URL$                                                                     
 * $Id$                                                                      
 ***************************************************************************/
/******************************************************************
 * Core Library Version 1.7, August 2004
 * Copyright (c) 1995-2002 Exact Computation Project
 * 
 * File: geombase.h
 * Synopsis:
 *      Code that is common to (and included by) geometry2d.h 
 *      and geometry3d.h
 *
 * Written by
 *       Shubin Zhao (shubinz@cs.nyu.edu) (2001)
 *
 * WWW URL: http://cs.nyu.edu/exact/
 * Email: exact@cs.nyu.edu
 *
 * $Id$
 *****************************************************************/

#ifndef CORE_GEOMETRY_H
#define CORE_GEOMETRY_H

#include <CGAL/CORE/CORE.h>


//base class for geom2d and geom3d classes
class GeomObj {

public:

  // Exceptions

  class Exception {
  	public:
  	  virtual void print_message( char* msg ) { std::cerr << msg <<std::endl; }
  };

  class NoIntersection : public Exception { };

  class IllegalOperation : public Exception { };

  virtual int dim() const { return -1; }

}; //class GeomObj

#endif
