// Copyright (c) 2006  GeometryFactory (France).
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
// Author(s)     : Andreas Fabri

//| This flag is set, if the executable does not properly parse an int followed by a comma
//| Can you believe it!!!

#include <iostream>
#include <sstream>

int main()
{
  int i;
  std::string one = "7812,"; 
  std::istringstream is(one);
  is >> i;
  std::cout << i << std::endl;
  if(i != 7812){ // bad luck if there is garbage that equals 7812. 
    return 1;
  }
  return 0;

}
