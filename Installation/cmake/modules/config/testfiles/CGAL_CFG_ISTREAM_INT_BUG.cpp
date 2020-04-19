// Copyright (c) 2006  GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
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
