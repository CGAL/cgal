// Copyright (c) 2006, 2007, 2009  Stanford University (USA),
// INRIA Sophia-Antipolis (France).  All rights reserved.
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
// Author(s)     : Daniel Russel, Sylvain Pion

// Tests if BOOST_PROGRAM_OPTIONS is available.

#include <iostream>
#include <string>
#include <boost/program_options.hpp>

namespace po = boost::program_options;

int main(int ac, char *av[])
{
  std::cout << "version=" << BOOST_VERSION/100000 << "."
            << ((BOOST_VERSION / 100) % 100) << "."
            << BOOST_VERSION % 100 << std::endl;

  std::string str;

  po::options_description desc("Allowed options");
  desc.add_options()
    ("help", "produce help message")
    ("input-file,f", po::value<std::string>(&str)->default_value("blabla.txt"),
     "name of file")
    ;

  po::variables_map vm;
  po::store(po::parse_command_line(ac, av, desc), vm);
  po::notify(vm);    

  if (vm.count("help")) {
    std::cout << "Help" << "\n";
    return 1;
  }
  return 0;
}
