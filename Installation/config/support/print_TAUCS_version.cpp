// Copyright (c) 2004  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
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
// $URL: svn+ssh://fcacciola@scm.gforge.inria.fr/svn/cgal/trunk/Installation/config/support/test_TAUCS.cpp $
// $Id: test_TAUCS.cpp 41686 2008-01-18 20:33:57Z spion $
//
// Author(s)     : Laurent Saboret

// Test if TAUCS is available


#include <iostream>


int main(int argc, char* argv[])
{
    // TAUCS provides no version number :-(
    // Version 1 is obsolete, thus we assume version 2 (latest is 2.2 on 03/2006)
    std::cout << "version=2.x" << std::endl;

    return 0;
}
