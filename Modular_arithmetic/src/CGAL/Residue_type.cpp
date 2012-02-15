// Copyright (c) Max-Planck-Institute Saarbruecken (Germany).
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
// Author(s)     : Michael Hemmer


#include <CGAL/Residue.h>

namespace CGAL{
#ifdef CGAL_HAS_THREADS
boost::thread_specific_ptr<int>    Residue::prime_int_;
boost::thread_specific_ptr<double> Residue::prime_;
boost::thread_specific_ptr<double> Residue::prime_inv_;
#else
int Residue::prime_int = 67111067;
double Residue::prime = 67111067.0;
double Residue::prime_inv =1/67111067.0;

#endif

const double Residue::CST_CUT = std::ldexp( 3., 51 );

}
