// Copyright (c) 1997-2001  
// ETH Zurich (Switzerland).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
// Author(s)     : Sven Schoenherr <sven@inf.ethz.ch>, Bernd Gaertner

namespace CGAL {

// Class implementation (continued)
// ================================

// I/O
// ---
template < class K_ >
std::ostream&
operator << ( std::ostream& os, const CGAL::Optimisation_ellipse_2<K_>& e)
{
    const char* const  empty       = "";
    const char* const  pretty_head = "CGAL::Optimisation_ellipse_2( ";
    const char* const  pretty_sep  = ", ";
    const char* const  pretty_tail = ")";
    const char* const  ascii_sep   = " ";

    const char*  head = empty;
    const char*  sep  = empty;
    const char*  tail = empty;

    switch ( CGAL::get_mode( os)) {
      case CGAL::IO::PRETTY:
        head = pretty_head;
        sep  = pretty_sep;
        tail = pretty_tail;
        break;
      case CGAL::IO::ASCII:
        sep  = ascii_sep;
        break;
      case CGAL::IO::BINARY:
        break;
      default:
        CGAL_optimisation_assertion_msg( false,
                                         "CGAL::get_mode( os) invalid!");
        break; }

    os << head << e.n_boundary_points;
    switch ( e.n_boundary_points) {
      case 0:
        break;
      case 1:
        os << sep << e.boundary_point1;
        break;
      case 2:
        os << sep << e.boundary_point1
           << sep << e.boundary_point2;
        break;
      case 3:
      case 5:
        os << sep << e.conic1;
        break;
      case 4:
        os << sep << e.conic1
           << sep << e.conic2;
        break; }
    os << tail;

    return( os);
}

template < class K_ >
std::istream&
operator >> ( std::istream& is, CGAL::Optimisation_ellipse_2<K_>& e)
{
    switch ( CGAL::get_mode( is)) {

      case CGAL::IO::PRETTY:
        std::cerr << std::endl;
        std::cerr << "Stream must be in ascii or binary mode" << std::endl;
        break;

      case CGAL::IO::ASCII:
      case CGAL::IO::BINARY:
        CGAL::read( is, e.n_boundary_points);
        switch ( e.n_boundary_points) {
          case 0:
            break;
          case 1:
            is >> e.boundary_point1;
            break;
          case 2:
            is >> e.boundary_point1
               >> e.boundary_point2;
            break;
          case 3:
          case 5:
            is >> e.conic1;
            break;
          case 4:
            is >> e.conic1
               >> e.conic2;
            break; }
        break;

      default:
        CGAL_optimisation_assertion_msg( false,
                                         "CGAL::get_mode( is) invalid!");
        break; }

    return( is);
}

} //namespace CGAL

// ===== EOF ==================================================================
