// Copyright (c) 2000  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
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
// Author(s)     : Susan Hert <hert@mpi-sb.mpg.de>

#ifndef CGAL_PARTITION_OPT_CVX_EDGE_H
#define CGAL_PARTITION_OPT_CVX_EDGE_H

#include <iostream>
#include <CGAL/Partition_2/Partition_opt_cvx_diagonal_list.h>

namespace CGAL {

enum Partition_opt_cvx_edge_validity {PARTITION_OPT_CVX_NOT_VALID, 
                                      PARTITION_OPT_CVX_START_VALID,
                                      PARTITION_OPT_CVX_END_VALID, 
                                      PARTITION_OPT_CVX_BOTH_VALID};

class Partition_opt_cvx_edge 
{
public:
   Partition_opt_cvx_edge(): _is_done(false), 
                             _validity(PARTITION_OPT_CVX_NOT_VALID), 
                             _is_visible(false), _value(0) {}

   bool is_done( ) const { return _is_done; }

   bool is_valid( ) const { return _validity != PARTITION_OPT_CVX_NOT_VALID; }

   Partition_opt_cvx_edge_validity validity( ) const { return _validity; }

   bool is_visible( ) const { return _is_visible; }

   int  value( ) const { return _value; }

   Partition_opt_cvx_diagonal_list  solution( ) const { return _solution; }

   void set_done(bool val)  { _is_done = val; }

   void set_valid(Partition_opt_cvx_edge_validity val) { _validity = val; }
   
   template <class Point_2_, class Traits>
   void set_valid(const Point_2_& p1, const Point_2_& p2, const Point_2_& p3,
                  const Point_2_& p4, const Point_2_& p5, const Point_2_& p6,
                  const Traits& traits)
   { 
      typedef typename Traits::Left_turn_2     Left_turn_2;
      Left_turn_2 left_turn = traits.left_turn_2_object();

      _validity = PARTITION_OPT_CVX_NOT_VALID;

      Turn_reverser<Point_2_, Left_turn_2>  right_turn(left_turn);
      if (right_turn(p1, p2, p3))
         _validity = PARTITION_OPT_CVX_START_VALID;
      if (right_turn(p4, p5, p6)) {
         if (_validity == PARTITION_OPT_CVX_START_VALID) 
            _validity = PARTITION_OPT_CVX_BOTH_VALID;
         else _validity = PARTITION_OPT_CVX_END_VALID;
      }
   }

   void set_visible(bool vis) { _is_visible = vis; }

   void set_value(int val) { _value = val; }

   void set_solution(const Partition_opt_cvx_diagonal_list& diag_list) 
   { _solution = diag_list; }


private:
   bool                            _is_done;
   Partition_opt_cvx_edge_validity _validity;
   bool                            _is_visible;
   int                             _value;
   Partition_opt_cvx_diagonal_list _solution;
};

inline
std::ostream& operator<<(std::ostream& os, const Partition_opt_cvx_edge& g)
{
   if (g.is_done()) os << "1";
   else os << "0";
   if (g.validity() == PARTITION_OPT_CVX_NOT_VALID) os << "0";
   else if (g.validity() == PARTITION_OPT_CVX_START_VALID) os << "1";
   else if (g.validity() == PARTITION_OPT_CVX_END_VALID) os << "2";
   else os << "3";
   if (g.is_visible()) os << "1";
   else os << "0";
   os << g.value();
   return os;
}

}

#endif // CGAL_PARTITION_OPT_CVX_EDGE_H
