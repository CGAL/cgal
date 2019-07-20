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
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Susan Hert <hert@mpi-sb.mpg.de>
#ifndef CGAL_CIRC_PAIR_H
#define CGAL_CIRC_PAIR_H

#include <CGAL/license/Partition_2.h>


#include <CGAL/enum.h>

namespace CGAL {

//
// "before" always means a point between front and back in the indicated
// direction
//
// "after" always means a point just beyond the range of [front,back]
//
template <class BidirectionalCirculator>
class Circ_pair 
{

public:
    Circ_pair(BidirectionalCirculator back, BidirectionalCirculator front) : 
        _front(front), _back(back), _direction(COUNTERCLOCKWISE) {}

    Circ_pair(BidirectionalCirculator back, BidirectionalCirculator front, 
              Orientation dir) : _front(front), _back(back), _direction(dir) {}

    Circ_pair(BidirectionalCirculator front_and_back, Orientation dir) : 
          _front(front_and_back), _back(front_and_back), _direction(dir) {}

    void initialize(BidirectionalCirculator new_back_and_front) 
    {
       _back = _front = new_back_and_front;
    }

    void push_back(BidirectionalCirculator new_back) 
    {
       _back = new_back;
    }

    void pop_back() 
    {
       _back = before_back();
    }

    void push_front(BidirectionalCirculator new_front)
    {
       _front = new_front;
    }

    void pop_front()
    {
       _front = before_front();
    }

    BidirectionalCirculator front() const
    {
       return _front;
    }

    BidirectionalCirculator back() const
    {
       return _back;
    }

    Orientation direction() const
    {
       return _direction;
    }

    void set_direction(Orientation direction) 
    {
       _direction = direction;
    }

    void change_dir()
    {
        if (_direction == CLOCKWISE)
           _direction = COUNTERCLOCKWISE;
        else
           _direction = CLOCKWISE;
    }

    BidirectionalCirculator before_back()
    {
        BidirectionalCirculator temp = _back;
        if (_direction == COUNTERCLOCKWISE)
           return ++temp;
        else
           return --temp;
    }
    BidirectionalCirculator after_back()
    {
        BidirectionalCirculator temp = _back;
        if (_direction == COUNTERCLOCKWISE)
           return --temp;
        else
           return ++temp;
    }
    BidirectionalCirculator before_front()
    {
        BidirectionalCirculator temp = _front;
        if (_direction == COUNTERCLOCKWISE)
           return --temp;
        else
           return ++temp;
    }
    BidirectionalCirculator after_front()
    {
        BidirectionalCirculator temp = _front;
        if (_direction == COUNTERCLOCKWISE)
           return ++temp;
        else
           return --temp;
    }
private:
    BidirectionalCirculator _front, _back;
    Orientation _direction;
};

}

#endif // CGAL_CIRC_PAIR_H
