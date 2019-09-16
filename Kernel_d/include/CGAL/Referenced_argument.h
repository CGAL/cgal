// Copyright (c) 2009-2014 INRIA Sophia-Antipolis (France).
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
// SPDX-License-Identifier: LGPL-3.0+
//
// Author(s)    : Samuel Hornus, Olivier Devillers

#ifndef REFERENCED_ARGUMENT_H
#define REFERENCED_ARGUMENT_H

template< typename T>
struct Referenced_argument
{
    Referenced_argument() : arg_() {}
    template<typename A>
    Referenced_argument(A a) : arg_(a) {}
    template<typename A, typename B>
    Referenced_argument(A a, B b) : arg_(a, b) {}
    T & arg() { return arg_; }

protected:
    T arg_;
};

#endif // REFERENCED_ARGUMENT_H
