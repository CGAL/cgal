// Copyright (c) 1999  Max-Planck-Institute Saarbrucken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Stefan Schirra 

#ifndef CGAL_STL_EXTENSIONS_H
#define CGAL_STL_EXTENSIONS_H

#include <CGAL/basic.h>
#include <iterator>
#include <utility>

CGAL_BEGIN_NAMESPACE

template <class ForwardIterator>
inline
ForwardIterator
successor( ForwardIterator it )
{
 return ++it;
}

template <class BidirectionalIterator>
inline
BidirectionalIterator
predecessor( BidirectionalIterator it )
{
 return --it;
}

CGAL_END_NAMESPACE

#include <CGAL/IO/Tee_for_output_iterator.h>

#endif // CGAL_STL_EXTENSIONS_H
