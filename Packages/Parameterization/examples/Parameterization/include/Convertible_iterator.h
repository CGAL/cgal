// Copyright (c) 2005  INRIA (France).
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
// $Revision$
// $Name$
//
// Author(s)     : Laurent Saboret, Bruno Levy, Pierre Alliez


#ifndef CGAL_CONVERTIBLE_ITERATOR_H
#define CGAL_CONVERTIBLE_ITERATOR_H


// Utility class for Mesh_adaptor_polyhedron_ex
// This class adds a conversion to handle/const handle to an iterator class
template <class Iter,			// base iterator
          class ConstHandle,	// const-handle type to convert to
          class Handle = void*>	// non-const-handle type to convert to 
								// (void* <=> none)
class Convertible_iterator : public Iter 
{
	typedef Iter											Base; 
	typedef Convertible_iterator<Iter, ConstHandle, Handle>	Self;

public:
	Convertible_iterator() {}
	Convertible_iterator(Base toCopy) : Base(toCopy) {}
	
	Convertible_iterator(const Self& toCopy) : Base(toCopy) {}
	Self& operator=(const Self& toCopy) { Base::operator=(toCopy); return *this; }

	Self & operator++()		{ Base::operator++(); return *this; }
	Self & operator--()		{ Base::operator--(); return *this; }
	Self operator++(int)	{ Self tmp(*this); ++(*this); return tmp; }
	Self operator--(int)	{ Self tmp(*this); --(*this); return tmp; }

	operator Handle()				{ return operator->(); }
	operator ConstHandle() const	{ return operator->(); }
};


#endif //CGAL_CONVERTIBLE_ITERATOR_H
