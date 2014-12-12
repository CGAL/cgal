// Copyright (c) 2014
// INRIA Saclay-Ile de France (France)
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
// Author(s)     : Marc Glisse

#ifndef CGAL_STATIC_INT_H
#define CGAL_STATIC_INT_H
#include <CGAL/constant.h>

namespace CGAL {
template <class NT> struct static_zero {
	operator NT() const { return constant<NT,0>(); }
};
template <class NT> struct static_one {
	operator NT() const { return constant<NT,1>(); }
};

template <class NT> static_zero<NT> operator-(static_zero<NT>) { return static_zero<NT>(); }

template <class NT> NT operator+(NT const& x, static_zero<NT>) { return x; }
template <class NT> NT operator+(static_zero<NT>, NT const& x) { return x; }
template <class NT> static_zero<NT> operator+(static_zero<NT>, static_zero<NT>) { return static_zero<NT>(); }
template <class NT> static_one<NT> operator+(static_zero<NT>, static_one<NT>) { return static_one<NT>(); }
template <class NT> static_one<NT> operator+(static_one<NT>, static_zero<NT>) { return static_one<NT>(); }

template <class NT> NT operator-(NT const& x, static_zero<NT>) { return x; }
template <class NT> NT operator-(static_zero<NT>, NT const& x) { return -x; }
template <class NT> static_zero<NT> operator-(static_zero<NT>, static_zero<NT>) { return static_zero<NT>(); }
template <class NT> static_zero<NT> operator-(static_one<NT>, static_one<NT>) { return static_zero<NT>(); }
template <class NT> static_one<NT> operator-(static_one<NT>, static_zero<NT>) { return static_one<NT>(); }

template <class NT> NT operator*(NT const& x, static_one<NT>) { return x; }
template <class NT> NT operator*(static_one<NT>, NT const& x) { return x; }
template <class NT> static_zero<NT> operator*(NT const&, static_zero<NT>) { return static_zero<NT>(); }
template <class NT> static_zero<NT> operator*(static_zero<NT>, NT const&) { return static_zero<NT>(); }
template <class NT> static_zero<NT> operator*(static_zero<NT>, static_zero<NT>) { return static_zero<NT>(); }
template <class NT> static_one<NT> operator*(static_one<NT>, static_one<NT>) { return static_one<NT>(); }
template <class NT> static_zero<NT> operator*(static_zero<NT>, static_one<NT>) { return static_zero<NT>(); }
template <class NT> static_zero<NT> operator*(static_one<NT>, static_zero<NT>) { return static_zero<NT>(); }

template <class NT> NT operator/(NT const& x, static_one<NT>) { return x; }
template <class NT> static_zero<NT> operator/(static_zero<NT>, NT const&) { return static_zero<NT>(); }
template <class NT> static_zero<NT> operator/(static_zero<NT>, static_one<NT>) { return static_zero<NT>(); }
template <class NT> static_one<NT> operator/(static_one<NT>, static_one<NT>) { return static_one<NT>(); }

}
#endif // CGAL_STATIC_INT_H
