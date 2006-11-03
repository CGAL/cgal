// Copyright (c) 2006 Inria Lorraine (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
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
// Author(s)     : Luis Peñaranda <penarand@loria.fr>

CGAL_BEGIN_NAMESPACE

#ifndef CGAL_ALGEBRAIC_1_IMPL_H
#define CGAL_ALGEBRAIC_1_IMPL_H

template<class T>bool Algebraic_1::operator==(const T &n2)const{
	if(contains(n2))
		if(is_point())
			return true;
		else
			overlap();
	return false;
};

template<class T>bool Algebraic_1::operator!=(const T &n2)const{
	return !(operator==(n2));
};

template<class T>bool Algebraic_1::operator<=(const T &n2)const{
	return ((operator==(n2))||(operator<(n2)));
};

template<class T>bool Algebraic_1::operator>=(const T &n2)const{
	return ((operator==(n2))||(operator>(n2)));
};

template<class T>bool operator==(const T &n1,const Algebraic_1 &n2){
	return (n2==n1);
}

template<class T>bool operator!=(const T &n1,const Algebraic_1 &n2){
	return (n2!=n1);
}

template<class T>bool operator<(const T &n1,const Algebraic_1 &n2){
	return (n2>n1);
}

template<class T>bool operator>(const T &n1,const Algebraic_1 &n2){
	return (n2<n1);
}

template<class T>bool operator<=(const T &n1,const Algebraic_1 &n2){
	return (n2>=n1);
}

template<class T>bool operator>=(const T &n1,const Algebraic_1 &n2){
	return (n2<=n1);
}

#endif // CGAL_ALGEBRAIC_1_IMPL_H

CGAL_END_NAMESPACE
