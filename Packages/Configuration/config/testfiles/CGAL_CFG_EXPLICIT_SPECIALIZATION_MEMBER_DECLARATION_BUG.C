// Copyright (c) 1997  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbrucken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : various

// CGAL_CFG_EXPLICIT_SPECIALIZATION_MEMBER_DECLARATION_BUG
// ---------------------------------------------------------------------
// A short test program to evaluate a C++ compiler / STL implementation
// whether it supports the new standard headers (i.e. without the .h suffix)
// This program is used by cgal_configure.
// The following documentation will be pasted in the generated configfile.
// ---------------------------------------------------------------------

//| The flag CGAL_CFG_EXPLICIT_SPECIALIZATION_MEMBER_DECLARATION_BUG is set,
//| if a compiler does not support the declaration of the specialized member.
//| This is useful when you need to put the definition in a .C file.
//| The solution for VC7(cl1300) is to remove the declaration and to make the
//| member inline to avoid linker errors (more definitions in multiple units).

#include <iostream>

template <class T>
struct A {
 A(T t)
  {
    std::cout << "In A<T>(const T& d)" << std::endl;
  }
};

template<>
A<double>::A(double d);

template<>
A<double>::A(double d)
{
  std::cout << "In A<double>(const double& d)" << std::endl;
} 

int main(){
  A<int> a(10);
  return 0;
}

//THIS WORKS:
//template <class T>
//struct A {
// A(T t)
//  {
//    std::cout << "In A<T>(const T& d)" << std::endl;
//  }
//};
//
//
//#include <iostream>
//
//template<>
//inline A<double>::A(double d)
//{
//  std::cout << "In A<double>(const double& d)" << std::endl;
//} 


