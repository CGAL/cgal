// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 
//
// file          : config/testfiles/CGAL_CFG_EXPLICIT_SPECIALIZATION_MEMBER_DECLARATION_BUG.C
// package       : Configuration (1.28)
// source        :
// revision      : 1.11
// revision_date : 29 Mar 1998
// author(s)     : various
//
// coordinator   : Utrecht University
//
// ======================================================================

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


