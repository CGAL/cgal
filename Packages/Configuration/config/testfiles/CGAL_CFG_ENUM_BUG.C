// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 
//
// file          : config/testfiles/CGAL_CFG_ENUM_BUG.C
// source        :
// revision      : 1.11
// revision_date : 29 Mar 1998
// author(s)     : various
//
// coordinator   : Utrecht University
//
// ======================================================================

// CGAL_CFG_ENUM_BUG.C
// ---------------------------------------------------------------------
// This program is used by cgal_configure.
// The following documentation will be pasted in the generated configfile.
// ---------------------------------------------------------------------

//| This flag is set, if the compiler does not promote enumeration types
//| (which depend on a template parameter) correctly when they are used 
//| as int template arguments. (e.g. Borland 5.5)

struct F {
  enum { a = 1 };
};

template < int i > struct B;
template <> struct B< 1 > {};

template < class T >
struct C {
  enum { ar = T::a };
  B< ar > b;
  // gives
  // Borland C++ 5.5.1 for Win32 Copyright (c) 1993, 2000 Borland
  // Error E2450 Undefined structure 'B<0>' in function main()
  
  // using 
  // B< T::a > b;
  // instead gives
  // Borland C++ 5.5.1 for Win32 Copyright (c) 1993, 2000 Borland
  // Error E2401 Invalid template argument list
  // Error E2040 Declaration terminated incorrectly
};

int main()
{
  C< F > c;
  return 0;
}

// EOF //
