// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 
//
// file          : config/testfiles/CGAL_CFG_NO_DYNAMIC_CAST.C
// source        :
// revision      : 1.11
// revision_date : 29 Mar 1998
// author(s)     : various
//
// coordinator   : Utrecht University
//
// ======================================================================

// CGAL_CFG_NO_DYNAMIC_CAST.C
// ---------------------------------------------------------------------
// A short test program to evaluate a C++ compiler.
// This program is used by cgal_configure.
// The following documentation will be pasted in the generated configfile.
// ---------------------------------------------------------------------

//| This flag is set if the compiler doesn't support the operator dynamic_cast.

#ifdef __GNUG__
#include <typeinfo>
#endif // __GNUG__


class Base
{
  public:
    virtual ~Base() {}
};

class Derived : public Base
{ };

bool all_assertions_correct = true;

int main()
{
  Base *p_Base_Derived = new Derived;
  Base *p_Base_Base = new Base;

  Derived* p_Derived_Derived = dynamic_cast<Derived *>(p_Base_Derived);
  all_assertions_correct &= (p_Derived_Derived != 0);
  Derived* p_Derived_Base = dynamic_cast<Derived *>(p_Base_Base);
  all_assertions_correct &= (p_Derived_Base == 0);

  return !all_assertions_correct;
}

// EOF //
