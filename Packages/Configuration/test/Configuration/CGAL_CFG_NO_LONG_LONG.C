// ======================================================================
//
// Copyright (c) 2002 The CGAL Consortium
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 
//
// file          : config/testfiles/CGAL_CFG_NO_LONG_LONG.C
// source        :
// revision      : 
// revision_date : 
// author(s)     : Sylvain Pion
//
// coordinator   : Utrecht University
//
// ======================================================================

// CGAL_CFG_NO_LONG_LONG.C
// ---------------------------------------------------------------------
// A short test program to evaluate a machine architecture.
// This program is used by cgal_configure.
// The following documentation will be pasted in the generated configfile.
// ---------------------------------------------------------------------

//| The long long built-in integral type is not part of the ISO C++ standard,
//| but many compilers support it nevertheless since it's part of the ISO
//| C standard.
//| The following definition is set if it is supported.

int main()
{
    long long int i = 1;
    (void) i;
    return 0;
}

