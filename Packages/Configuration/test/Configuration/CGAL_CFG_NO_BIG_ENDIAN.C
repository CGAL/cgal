// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 
//
// file          : config/testfiles/CGAL_CFG_NO_BIG_ENDIAN.C
// source        :
// revision      : 1.11
// revision_date : 29 Mar 1998
// author(s)     : various
//
// coordinator   : Utrecht University
//
// ======================================================================

// CGAL_CFG_NO_BIG_ENDIAN.C
// ---------------------------------------------------------------------
// A short test program to evaluate a machine architecture.
// This program is used by cgal_configure.
// The following documentation will be pasted in the generated configfile.
// ---------------------------------------------------------------------

//| The byte order of a machine architecture distinguishes into
//| big-endian and little-endian machines.
//| The following definition is set if it is a little-endian machine.

union {
  int       testWord;
  char      testByte[sizeof(int)];
} endianTest;

int main() {
    endianTest.testWord = 1;
    if (endianTest.testByte[0] == 1)
	return 1;  // little-endian
    return 0;      // big-endian
}

// EOF //

