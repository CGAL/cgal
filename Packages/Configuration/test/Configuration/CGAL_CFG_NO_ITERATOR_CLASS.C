// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 
//
// file          : config/testfiles/CGAL_CFG_NO_TYPENAME.C
// source        :
// revision      : 1.0
// revision_date : 29 June 1999
// author(s)     : Geert-Jan Giezeman
//
// coordinator   : Utrecht University
//
// ======================================================================

//| Some STL libraries don't supply the class iterator yet.
//| These libraries usually supply the old and non-standard HP classes
//| input_iterator, forward_iterator, etc.


#include <iterator>

class my_forward_iterator :
	public std::iterator<std::forward_iterator_tag,double>
	//public std::forward_iterator<double,ptrdiff_t>
{};

int main()
{
    my_forward_iterator it;
    return 0;
}

