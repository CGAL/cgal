// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// ----------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : config/testfiles/CGAL_CFG_DYNAMIC_CAST_BUG.C
// source        :
// revision      : 1.11
// revision_date : 29 Mar 1998
// author(s)     : various
//
// coordinator   : Utrecht University
//
// ======================================================================

// CGAL_CFG_DYNAMIC_CAST_BUG.C
// ---------------------------------------------------------------------
// A short test program to evaluate a C++ compiler.
// This program is used by cgal_configure.
// The following documentation will be pasted in the generated configfile.
// ---------------------------------------------------------------------

//| When a dynamic cast involves a pointer to a not yet instantiated 
//| template class, some compilers give an error.
//| This program is used to detect this problem.

class B
{
public:
  virtual ~B() {}
};

template <class T>
class W : public B
{
public:
  T obj;
};

template <class T>
class L
{
public:
  T nt;
};

template <class T>
class P
{
public:
  T nt;
};

int
main()
{
  W< P<int> > wp;
  W< L<int> >* wl_ptr = dynamic_cast<W< L<int> >* >( &wp );
  (void) wl_ptr;

  return 0;
}
  
//| "CGAL_CFG_DYNAMIC_CAST_BUG.C", line 45: error(3105): the type in a
//|           dynamic_cast must be a pointer or reference to a complete class
//|           type, or void *
//|     W< L<int> >* wl_ptr = dynamic_cast<W< L<int> >* >( &wp );
//|                                        ^

