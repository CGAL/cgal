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
// release       : $CGAL_Revision: V $
// release_date  : $CGAL_Date: 1998, July 17 $
//
// file          : config/testfiles/CGAL_CFG_INCOMPLETE_TYPE_BUG_1.C
// source        :
// revision      : 1.11
// revision_date : 29 Mar 1998
// author(s)     : various
//
// coordinator   : Utrecht University
//
// ======================================================================

// CGAL_CFG_INCOMPLETE_TYPE_BUG_1.C
// ---------------------------------------------------------------------
// A short test program to evaluate a C++ compiler.
// This program is used by cgal_configure.
// The following documentation will be pasted in the generated configfile.
// ---------------------------------------------------------------------

//| When a class (Cls_P / Cls_Pl) refers to a not yet defined class (Cls_At),
//| some compilers give an "incomplete type error".
//| This program is used to detect a special case of this problem
//| where templates are involved.

template <class T> class Cls_P;
template <class T> class Cls_Pl;
template <class T> class Cls_At;

template <class T>
class Cls_P
{
 public:
             Cls_P(const T& t);

    Cls_P<T> transform(const Cls_At<T>& a) const;

    T  _T;
};

template <class T>
Cls_P<T>::Cls_P(const T& t) : _T(t)
{}

template <class T>
Cls_P<T>
Cls_P<T>::transform(const Cls_At<T>& a) const
{
 return a.transform( *this );
}

template <class T>
class Cls_Pl
{
 public:
             Cls_Pl(const Cls_P<T>& t);

    Cls_Pl<T> transform(const Cls_At<T>& a) const;
    Cls_P<T>  point() const;

    Cls_P<T>  _P;
};

template <class T>
Cls_Pl<T>::Cls_Pl(const Cls_P<T>& t) : _P(t)
{}

template <class T>
Cls_Pl<T>
Cls_Pl<T>::transform(const Cls_At<T>& a) const
{
 return a.transform( *this );
}

template <class T>
Cls_P<T>
Cls_Pl<T>::point() const
{
 return _P;
}

template <class T>
class Cls_At
{
 public:
      Cls_At(const T& t);
 
      Cls_P<T>    transform(const Cls_P<T>& p) const;
      Cls_Pl<T>   transform(const Cls_Pl<T>& pl) const;

   T  _aT;
};

template <class T>
Cls_At<T>::Cls_At(const T& t) : _aT(t)
{}

template <class T>
Cls_P<T>
Cls_At<T>::transform(const Cls_P<T>& p) const
{
 return p;
}

template <class T>
Cls_Pl<T>
Cls_At<T>::transform(const Cls_Pl<T>& pl) const
{
 Cls_P<T>  hp = pl.point();
 return pl;
}

 
int main()
{
  // The presence of the following line makes a difference
  // Cls_Pl<int>*  uninitptr;     
  int i1 = 12;
  int i2 =  8;
  int i3 = 17;

  Cls_P<int>  p1(i1);
  Cls_P<int>  p2(i2);

  Cls_Pl<int> eb1( p1 );
  Cls_Pl<int> eb2( p2 );

  Cls_At<int> at( i3 );

  eb2 = eb1.transform( at );
  eb2 = at.transform( eb1 );

  return 0;
}  

// EOF //
