// Copyright (c) 2002  Max Planck Institut fuer Informatik (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
// Author(s)     : Radu Ursu

#include "cgal_types.h"
#include<vector>
#include<boost/shared_ptr.hpp>

class FigureBase
{
protected:

  FigureBase ( CGAL::Color aColor ) : mColor(aColor), mIsVisible(true) {}
  
  virtual void Render ( CGAL::Qt_widget *widget ) = 0 ;
  
public:

  void render ( CGAL::Qt_widget *widget ) 
  {
    if ( mIsVisible )
    {
      *widget << mColor ;
      Render(widget);
    }
  }
  
  
  virtual ~FigureBase() {}
  
  bool is_visible() const { return mIsVisible ; }
  
  void set_is_visible( bool aState ) { mIsVisible = aState ; }
  
protected:  

  CGAL::Color mColor ;
  bool        mIsVisible ;
} ;
typedef boost::shared_ptr<FigureBase> FigureBasePtr ;

template<class T>
class Figure : public FigureBase
{
public:
  
  Figure( CGAL::Color aColor, T const& aF ) : FigureBase(aColor), mF(aF) {}
  
  virtual void Render ( CGAL::Qt_widget *widget ) { *widget << mF ; }
  
private:

  T mF ;  
} ;

typedef std::vector<FigureBasePtr> FigureList ;
