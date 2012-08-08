// Copyright (c) 2002-2004  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Radu Ursu and Laurent Rineau

#include <CGAL/basic.h>


#include <CGAL/IO/Qt_widget_history.h>

namespace CGAL {

Qt_widget_history::Qt_widget_history(Qt_widget* parent, const char* name):
  QObject(parent, name), widget(parent)
{
  it = history_list.begin();

  connect(widget, SIGNAL(rangesChanged()),
	  this, SLOT(save()));
  
  // backward compatibility with CGAL-2.4
  connect(parent, SIGNAL(internal_back()),
	  this, SLOT(backward()));
  connect(parent, SIGNAL(internal_forth()),
	  this, SLOT(forward()));
  connect(parent, SIGNAL(internal_add_to_history()),
	  this, SLOT(save()));
  connect(parent, SIGNAL(internal_clear_history()),
	  this, SLOT(clear()));
}

void Qt_widget_history::backward()
{
  if( (! history_list.empty()) && it!=history_list.begin() )
    {
      it--;
      restore();
    }
  emit forwardAvaillable(true);
  if(it == history_list.begin())
    emit backwardAvaillable(false);
}

void Qt_widget_history::forward()
{
  if( it != history_list.end() ) 
  {
    if( ++it != history_list.end() )
      restore();
    else 
      --it;
  }
  emit backwardAvaillable(true);
  if(it == --history_list.end())
    emit forwardAvaillable(false);
}

void Qt_widget_history::save()
{
  if( it != history_list.end() )
  {
    ++it;
    history_list.erase(it, history_list.end());
  }
  
  History_atom* atom = new History_atom();
  atom->save(*widget);
  history_list.push_back(atom);
  
  it = history_list.end();
  it--;

  if( it != history_list.begin() )
    emit backwardAvaillable(true);
  emit forwardAvaillable(false);
}

} // end namespace
#include "Qt_widget_history.moc"

