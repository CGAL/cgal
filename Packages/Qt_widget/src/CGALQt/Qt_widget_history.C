// ======================================================================
//
// Copyright (c) 1997-2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// file          : src/CGALQt/Qt_widget_history.C
// package       : Qt_widget (1.3.21)
// maintainer    : Laurent Rineau <rineau@clipper.ens.fr>
// author(s)     : Laurent Rineau
// release       : 
// release_date  : 
//
// coordinator   : Laurent Rineau <rineau@clipper.ens.fr>
//
// ======================================================================

#ifdef CGAL_USE_QT

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
};

void Qt_widget_history::backward()
{
  if( (! history_list.empty()) && it!=history_list.begin() )
    {
      it--;
      restore();
    }
  emit(forwardAvaillable(true));
  if(it == history_list.begin())
    emit(backwardAvaillable(false));
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
  emit(backwardAvaillable(true));
  if( it == --history_list.end())
    emit(forwardAvaillable(false));
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
    emit(backwardAvaillable(true));
  emit(forwardAvaillable(false));
}

}; // end namespace
#include "Qt_widget_history.moc"

#endif // CGAL_USE_QT
