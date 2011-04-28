// Copyright (c) 2005  Stanford University (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
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
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#include <CGAL/basic.h>

#include <CGAL/Kinetic/IO/internal/Qt_window_2.h>
#include <CGAL/Kinetic/IO/internal/pixmaps.h>
#include "Kinetic_Qt_window_2.moc"

namespace CGAL { namespace Kinetic { namespace internal {
        
Qt_window_2::Qt_window_2(int xmin, int xmax, int ymin, int ymax) {
  widget_= new Qt_widget_2_core(this);
  setCentralWidget(widget_);
  resize(800,800);
  widget()->show();
  widget()->set_window(xmin, xmax, ymin, ymax);

  //How to attach the standard toolbar
  _std_toolbar = new CGAL::Qt_widget_standard_toolbar(widget(), this, "toolbar");
  this->addToolBar(_std_toolbar->toolbar(), Top, FALSE);

  /*connect(_widget, SIGNAL(custom_redraw()),
    this, SLOT(redraw_win()) );*/

  QToolButton *play_button;         //the toolbar button
  play_button =  new QToolButton(QPixmap( (const char**)play_xpm ),
				 "Play",
				 0,
				 &core_,
				 SLOT(play_button()),
				 _std_toolbar,
				 "Play");

  QToolButton *pause_button;        //the toolbar button
  pause_button =  new QToolButton(QPixmap( (const char**)pause_xpm ),
				  "Pause",
				  0,
				  &core_,
				  SLOT(pause_button()),
				  _std_toolbar,
				  "Pause");

  QToolButton *stop_button;         //the toolbar button
  stop_button =  new QToolButton(QPixmap( (const char**)stop_xpm ),
				 "Stop",
				 0,
				 &core_,
				 SLOT(stop_button()),
				 _std_toolbar,
				 "Stop");

  QToolButton *play_to_button;      //the toolbar button
  play_to_button =  new QToolButton(QPixmap( (const char**)play_to_xpm ),
				    "Play to",
				    0,
				    &core_,
				    SLOT(play_to_button()),
				    _std_toolbar,
				    "Play to");

  QToolButton *play_through_button; //the toolbar button
  play_through_button =  new QToolButton(QPixmap( (const char**)play_through_xpm ),
					 "Play through",
					 0,
					 &core_,
					 SLOT(play_through_button()),
					 _std_toolbar,
					 "Play through");

  QToolButton *reverse_button;      //the toolbar button
  reverse_button =  new QToolButton(QPixmap( (const char**)reverse_xpm ),
				    "Reverse",
				    0,
				    &core_,
				    SLOT(reverse_button()),
				    _std_toolbar,
				    "Reverse");
  QToolButton *faster_button;       //the toolbar button
  faster_button =  new QToolButton(QPixmap( (const char**)faster_xpm ),
				   "Faster",
				   0,
				   &core_,
				   SLOT(faster_button()),
				   _std_toolbar,
				   "Faster");
  QToolButton *slower_button;       //the toolbar button
  slower_button =  new QToolButton(QPixmap( (const char**)slower_xpm ),
				   "Slower",
				   0,
				   &core_,
				   SLOT(slower_button()),
				   _std_toolbar,
				   "Slower");

  QToolButton *filePrintAction;
  filePrintAction = new QToolButton(QPixmap( (const char**)print_xpm ),
				    "Print", 0,
				    widget_,
				    SLOT(print_to_ps()),
				    _std_toolbar,
				    "Print");


}
} } } //namespace CGAL::Kinetic::internal
