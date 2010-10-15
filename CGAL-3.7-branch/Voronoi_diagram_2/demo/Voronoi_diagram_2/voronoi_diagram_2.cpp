// Copyright (c) 2006 Foundation for Research and Technology-Hellas (Greece).
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
// $URL$
// $Id$
//
//
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>

#include <CGAL/basic.h>

#include <iostream>
#include <fstream>

#include <qapplication.h>
#include <qmainwindow.h>

#include <CGAL/Timer.h>

#include <list>
#include <vector>

#include "typedefs.h"

//************************************
// global variables
//************************************

int num_selected;
std::vector<Site_2> sitelist;

#include "my_window.h"

#include "qt_file_toolbar.moc"
#include "qt_layers_toolbar.moc"
#include "my_window.moc"



int
main(int argc, char* argv[])
{

  int size = 750;

  QApplication app( argc, argv );

  My_Window W(size,size);
  app.setMainWidget( &W );
#if !defined (__POWERPC__)
  QPixmap cgal_icon = QPixmap((const char**)demoicon_xpm);
  W.setIcon(cgal_icon);
#endif
  W.show();
  W.set_window(0,size,0,size);
  W.setCaption(W.title());


  return app.exec();
}
