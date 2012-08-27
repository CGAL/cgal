// Copyright (c) 2005  Stanford University (USA).
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
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_KINETIC_IO_INTERNAL_COIN_SIMULATOR_GUI_H
#define CGAL_KINETIC_IO_INTERNAL_COIN_SIMULATOR_GUI_H

#include <CGAL/Kinetic/basic.h>

//#ifdef HAVE_INVENTOR_QT_SOQT_H
//#include <Inventor/Qt/SoQt.h>
#include <Inventor/Qt/viewers/SoQtExaminerViewer.h>
#include <Inventor/nodes/SoNode.h>
//#include <Inventor/nodes/SoSeparator.h>
//#else
//#ifdef HAVE_COIN2_INVENTOR_QT_SOQT_H
//#include <Coin2/Inventor/Qt/SoQt.h>
//#include <Coin2/Inventor/Qt/viewers/SoQtExaminerViewer.h>
//#include <Coin2/Inventor/nodes/SoSeparator.h>
//#endif
//#endif

#include <CGAL/Kinetic/IO/internal/Qt_core.h>
#include <qapplication.h>
//#include <qmainwindow.h>
#include <map>

class QPushButton;
class SoSeparator;

// I think I need these here explicitly for MOC to work
namespace CGAL
{
  namespace Kinetic
  {
    /*  Usage main_window_= SoQt::init(argc, argv, argv[0]);
	viewer_= new Coin_simulator_viewer(main_window_);
	SoQt::show(main_window_); SoQt::mainLoop();

	Note that for some reason QObject must come before
	SoQtExaminerViewer in the inheritence list.
    */
    class SoQt_examiner_viewer : public QObject, public SoQtExaminerViewer
    {
      Q_OBJECT
    public:
      SoQt_examiner_viewer(QWidget * parent);

      virtual ~SoQt_examiner_viewer(){}

      void new_subgraph(SoNode *p);

      void delete_subgraph(SoNode *p);

      typedef internal::Qt_core Button_handler;
      Button_handler *button_handler() {
	return &core_;
      }
    protected:

      virtual void createViewerButtons(QWidget * parent, SbPList * buttonlist);

      SoSeparator* root_;
      QPushButton *play_button_, *pause_button_, *stop_button_;
      QPushButton *play_to_button_, *play_through_button_, *reverse_button_;
      QPushButton *faster_button_, *slower_button_;
      internal::Qt_core core_;
    };

  };
};
#endif                                            // guard
