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

//#ifdef CGAL_USE_SOQT

#include "include/control_coin.h"
#ifdef CGAL_USE_COIN 
#include "include/SoQt_examiner_viewer.h"
#include <qpushbutton.h>
#include <Inventor/nodes/SoSeparator.h>
#include <CGAL/Kinetic/IO/internal/pixmaps.h>

#include "SoQt_examiner_viewer.moc"

namespace CGAL
{
  namespace Kinetic
  {
    SoQt_examiner_viewer::SoQt_examiner_viewer(QWidget * parent):
      SoQtExaminerViewer(parent, NULL, TRUE,
			 SoQtFullViewer::BUILD_ALL,
			 SoQtFullViewer::BROWSER,
			 FALSE) {
      // Explicitly trigger the construction of viewer decorations.
      QWidget * widget = this->buildWidget(this->getParentWidget());
      this->setBaseWidget(widget);
      root_= new SoSeparator;
      this->setSceneGraph(root_);
    }

    void SoQt_examiner_viewer::new_subgraph(SoNode *p) {
      root_->addChild(p);
    }

    void SoQt_examiner_viewer::delete_subgraph(SoNode *p) {
      root_->removeChild(p);
    }

#define SETUP_QT_BUTTON(name)   name ## _button_ = new QPushButton(parent); \
	name##_button_->setFocusPolicy(QWidget::NoFocus);		\
	name##_button_->setPixmap(QPixmap((const char **) internal::name##_xpm)); \
	name##_button_->adjustSize();					\
	QObject::connect(name##_button_, SIGNAL(clicked()),		\
			 &core_, SLOT(name##_button()));		\
	buttonlist->append(name##_button_);

    void SoQt_examiner_viewer::createViewerButtons(QWidget * parent, SbPList * buttonlist) {
      SoQtExaminerViewer::createViewerButtons(parent, buttonlist);
      // [now add your own button(s) to the buttonlist]
      SETUP_QT_BUTTON(play);
      SETUP_QT_BUTTON(pause);
      SETUP_QT_BUTTON(stop);
      SETUP_QT_BUTTON(play_to);
      SETUP_QT_BUTTON(play_through);
      SETUP_QT_BUTTON(reverse);
      SETUP_QT_BUTTON(faster);
      SETUP_QT_BUTTON(slower);
    }

  };
}
//#else

//static bool SoQt_examiner_viewer_compiled_without_CGAL_USE_SOQT_defined;

//#endif
#endif // CGAL_USE_COIN
