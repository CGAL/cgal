// Copyright (c) 2008  GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
// Author(s)     : Andreas Fabri <Andreas.Fabri@geometryfactory.com>
//                 Laurent Rineau <Laurent.Rineau@geometryfactory.com>

#ifndef CGAL_QT_DEMOS_MAIN_WINDOW_H
#define CGAL_QT_DEMOS_MAIN_WINDOW_H
#include <iostream>
#include <QVector>
#include <QMainWindow>
#include <QDragEnterEvent>
#include <QDropEvent>
#include <CGAL/auto_link/Qt4.h>
#include <CGAL/export/Qt4.h>
#include <CGAL/Qt/resources.h>

// forward declaration
class QLabel;
class QGraphicsView;
class QAction;
class QMenu;

namespace CGAL {
namespace Qt {

// forward declaration
class GraphicsViewNavigation;

class CGAL_QT4_EXPORT DemosMainWindow : public QMainWindow 
{
  Q_OBJECT

public:
  enum Option { NoOption = 0x0,
		Size     = 0x1,
		Position = 0x2,
		State    = 0x4};

  Q_DECLARE_FLAGS(Options, Option)

  virtual void dragEnterEvent(QDragEnterEvent *event);
  virtual void dropEvent(QDropEvent *event);

  virtual void open(QString)
  {
    std::cerr << "You should implement open(QString);" << std::endl; 
  }

public:
  unsigned int maxNumberOfRecentFiles() const ;

public slots:
  void setMaxNumberOfRecentFiles(const unsigned int);

private:
  QMenu* getMenu(QString objectName, QString title);
  void popupAboutBox(QString title, QString html_resource_name);
  QMenu* getHelpMenu();

protected:
  DemosMainWindow (QWidget * parent = 0, ::Qt::WindowFlags flags = 0 );
  void setupStatusBar();
  void addNavigation(QGraphicsView*);
  void setupOptionsMenu(QMenu* menu  = 0);
  void addAboutCGAL(QMenu* menu  = 0);
  void addAboutDemo(QString htmlResourceName, QMenu* menu  = 0);

  void addRecentFiles(QMenu* menu, QAction* insertBefore = 0);

  void writeState(QString groupname = "MainWindow");
  void readState(QString groupname = "MainWindow",
		 Options what_to_save = Options(Size|State));

protected slots:
  void setUseAntialiasing(bool checked);
  void setUseOpenGL(bool checked);
  void popupAboutCGAL();
  void popupAboutDemo();

  void openRecentFile_aux();
  void addToRecentFiles(QString fileName);
  void updateRecentFileActions();

signals:
  void openRecentFile(QString filename);

protected:
  QGraphicsView* view;
  GraphicsViewNavigation* navigation;
  QLabel* xycoord ;

  QAction *actionUse_OpenGL;
  QAction *actionUse_Antialiasing;
  QAction *actionAbout;
  QAction *actionAboutCGAL;

  QString aboutHtmlResource;

  QAction* recentFilesSeparator;
  unsigned int maxNumRecentFiles;
  QVector<QAction*> recentFileActs;
}; // end class DemosMainWindow

} // namespace Qt
} // namespace CGAL

Q_DECLARE_OPERATORS_FOR_FLAGS(CGAL::Qt::DemosMainWindow::Options)

#endif // CGAL_QT_DEMOS_MAIN_WINDOW_H
