// Copyright (c) 2008  GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Andreas Fabri <Andreas.Fabri@geometryfactory.com>
//                 Laurent Rineau <Laurent.Rineau@geometryfactory.com>

#ifndef CGAL_QT_DEMOS_MAIN_WINDOW_H
#define CGAL_QT_DEMOS_MAIN_WINDOW_H

#include <CGAL/license/GraphicsView.h>


#include <iostream>
#include <QVector>
#include <QMainWindow>
#include <QDragEnterEvent>
#include <QDropEvent>

#include <CGAL/auto_link/Qt.h>
#include <CGAL/export/Qt.h>
#include <CGAL/Qt/resources.h>
#include <QSettings>

// forward declaration
class QLabel;
class QGraphicsView;
class QAction;
class QMenu;

namespace CGAL {
namespace Qt {

// forward declaration
class GraphicsViewNavigation;

class CGAL_QT_EXPORT DemosMainWindow : public QMainWindow
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

public Q_SLOTS:
  void setMaxNumberOfRecentFiles(const unsigned int);

private:
  QMenu* getMenu(QString objectName, QString title);
  void popupAboutBox(QString title, QString html_resource_name);
  QMenu* getHelpMenu();

protected:
  DemosMainWindow (QWidget * parent = nullptr, ::Qt::WindowFlags flags = ::Qt::WindowType(0) );
  ~DemosMainWindow();
  void setupStatusBar();
  void addNavigation(QGraphicsView*);
  void setupOptionsMenu(QMenu* menu  = nullptr);
  void addAboutCGAL(QMenu* menu  = nullptr);
  void addAboutDemo(QString htmlResourceName, QMenu* menu  = nullptr);
#if QT_SVG_LIB
  void setupExportSVG(QAction*, QGraphicsView*);
#endif

  void addRecentFiles(QMenu* menu, QAction* insertBefore = nullptr);

  void writeState(QString groupname = "MainWindow");
  void readState(QString groupname = "MainWindow",
                 Options what_to_save = Options(Size|State));

protected Q_SLOTS:
  void setUseAntialiasing(bool checked);
  void popupAboutCGAL();
  void popupAboutDemo();

#if QT_SVG_LIB
  void exportSVG();
#endif
  void openRecentFile_aux();
  void addToRecentFiles(QString fileName);
  void updateRecentFileActions();

Q_SIGNALS:
  void openRecentFile(QString filename);

protected:
  QGraphicsView* view;
  GraphicsViewNavigation* navigation;
  QLabel* xycoord ;

  QAction *actionUse_Antialiasing;
  QAction *actionAbout;
  QAction *actionAboutCGAL;

  QString aboutHtmlResource;

  QAction* recentFilesSeparator;
  unsigned int maxNumRecentFiles;
  QVector<QAction*> recentFileActs;
  QSettings settings;
}; // end class DemosMainWindow

} // namespace Qt
} // namespace CGAL

Q_DECLARE_OPERATORS_FOR_FLAGS(CGAL::Qt::DemosMainWindow::Options)

#ifdef CGAL_HEADER_ONLY
#include <CGAL/Qt/DemosMainWindow_impl.h>
#endif // CGAL_HEADER_ONLY

#endif // CGAL_QT_DEMOS_MAIN_WINDOW_H

