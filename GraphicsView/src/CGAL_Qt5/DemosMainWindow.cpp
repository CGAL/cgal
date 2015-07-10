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

#include <CGAL/Qt/GraphicsViewNavigation.h>
#include <QApplication>
#include <QLabel>
#include <QFile>
#include <QFileDialog>
#include <QFileInfo>
#include <QMenu>
#include <QMenuBar>
#include <QAction>
#include <QMessageBox>
#include <QStatusBar>
#include <QGraphicsView>
#include <QGLWidget>
#include <QTextStream>
#include <QSettings>
#include <QUrl>
#include <QDesktopWidget>
#include <QRegExp>
#include <QSvgGenerator>
#include <QtCore>
#include <QtOpenGL>


#include <CGAL/config.h> // needed to get CGAL_VERSION_STR
#include <CGAL/Qt/DemosMainWindow.h>
#include <iostream>

namespace CGAL {
namespace Qt {

DemosMainWindow::DemosMainWindow(QWidget * parent, ::Qt::WindowFlags flags)
  : QMainWindow(parent, flags),
    maxNumRecentFiles(10),
    recentFileActs(maxNumRecentFiles)
{
  // prepare the QLabel xycoord for inclusion in the statusBar()
  xycoord = new QLabel(" -0.00000 , -0.00000 ", this);
  xycoord->setAlignment(::Qt::AlignHCenter);
  xycoord->setMinimumSize(xycoord->sizeHint());
  xycoord->clear();

  actionUse_OpenGL = new QAction(this);
  actionUse_OpenGL->setObjectName("actionUse_OpenGL");
  actionUse_OpenGL->setCheckable(true);
  actionUse_OpenGL->setText(tr("Use &OpenGL"));
  actionUse_OpenGL->setStatusTip(tr("Make Qt use OpenGL to display the graphical items, instead of its native painting system."));
  actionUse_OpenGL->setShortcut(tr("Ctrl+G"));

  actionUse_Antialiasing = new QAction(this);
  actionUse_Antialiasing->setObjectName("actionUse_Antialiasing");
  actionUse_Antialiasing->setCheckable(true);
  actionUse_Antialiasing->setText(tr("Use &anti-aliasing"));
  actionUse_Antialiasing->setStatusTip(tr("Make Qt use anti-aliasing when displaying the graphical items."));
  actionUse_Antialiasing->setShortcut(tr("Ctrl+A"));

  actionAboutCGAL = new QAction(this);
  actionAboutCGAL->setObjectName("actionAboutCGAL");
  actionAboutCGAL->setText(tr("About &CGAL..."));

  actionAbout = new QAction(this);
  actionAbout->setObjectName("actionAbout");
  actionAbout->setText(tr("&About..."));

  setAcceptDrops(true);
}

DemosMainWindow::~DemosMainWindow()
{
}

void 
DemosMainWindow::dragEnterEvent(QDragEnterEvent *event)
{
  if (event->mimeData()->hasFormat("text/uri-list"))
    event->acceptProposedAction();
}

void 
DemosMainWindow::dropEvent(QDropEvent *event)
{
  Q_FOREACH(QUrl url, event->mimeData()->urls()) {
    QString filename = url.toLocalFile();
    this->open(filename);
  }
  event->acceptProposedAction();
}

void
DemosMainWindow::addNavigation(QGraphicsView* graphicsView)
{
  navigation = new CGAL::Qt::GraphicsViewNavigation();
  graphicsView->viewport()->installEventFilter(navigation);
  graphicsView->installEventFilter(navigation);
  QObject::connect(navigation, SIGNAL(mouseCoordinates(QString)),
		   xycoord, SLOT(setText(QString)));
  view = graphicsView;
}

void
DemosMainWindow::setupStatusBar()
{
  this->statusBar()->addWidget(new QLabel(this), 1);
  this->statusBar()->addWidget(xycoord, 0);
}

void
DemosMainWindow::setupOptionsMenu(QMenu* menuOptions)
{
  // search for the Options menu
  if(!menuOptions) {
    menuOptions = getMenu("menuOptions", tr("&Options"));
  }

  // if not found, then create it
  if(!menuOptions) {
    menuOptions = new QMenu(this->menuBar());
    menuOptions->setTitle(tr("&Options"));
    this->menuBar()->addAction(menuOptions->menuAction());
    menuOptions->setObjectName("menuOptions");
  }

  if(!menuOptions->isEmpty()) {
    menuOptions->addSeparator();
  }
  menuOptions->addAction(actionUse_OpenGL);
  menuOptions->addAction(actionUse_Antialiasing);
  connect(actionUse_Antialiasing, SIGNAL(toggled(bool)),
          this, SLOT(setUseAntialiasing(bool)));
  connect(actionUse_OpenGL, SIGNAL(toggled(bool)),
          this, SLOT(setUseOpenGL(bool)));
  actionUse_Antialiasing->setChecked(true);
}

void
DemosMainWindow::setupExportSVG(QAction* action, QGraphicsView* view)
{
  this->view = view;
  connect(action, SIGNAL(triggered(bool)),
	  this, SLOT(exportSVG()));
}

void DemosMainWindow::exportSVG()
{
  QString fileName = QFileDialog::getSaveFileName(this,
						  tr("Export to SVG"),
						  ".",
						  tr("SVG (*.svg)\n"));

  QSvgGenerator svg;
  svg.setFileName(fileName);

  svg.setSize(this->view->size());
  svg.setViewBox(this->view->sceneRect());
  svg.setTitle(tr("%1 drawing").arg(qApp->applicationName()));
  svg.setDescription(tr("Generated using %1").arg(qApp->applicationName()));

  QPainter painter;
  painter.begin(&svg);
  this->view->render(&painter);
  painter.end();
}

void
DemosMainWindow::setUseAntialiasing(bool checked)
{
  view->setRenderHint(QPainter::Antialiasing, checked);
  view->setRenderHint(QPainter::HighQualityAntialiasing, checked);

  statusBar()->showMessage(tr("Antialiasing %1activated").arg(checked?"":"de-"),
                           1000);
}

void
DemosMainWindow::setUseOpenGL(bool checked)
{ 
  if(checked) {
    QGLWidget* new_viewport = new QGLWidget;

    // Setup the format to allow antialiasing with OpenGL:
    // one need to activate the SampleBuffers, if the graphic driver allows
    // this.
    QGLFormat glformat = new_viewport->format();
    glformat.setOption(QGL::SampleBuffers);
    new_viewport->setFormat(glformat);

    view->setViewport(new_viewport);
  }
  else {
    view->setViewport(new QWidget);
  }
  statusBar()->showMessage(tr("OpenGL %1activated").arg(checked?"":"de-"),
                           1000);
  view->viewport()->installEventFilter(navigation);
  view->setFocus();
}

QMenu* 
DemosMainWindow::getMenu(QString objectName, QString title)
{
  QMenu* menu = NULL;

  QString title2 = title;
  title2.remove('&');
  // search if a menu has objectName()==objectName
  menu = this->findChild<QMenu*>(objectName);

  // then search if a menu has title()==title
  if(menu) {
    return menu;
  } else {
    Q_FOREACH(menu, this->findChildren<QMenu*>()) {
      if(menu->title() == title ||
         menu->title() == title2) {
        return menu;
      }
    }
  }
  return NULL;
}

void 
DemosMainWindow::popupAboutBox(QString title, QString html_resource_name)
{
  QFile about_CGAL(html_resource_name);
  about_CGAL.open(QIODevice::ReadOnly);
  QString about_CGAL_txt = QTextStream(&about_CGAL).readAll();
#ifdef CGAL_VERSION_STR
  QString cgal_version(CGAL_VERSION_STR);
#  ifdef CGAL_FAKE_PUBLIC_RELEASE
  cgal_version.replace(QRegExp("-Ic?.*"), "");
#  endif
  about_CGAL_txt.replace("<!--CGAL_VERSION-->",
                         QString(" (version %1)")
                         .arg(cgal_version));
#endif
  QMessageBox mb(QMessageBox::NoIcon,
                 title,
                 about_CGAL_txt,
                 QMessageBox::Ok,
                 this);

  QLabel* mb_label = mb.findChild<QLabel*>("qt_msgbox_label");
  if(mb_label) {
    mb_label->setTextInteractionFlags(mb_label->textInteractionFlags() | 
                                      ::Qt::LinksAccessibleByMouse | 
                                      ::Qt::LinksAccessibleByKeyboard);
  }
  else {
    std::cerr << "Cannot find child \"qt_msgbox_label\" in QMessageBox\n"
              << "  with Qt version " << QT_VERSION_STR << "!\n";
  }
  mb.exec();
}

QMenu* DemosMainWindow::getHelpMenu()
{
  QMenu* menuHelp = getMenu("menuHelp", tr("&Help"));
  if(!menuHelp) {
    menuHelp = new QMenu(this->menuBar());
    menuHelp->setTitle(tr("&Help"));
    this->menuBar()->addAction(menuHelp->menuAction());
    menuHelp->setObjectName("menuHelp");
  }
  return menuHelp;
}

void 
DemosMainWindow::addAboutCGAL(QMenu* menuHelp)
{
  if(!menuHelp) {
    menuHelp = getHelpMenu();
  }
  menuHelp->addAction(actionAboutCGAL);

  connect(actionAboutCGAL, SIGNAL(triggered()),
          this, SLOT(popupAboutCGAL()));
}

void 
DemosMainWindow::addAboutDemo(QString htmlResourceName, QMenu* menuHelp)
{
  if(!menuHelp) {
    menuHelp = getHelpMenu();
  }
  menuHelp->addAction(actionAbout);
  aboutHtmlResource = htmlResourceName;

  connect(actionAbout, SIGNAL(triggered()),
          this, SLOT(popupAboutDemo()));
}

void
DemosMainWindow::popupAboutCGAL()
{
  popupAboutBox(tr("About CGAL..."),
                ":/cgal/help/about_CGAL.html");
}

void
DemosMainWindow::popupAboutDemo()
{
  popupAboutBox(tr("About the demo..."),
                aboutHtmlResource);
}

void
DemosMainWindow::setMaxNumberOfRecentFiles(const unsigned int i)
{
  maxNumRecentFiles = i;
  recentFileActs.resize(maxNumRecentFiles);
}

unsigned int 
DemosMainWindow::maxNumberOfRecentFiles() const
{
  return maxNumRecentFiles;
}

void 
DemosMainWindow::openRecentFile_aux()
{
  QAction *action = qobject_cast<QAction *>(sender());
  if (action)
    emit openRecentFile(action->data().toString());
}

void 
DemosMainWindow::addToRecentFiles(QString fileName)
{
  QSettings settings;
  QStringList files = settings.value("recentFileList").toStringList();
  files.removeAll(fileName);
  files.prepend(fileName);
  while (files.size() > (int)maxNumberOfRecentFiles())
    files.removeLast();

  settings.setValue("recentFileList", files);

  updateRecentFileActions();
}

void
DemosMainWindow::addRecentFiles(QMenu* menu, QAction* insertBeforeAction)
{
  if(!insertBeforeAction) {
    recentFilesSeparator = menu->addSeparator();
  }

  for (unsigned int i = 0; i < maxNumberOfRecentFiles(); ++i) {
    recentFileActs[i] = new QAction(this);
    recentFileActs[i]->setVisible(false);
    connect(recentFileActs[i], SIGNAL(triggered()),
            this, SLOT(openRecentFile_aux()));
    if(insertBeforeAction)
      menu->insertAction(insertBeforeAction, recentFileActs[i]);
    else
      menu->addAction(recentFileActs[i]);
  }

  if(insertBeforeAction) {
    recentFilesSeparator = menu->insertSeparator(insertBeforeAction);
  }

  recentFilesSeparator->setVisible(false);

  updateRecentFileActions();
}

void 
DemosMainWindow::updateRecentFileActions()
{
  QSettings settings;
  QStringList files = settings.value("recentFileList").toStringList();

  int numRecentFiles = qMin(files.size(), (int)this->maxNumberOfRecentFiles());
  
  for (int i = 0; i < numRecentFiles; ++i) {
    QString strippedName = QFileInfo(files[i]).fileName();
    QString text = tr("&%1 %2").arg(i).arg(strippedName);
    recentFileActs[i]->setText(text);
    recentFileActs[i]->setData(files[i]);
    recentFileActs[i]->setVisible(true);
  }
  for (unsigned int j = numRecentFiles; j < maxNumberOfRecentFiles(); ++j)
    recentFileActs[j]->setVisible(false);
  
  recentFilesSeparator->setVisible(numRecentFiles > 0);
}

void DemosMainWindow::writeState(QString groupname)
{
  QSettings settings;

  settings.beginGroup(groupname);
  settings.setValue("size", size());
  settings.setValue("pos", pos());
  settings.setValue("state", saveState());
  settings.endGroup();
}

void DemosMainWindow::readState(QString groupname, Options /*what_to_save*/)
{
  QSettings settings;
  
  settings.beginGroup(groupname);
  resize(settings.value("size", this->size()).toSize());

  QDesktopWidget* desktop = qApp->desktop();
  QPoint pos = settings.value("pos", this->pos()).toPoint();
  if(desktop->availableGeometry(pos).contains(pos)) {
    move(pos);
  }
  QByteArray mainWindowState = settings.value("state").toByteArray();
  if(!mainWindowState.isNull()) {
    this->restoreState(mainWindowState);
  }
  settings.endGroup();
}


} // namespace Qt
} // namespace CGAL
