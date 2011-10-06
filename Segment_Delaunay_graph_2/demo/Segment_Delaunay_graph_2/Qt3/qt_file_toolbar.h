// Copyright (c) 2003,2004,2005  INRIA Sophia-Antipolis (France).
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

#ifndef QT_FILE_TOOLBAR_H
#define QT_FILE_TOOLBAR_H

// include files for QT
#include <qapplication.h>
#include <qmainwindow.h>
#include <qaction.h>
#include <qmenubar.h>
#include <qpopupmenu.h>
#include <qtoolbar.h>
#include <qtoolbutton.h>
#include <qstatusbar.h>
#include <qwhatsthis.h>
#include <qstring.h>
#include <qpixmap.h>
#include <qmessagebox.h>
#include <qfiledialog.h>
#include <qprinter.h>
#include <qpainter.h>

#include "filesave.xpm"
#include "fileopen.xpm"
#include "filenew.xpm"
#include "fileprint.xpm"


class File_toolbar : public QToolBar
{
  Q_OBJECT
private:
  QMainWindow  *window;

  QAction *fileNew;
  QAction *fileOpen;
  QAction *fileSave;
  QAction *fileSaveAs;
  QAction *fileClose;
  QAction *filePrint;
public:
  /** construtor */
  File_toolbar(const QString& label, QMainWindow* mainWindow,
	       QWidget* parent, bool newLine = FALSE,
	       const char* name = 0, WFlags f = 0 )
    : QToolBar(label, mainWindow, parent, newLine, name, f),
      window(mainWindow)
  {
    initActions();
    initToolBar();
  }

  ~File_toolbar()
    {
      delete fileNew;
      delete fileOpen;
      delete fileSave;
      delete fileSaveAs;
      delete fileClose;
      delete filePrint;
    }


  inline QToolBar* toolbar() { return this; }

signals:
  void fileToRead(const QString& f);
  void fileToWrite(const QString& f);
  void printScreen();
  void clearAll();

private:
  void initActions()
    {
      QString str;
      QPixmap openIcon, saveIcon, newIcon, printIcon;
      newIcon = QPixmap(filenew);
      openIcon = QPixmap(fileopen);
      saveIcon = QPixmap(filesave);
      printIcon = QPixmap(fileprint);

      fileNew = new QAction(tr("Clear diagram"), newIcon, tr("&New"),
			    0, this);
      fileNew->setStatusTip(tr("Clears the current Delaunay graph"));
      fileNew->setWhatsThis
	(tr("Clear diagram\n\nClears the current Delaunay graph"));
      connect(fileNew, SIGNAL(activated()), this, SLOT(slotFileNew()));

      fileOpen = new QAction(tr("Open file"), openIcon,
			     tr("&Open..."), 0, this);
      fileOpen->setStatusTip(tr("Opens a file with input data"));
      fileOpen->setWhatsThis(tr("Open file\n\nOpens a file with input data"));
      connect(fileOpen, SIGNAL(activated()), this, SLOT(slotFileOpen()));

      fileSave = new QAction(tr("Save to file"), saveIcon, tr("&Save"),
			     0, this);
      fileSave->setStatusTip(tr("Saves the input sites to a file"));
      str = "Save to file\n\nSaves the input sites to a file";
      fileSave->setWhatsThis(tr(str));
      connect(fileSave, SIGNAL(activated()), this, SLOT(slotFileSave()));

      fileSaveAs = new QAction(tr("Save to new file"), saveIcon,
			       tr("Save &as..."), 0, this);
      str = "Saves the output sites to a new file";
      fileSaveAs->setStatusTip(tr(str));
      str = "Save as\n\nSaves the output sites to a new file";
      fileSaveAs->setWhatsThis(tr(str));
      connect(fileSaveAs, SIGNAL(activated()), this, SLOT(slotFileSaveAs()));

      fileClose = new QAction(tr("Close file"), tr("&Close"),
			      0, this);
      fileClose->setStatusTip(tr("Closes the actual document"));
      fileClose->setWhatsThis(tr("Close file\n\nCloses the actual document"));
      connect(fileClose, SIGNAL(activated()), this, SLOT(slotFileClose()));

      filePrint = new QAction(tr("Print canvas"), printIcon, tr("&Print"),
			      0, this);
      filePrint->setStatusTip(tr("Prints out the current canvas"));
      str = "Print canvas\n\nPrints out the current canvas";
      filePrint->setWhatsThis(tr(str));
      connect(filePrint, SIGNAL(activated()), this, SLOT(slotFilePrint()));
    }

  void initToolBar()
    {
      fileNew->addTo(this);
      fileOpen->addTo(this);
      fileSaveAs->addTo(this);
      filePrint->addTo(this);
      this->addSeparator();
      QWhatsThis::whatsThisButton(this);
    }

  /** setup the statusbar */

public slots:
  /** generate a new document in the actual view */
  void slotFileNew()
    {
      emit clearAll();
    }
  /** open a document */
  void slotFileOpen()
    {
      QString fileName =
	QFileDialog::getOpenFileName(QString::null, QString::null,
				     window, "Open file...");

      if ( !fileName.isNull() ) {
	emit fileToRead(fileName);
      }
    }
  /** save a document */
  void slotFileSave()
    {

    }
  /** save a document under a different filename*/
  void slotFileSaveAs()
    {
      QString fileName =
	QFileDialog::getSaveFileName(tr("data.cin"), QString::null,
				     window, "Save data as...");

      if ( !fileName.isNull() ) {
	emit fileToWrite(fileName);
      }
    }
  /** close the actual file */
  void slotFileClose() {}
  /** print the actual file */
  void slotFilePrint() {
    emit printScreen();
  }
  /** put the marked text/object into the clipboard and remove
   * it from the document */
};

#endif // QT_FILE_TOOLBAR_H
