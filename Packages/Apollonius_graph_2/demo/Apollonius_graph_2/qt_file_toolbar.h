#ifndef QT_FILE_TOOLBAR_H
#define QT_FILE_TOOLBAR_H

// include files for QT
#include <qapp.h>
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
#include <qmsgbox.h>
#include <qfiledialog.h>
#include <qprinter.h>
#include <qpainter.h>

#include "filesave.xpm"
#include "fileopen.xpm"
#include "filenew.xpm"
#include "fileprint.xpm"


class File_toolbar : public QObject
{
  Q_OBJECT
private:
  QToolBar     *fileToolbar;
  QMainWindow  *window;

  QAction *fileNew;
  QAction *fileOpen;
  QAction *fileSave;
  QAction *fileSaveAs;
  QAction *fileClose;
  QAction *filePrint;
public:
  /** construtor */
  File_toolbar(QMainWindow *mw)
    {
      window = mw;
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

      delete fileToolbar;
    }


  inline QToolBar* toolbar() { return fileToolbar; }

signals:
  void fileToRead(const QString& f);
  void clearAll();
  void printScreen();

private:
  void initActions()
    {
      QPixmap openIcon, saveIcon, newIcon, printIcon;
      newIcon = QPixmap(filenew);
      openIcon = QPixmap(fileopen);
      saveIcon = QPixmap(filesave);
      printIcon = QPixmap(fileprint);

      fileNew = new QAction(tr("New File"), newIcon, tr("&New"),
			    0, this);
      fileNew->setStatusTip(tr("Creates a new document"));
      fileNew->setWhatsThis(tr("New File\n\nCreates a new document"));
      connect(fileNew, SIGNAL(activated()), this, SLOT(slotFileNew()));

      fileOpen = new QAction(tr("Open File"), openIcon,
			     tr("&Open..."), 0, this);
      fileOpen->setStatusTip(tr("Opens an existing document"));
      fileOpen->setWhatsThis(tr("Open File\n\nOpens an existing document"));
      connect(fileOpen, SIGNAL(activated()), this, SLOT(slotFileOpen()));

      fileSave = new QAction(tr("Save File"), saveIcon, tr("&Save"),
			     0, this);
      fileSave->setStatusTip(tr("Saves the actual document"));
      fileSave->setWhatsThis(tr("Save File.\n\nSaves the actual document"));
      connect(fileSave, SIGNAL(activated()), this, SLOT(slotFileSave()));

      fileSaveAs = new QAction(tr("Save File As"), saveIcon, tr("Save &as..."),
			       0, this);
      fileSaveAs->setStatusTip(tr("Saves the actual document under a new filename"));
      QString _save_as =
	tr("Save As\n\nSaves the actual document under a new filename");
      fileSaveAs->setWhatsThis(_save_as);
      connect(fileSaveAs, SIGNAL(activated()), this, SLOT(slotFileSaveAs()));
#if 0
      fileSaveAs = new QAction(tr("Save File As"), tr("Save &as..."), 0, this);
      fileSaveAs->setStatusTip(tr("Saves the actual document under a new filename"));
      fileSaveAs->setWhatsThis(tr("Save As\n\nSaves the actual document under a new filename"));
      connect(fileSaveAs, SIGNAL(activated()), this, SLOT(slotFileSave()));
#endif

      fileClose = new QAction(tr("Close File"), tr("&Close"),
			      0, this);
      fileClose->setStatusTip(tr("Closes the actual document"));
      fileClose->setWhatsThis(tr("Close File\n\nCloses the actual document"));
      connect(fileClose, SIGNAL(activated()), this, SLOT(slotFileClose()));

      filePrint = new QAction(tr("Print File"), printIcon, tr("&Print"),
			      0, this);
      filePrint->setStatusTip(tr("Prints out the actual document"));
      filePrint->setWhatsThis(tr("Print File\n\nPrints out the actual document"));
      connect(filePrint, SIGNAL(activated()), this, SLOT(slotFilePrint()));
    }

  void initToolBar()
    {
      fileToolbar = new QToolBar(window, "file operations");
      fileNew->addTo(fileToolbar);
      fileOpen->addTo(fileToolbar);
      fileSaveAs->addTo(fileToolbar);
      filePrint->addTo(fileToolbar);
      fileToolbar->addSeparator();
      QWhatsThis::whatsThisButton(fileToolbar);
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
	QFileDialog::getSaveFileName(tr("data.out"), QString::null,
				     window, "Save file As...");
						      
      if ( !fileName.isNull() ) {
	std::cout << fileName << std::endl;
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

