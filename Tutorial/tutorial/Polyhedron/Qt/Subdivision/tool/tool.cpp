/***************************************************************************

    begin                : jan 02
    copyright            : (C) 2002 by Pierre Alliez
    email                : pierre.alliez@sophia.inria.fr
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

// Qt includes
#include <qvbox.h>
#include <qaccel.h>

// application specific includes
#include "toolview.h"
#include "tooldoc.h"
#include "tool.h"

#include "icons/cgal.xpm"
#include "icons/filenew.xpm"
#include "icons/fileopen.xpm"
#include "icons/filesave.xpm"
#include "icons/fileprint.xpm"
#include "icons/param.xpm"
#include "icons/directions.xpm"

#include "icons/points.xpm"
#include "icons/wireframe.xpm"
#include "icons/fill.xpm"
#include "icons/superimposed_edges.xpm"
#include "icons/superimposed_vertices.xpm"
#include "icons/light.xpm"
#include "icons/antialiasing.xpm"
#include "icons/smooth.xpm"

#include "icons/quadtriangle.xpm"
#include "icons/sqrt3.xpm"
#include "icons/sqrt3_twice.xpm"

ToolApp::ToolApp()
{
  // setCaption("Subdivision tool", "1.5");

  printer = new QPrinter;
  untitledCount=0;
  pDocList = new QList<ToolDoc>();
  pDocList->setAutoDelete(true);

  // modeless dialogs
  m_pDlgParamDirection = NULL;

  ///////////////////////////////////////////////////////////////////
  // call inits to invoke all other construction parts
  initView();
  initActions();
  initMenuBar();
  initToolBar();
  initStatusBar();
  resize( 450, 400 );

  viewToolBar->setOn(true);
  viewRenderingToolbar->setOn(false);
  renderingToolbar->hide();
  viewStatusBar->setOn(true);
  viewRenderingToolbar->setEnabled(false);
  renderActionGroup->setEnabled(false);
  subdivisionActionGroup->setEnabled(false);
}

ToolApp::~ToolApp()
{
  delete printer;
}

void ToolApp::initActions()
{
  // pixmap
  QPixmap newIcon = QPixmap(filenew);
  QPixmap openIcon = QPixmap(fileopen);
  QPixmap saveIcon = QPixmap(filesave);
  QPixmap printIcon = QPixmap(fileprint);
  QPixmap paramIcon = QPixmap(param);
  QPixmap directionIcon = QPixmap(directions);

  QPixmap pointsIcon = QPixmap(points_xpm);
  QPixmap wireframeIcon = QPixmap(wireframe_xpm);
  QPixmap fillIcon = QPixmap(fill_xpm);
  QPixmap superimposededgesIcon = QPixmap(superimposed_edges_xpm);
  QPixmap superimposedverticesIcon = QPixmap(superimposed_vertices_xpm);
  QPixmap lightIcon = QPixmap(light_xpm);
  QPixmap antialiasingIcon = QPixmap(antialiasing_xpm);
  QPixmap smoothIcon = QPixmap(smooth_xpm);

  // subdivision pixmaps
  QPixmap quadtriangleIcon = QPixmap(quadtriangle_xpm);
  QPixmap sqrt3Icon = QPixmap(sqrt3_xpm);
  QPixmap sqrt3_twiceIcon = QPixmap(sqrt3_twice_xpm);

  fileOpen = new QAction(tr("Open File"), openIcon, tr("&Open..."), 0, this);
  fileOpen->setStatusTip(tr("Opens an existing document"));
  fileOpen->setWhatsThis(tr("Open File\n\nOpens an existing document"));
  connect(fileOpen, SIGNAL(activated()), this, SLOT(slotFileOpen()));

  fileSave = new QAction(tr("Save File"), saveIcon, tr("&Save"),
 QAccel::stringToKey(tr("Ctrl+S")), this);
  fileSave->setStatusTip(tr("Saves the actual document"));
  fileSave->setWhatsThis(tr("Save File.\n\nSaves the actual document"));
  connect(fileSave, SIGNAL(activated()), this, SLOT(slotFileSave()));

  fileSaveAs = new QAction(tr("Save File As"), tr("Save &as..."), 0, this);
  fileSaveAs->setStatusTip(tr("Saves the actual document under a new filename"));
  fileSaveAs->setWhatsThis(tr("Save As\n\nSaves the actual document under a new filename"));
  connect(fileSaveAs, SIGNAL(activated()), this, SLOT(slotFileSave()));

  fileClose = new QAction(tr("Close File"), tr("&Close"),
 QAccel::stringToKey(tr("Ctrl+W")), this);
  fileClose->setStatusTip(tr("Closes the actual document"));
  fileClose->setWhatsThis(tr("Close File\n\nCloses the actual document"));
  connect(fileClose, SIGNAL(activated()), this, SLOT(slotFileClose()));

  filePrint = new QAction(tr("Print File"), printIcon, tr("&Print"),
 QAccel::stringToKey(tr("Ctrl+P")), this);
  filePrint->setStatusTip(tr("Prints out the actual document"));
  filePrint->setWhatsThis(tr("Print File\n\nPrints out the actual document"));
  connect(filePrint, SIGNAL(activated()), this, SLOT(slotFilePrint()));

  fileQuit = new QAction(tr("Exit"), tr("E&xit"),
 QAccel::stringToKey(tr("Ctrl+Q")), this);
  fileQuit->setStatusTip(tr("Quits the application"));
  fileQuit->setWhatsThis(tr("Exit\n\nQuits the application"));
  connect(fileQuit, SIGNAL(activated()), this, SLOT(slotFileQuit()));

  editCut = new QAction(tr("Cut"), tr("Cu&t"),
 QAccel::stringToKey(tr("Ctrl+X")), this);
  editCut->setStatusTip(tr("Cuts the selected section and puts it to the clipboard"));
  editCut->setWhatsThis(tr("Cut\n\nCuts the selected section and puts it to the clipboard"));
  connect(editCut, SIGNAL(activated()), this, SLOT(slotEditCut()));

  editCopy = new QAction(tr("Copy"), tr("&Copy"),
 QAccel::stringToKey(tr("Ctrl+C")), this);
  editCopy->setStatusTip(tr("Copies the selected section to the clipboard"));
  editCopy->setWhatsThis(tr("Copy\n\nCopies the selected section to the clipboard"));
  connect(editCopy, SIGNAL(activated()), this, SLOT(slotEditCopy()));

  editUndo = new QAction(tr("Undo"), tr("&Undo"),
 QAccel::stringToKey(tr("Ctrl+Z")), this);
  editUndo->setStatusTip(tr("Reverts the last editing action"));
  editUndo->setWhatsThis(tr("Undo\n\nReverts the last editing action"));
  connect(editUndo, SIGNAL(activated()), this, SLOT(slotEditUndo()));

  editPaste = new QAction(tr("Paste"), tr("&Paste"),
 QAccel::stringToKey(tr("Ctrl+V")), this);
  editPaste->setStatusTip(tr("Pastes the clipboard contents to actual position"));
  editPaste->setWhatsThis(tr("Paste\n\nPastes the clipboard contents to actual position"));
  connect(editPaste, SIGNAL(activated()), this, SLOT(slotEditPaste()));

  viewToolBar = new QAction(tr("Toolbar"), tr("Tool&bar"), 0, this, 0, true);
  viewToolBar->setStatusTip(tr("Enables/disables the toolbar"));
  viewToolBar->setWhatsThis(tr("Toolbar\n\nEnables/disables the toolbar"));
  connect(viewToolBar, SIGNAL(toggled(bool)), this, SLOT(slotViewToolBar(bool)));

  viewRenderingToolbar = new QAction(tr("Rendering toolbar"), tr("Rendering tool&bar"), 0, this, 0, true);
  viewRenderingToolbar->setStatusTip(tr("Enables/disables the rendering toolbar"));
  viewRenderingToolbar->setWhatsThis(tr("Toolbar\n\nEnables/disables the rendering toolbar"));
  connect(viewRenderingToolbar, SIGNAL(toggled(bool)), this, SLOT(slotViewRenderingToolBar(bool)));

  viewStatusBar = new QAction(tr("Statusbar"), tr("&Statusbar"), 0, this, 0,
 true);
  viewStatusBar->setStatusTip(tr("Enables/disables the statusbar"));
  viewStatusBar->setWhatsThis(tr("Statusbar\n\nEnables/disables the statusbar"));
  connect(viewStatusBar, SIGNAL(toggled(bool)), this,
 SLOT(slotViewStatusBar(bool)));

  // ** copy/paste trackball **
  viewCopyViewpoint = new QAction("Copy viewpoint","&Copy viewpoint",0,this);
  viewCopyViewpoint->setStatusTip("Copy the current viewpoint");
  connect(viewCopyViewpoint, SIGNAL(activated()), this, SLOT(slotViewCopyViewpoint()));

  viewPasteViewpoint = new QAction("Paste viewpoint","&Paste viewpoint",0,this);
  viewPasteViewpoint->setStatusTip("Copy the current viewpoint");
  connect(viewPasteViewpoint, SIGNAL(activated()), this, SLOT(slotViewPasteViewpoint()));

  windowNewWindow = new QAction(tr("New Window"), tr("&New Window"), 0, this);
  windowNewWindow->setStatusTip(tr("Opens a new view for the current document"));
  windowNewWindow->setWhatsThis(tr("New Window\n\nOpens a new view for the current document"));
  connect(windowNewWindow, SIGNAL(activated()), this, SLOT(slotWindowNewWindow()));

  windowCascade = new QAction(tr("Cascade"), tr("&Cascade"), 0, this);
  windowCascade->setStatusTip(tr("Cascades all windows"));
  windowCascade->setWhatsThis(tr("Cascade\n\nCascades all windows"));
  connect(windowCascade, SIGNAL(activated()), pWorkspace, SLOT(cascade()));

  windowTile = new QAction(tr("Tile"), tr("&Tile"), 0, this);
  windowTile->setStatusTip(tr("Tiles all windows"));
  windowTile->setWhatsThis(tr("Tile\n\nTiles all windows"));
  connect(windowTile, SIGNAL(activated()), pWorkspace, SLOT(tile()));


  windowAction = new QActionGroup(this, 0, false);
  windowAction->insert(windowNewWindow);
  windowAction->insert(windowCascade);
  windowAction->insert(windowTile);

  helpAboutApp = new QAction(tr("About"), tr("&About..."), 0, this);
  helpAboutApp->setStatusTip(tr("About the application"));
  helpAboutApp->setWhatsThis(tr("About\n\nAbout the application"));
  connect(helpAboutApp, SIGNAL(activated()), this, SLOT(slotHelpAbout()));


  // RENDER
  renderWireframe = new QAction(tr("Wireframe"), wireframeIcon, tr("&Wireframe"),QAccel::stringToKey(tr("w")),this, 0, true);
  connect(renderWireframe,SIGNAL(activated()), this, SLOT(slotRenderWireframe()));

  renderVertex = new QAction(tr("Vertex"), pointsIcon, tr("&Vertex"),QAccel::stringToKey(tr("v")),this, 0, true);
  connect(renderVertex,SIGNAL(activated()), this, SLOT(slotRenderVertex()));

  renderFill = new QAction(tr("Fill"), fillIcon, tr("&Fill"),QAccel::stringToKey(tr("f")),this, 0, true);
  connect(renderFill,SIGNAL(activated()), this, SLOT(slotRenderFill()));

  renderEdges = new QAction(tr("Edges"), tr("&Edges"),QAccel::stringToKey(tr("e")),this, 0, true);
  connect(renderEdges,SIGNAL(activated()), this, SLOT(slotRenderEdges()));

  renderSmooth = new QAction(tr("Smooth"), smoothIcon, tr("&Smooth"),QAccel::stringToKey(tr("z")),this, 0, true);
  connect(renderSmooth,SIGNAL(activated()), this, SLOT(slotRenderSmooth()));

  renderCulling = new QAction(tr("Culling"), tr("&Culling"),QAccel::stringToKey(tr("c")),this, 0, true);
  connect(renderCulling,SIGNAL(activated()), this, SLOT(slotRenderCulling()));

  renderLighting = new QAction(tr("Lighting"), lightIcon, tr("&Lighting"),QAccel::stringToKey(tr("l")),this, 0, true);
  connect(renderLighting,SIGNAL(activated()), this, SLOT(slotRenderLighting()));

  renderAntialiasing = new QAction(tr("Antialiasing"), antialiasingIcon, tr("&Antialiasing"),QAccel::stringToKey("a"),this, 0, true);
  connect(renderAntialiasing,SIGNAL(activated()), this, SLOT(slotRenderAntialiasing()));

  renderSuperimpose = new QAction(tr("Superimpose edges"), superimposededgesIcon, tr("&Superimpose edges"),QAccel::stringToKey("s"),this, 0, true);
  connect(renderSuperimpose,SIGNAL(activated()), this, SLOT(slotRenderSuperimposing()));

  renderSuperimposeV = new QAction(tr("Superimpose vertices"), superimposedverticesIcon, tr("&Superimpose vertices"),QAccel::stringToKey("s"),this, 0, true);
  connect(renderSuperimposeV, SIGNAL(activated()), this, SLOT(slotRenderSuperimposedV()));

  renderActionGroup = new QActionGroup(this, 0, false);
  renderActionGroup->add(renderWireframe);
  renderActionGroup->add(renderVertex);
  renderActionGroup->add(renderFill);
  renderActionGroup->add(renderSmooth);
  renderActionGroup->add(renderCulling);
  renderActionGroup->add(renderLighting);
  renderActionGroup->add(renderAntialiasing);
  renderActionGroup->add(renderSuperimpose);
  renderActionGroup->add(renderSuperimposeV);

  // Subdivision schemes
  algoSubdivisionStamLoop = new QAction(tr("Quad/Triangle"), quadtriangleIcon, tr("&Quad/Triangle"), 0, this);
  connect(algoSubdivisionStamLoop, SIGNAL(activated()), this, SLOT(slotSubdivisionStamLoop()));
  
  algoSubdivisionSqrt3 = new QAction(tr("Sqrt3"),  sqrt3Icon, tr("&Sqrt3"), 0, this);
  connect(algoSubdivisionSqrt3, SIGNAL(activated()), this, SLOT(slotSubdivisionSqrt3()));

  algoSubdivisionSqrt3_twice = new QAction(tr("Sqrt3 (apply twice)"), sqrt3_twiceIcon, tr("&Sqrt3 (apply twice)"), 0, this);
  connect(algoSubdivisionSqrt3_twice, SIGNAL(activated()), this, SLOT(slotSubdivisionSqrt3_twice()));

  algoSubdivisionLoop = new QAction("Loop", "&Loop", 0, this);
  connect(algoSubdivisionLoop, SIGNAL(activated()), this, SLOT(slotSubdivisionLoop()));

  algoSubdivisionCatmullClark = new QAction("Catmull-Clark", "&Catmull-Clark", 0, this);
  connect(algoSubdivisionCatmullClark, SIGNAL(activated()), this, SLOT(slotSubdivisionCatmullClark()));

  algoSubdivisionDooSabin = new QAction("Doo-Sabin", "&Doo-Sabin", 0, this);
  connect(algoSubdivisionDooSabin, SIGNAL(activated()), this, SLOT(slotSubdivisionDooSabin()));

  subdivisionActionGroup = new QActionGroup(this, 0, false);
  subdivisionActionGroup->add(algoSubdivisionStamLoop);
  subdivisionActionGroup->add(algoSubdivisionSqrt3);
  subdivisionActionGroup->add(algoSubdivisionSqrt3_twice);
  subdivisionActionGroup->add(algoSubdivisionLoop);
  subdivisionActionGroup->add(algoSubdivisionCatmullClark);
  subdivisionActionGroup->add(algoSubdivisionDooSabin);

}

void ToolApp::initMenuBar()
{
  // MENUBAR

  // menuBar entry pFileMenu
  pFileMenu = new QPopupMenu();
  fileOpen->addTo(pFileMenu);
  fileClose->addTo(pFileMenu);
  pFileMenu->insertSeparator();
  fileSave->addTo(pFileMenu);
  fileSaveAs->addTo(pFileMenu);
  pFileMenu->insertSeparator();
  filePrint->addTo(pFileMenu);
  pFileMenu->insertSeparator();
  fileQuit->addTo(pFileMenu);

  // menuBar entry editMenu
  pEditMenu=new QPopupMenu();
  editUndo->addTo(pEditMenu);
  pEditMenu->insertSeparator();
  editCut->addTo(pEditMenu);
  editCopy->addTo(pEditMenu);
  editPaste->addTo(pEditMenu);

  // RENDER
  pRenderMenu = new QPopupMenu();
  pRenderModeMenu = new QPopupMenu();
  pRenderModeMenu->setCheckable(true);
  QActionGroup *ag1 = new QActionGroup(this, "mode menu", 1);
  ag1->add(renderWireframe);
  ag1->add(renderVertex);
  ag1->add(renderFill);
  ag1->addTo(pRenderModeMenu);
  
  pRenderSuperimposeMenu = new QPopupMenu();
  pRenderSuperimposeMenu->setCheckable(true);
  renderSuperimpose->addTo(pRenderSuperimposeMenu);
  renderSuperimposeV->addTo(pRenderSuperimposeMenu);
  renderAntialiasing->addTo(pRenderMenu);
  renderCulling->addTo(pRenderMenu);
  renderSmooth->addTo(pRenderMenu);
  renderLighting->addTo(pRenderMenu);

  ///////////////////////////////////////////////////////////////////
  // menuBar entry viewMenu
  pViewMenu=new QPopupMenu();
  pViewMenu->setCheckable(true);
  viewToolBar->addTo(pViewMenu);
  viewRenderingToolbar->addTo(pViewMenu);
  viewStatusBar->addTo(pViewMenu);
  viewCopyViewpoint->addTo(pViewMenu);
  viewPasteViewpoint->addTo(pViewMenu);

  ///////////////////////////////////////////////////////////////////
  // EDIT YOUR APPLICATION SPECIFIC MENUENTRIES HERE

  ///////////////////////////////////////////////////////////////////
  // menuBar entry windowMenu
  pWindowMenu = new QPopupMenu(this);
  pWindowMenu->setCheckable(true);
  connect(pWindowMenu, SIGNAL(aboutToShow()), this,
 SLOT(windowMenuAboutToShow()));

  ///////////////////////////////////////////////////////////////////
  // menuBar entry helpMenu
  pHelpMenu=new QPopupMenu();
  helpAboutApp->addTo(pHelpMenu);
  pHelpMenu->insertSeparator();
  pHelpMenu->insertItem(tr("What's &This"), this, SLOT(whatsThis()),
 SHIFT+Key_F1);

  menuBar()->insertItem(tr("&File"),       pFileMenu);
  menuBar()->insertItem(tr("&Edit"),       pEditMenu);

  // SUBDIVISION
  pAlgoMenu = new QPopupMenu();
  menuBar()->insertItem(tr("&Subdivision"), pAlgoMenu);
  algoSubdivisionStamLoop->addTo(pAlgoMenu);
  algoSubdivisionSqrt3->addTo(pAlgoMenu);
  algoSubdivisionSqrt3_twice->addTo(pAlgoMenu);
  algoSubdivisionDooSabin->addTo(pAlgoMenu);
  algoSubdivisionCatmullClark->addTo(pAlgoMenu);
  algoSubdivisionLoop->addTo(pAlgoMenu);

  // RENDERING OPTIONS
  menuBar()->insertItem(tr("&Render"),        pRenderMenu);
  pRenderMenu->insertItem(tr("&Mode"),        pRenderModeMenu);
  pRenderMenu->insertItem(tr("&Superimpose"), pRenderSuperimposeMenu);

  menuBar()->insertItem(tr("&View"),       pViewMenu);
  menuBar()->insertItem(tr("&Window"),     pWindowMenu);
  menuBar()->insertItem(tr("&Help"),       pHelpMenu);
}

void ToolApp::initToolBar()
{
  ///////////////////////////////////////////////////////////////////
  // TOOLBAR FILE OPERATIONS
  //////////////////////////////////////////////////////////////////
  fileToolbar = new QToolBar(this, "file operations");
  fileOpen->addTo(fileToolbar);
  fileSave->addTo(fileToolbar);
  filePrint->addTo(fileToolbar);
  fileToolbar->addSeparator();

  /////////////////////////////////////////////////////////////////
  // TOOLBAR RENDERING
  //////////////////////////////////////////////////////////////////

  renderingToolbar = new QToolBar(this, "rendering toolbar");
  QLabel *l1 = new QLabel(renderingToolbar, "BgColor");
  l1->setText("BgCol");
  bgcolor_button = new QToolButton(renderingToolbar, "bgcolor");
  bgcolor_button->setText("BG");
  bgcolor_button->setTextLabel("Change background color");
  bgcolor_button->setPaletteBackgroundColor(QColor(0, 0, 0));
  bgcolor_button->setPaletteForegroundColor(QColor(0, 0, 0));
  connect(bgcolor_button, SIGNAL(clicked()), this, SLOT(slotSelectBgColor()));
  QLabel *l2 = new QLabel(renderingToolbar, "FG");
  l2->setText("FaceCol");
  fgcolor_button = new QToolButton(renderingToolbar, "fgcolor");
  fgcolor_button->setText("FG");
  fgcolor_button->setTextLabel("Change face color");
  fgcolor_button->setPaletteBackgroundColor(QColor(255, 255, 255));
  fgcolor_button->setPaletteForegroundColor(QColor(255, 255, 255));
  connect(fgcolor_button, SIGNAL(clicked()), this, SLOT(slotSelectFgColor()));

  renderingToolbar->addSeparator();
  renderVertex->addTo(renderingToolbar);
  renderWireframe->addTo(renderingToolbar);
  renderFill->addTo(renderingToolbar);
  renderingToolbar->addSeparator();
  renderSuperimpose->addTo(renderingToolbar);
  renderSuperimposeV->addTo(renderingToolbar);
  renderingToolbar->addSeparator();
  renderSmooth->addTo(renderingToolbar);
  renderLighting->addTo(renderingToolbar);
  renderAntialiasing->addTo(renderingToolbar);

  /////////////////////////////////////////////////////////////////
  // SUBDIVISION TOOLBAR
  //////////////////////////////////////////////////////////////////

  subdivisionToolbar = new QToolBar(this,"subdivision toolbar");
  algoSubdivisionStamLoop->addTo(subdivisionToolbar);
  algoSubdivisionSqrt3->addTo(subdivisionToolbar);
  algoSubdivisionSqrt3_twice->addTo(subdivisionToolbar);
  algoSubdivisionDooSabin->addTo(subdivisionToolbar);
  //algoSubdivisionDooSabin->addTo(subdivisionToolbar);
  algoSubdivisionCatmullClark->addTo(subdivisionToolbar);
  algoSubdivisionLoop->addTo(subdivisionToolbar);
}

void ToolApp::initStatusBar()
{
  ///////////////////////////////////////////////////////////////////
  //STATUSBAR
  statusBar()->message(tr("Ready."));
}

void ToolApp::initView()
{
  // set the main widget here
  QVBox* view_back = new QVBox( this );
  view_back->setFrameStyle( QFrame::StyledPanel | QFrame::Sunken );
  pWorkspace = new QWorkspace( view_back );
  setCentralWidget(view_back);
}

void ToolApp::createClient(ToolDoc* doc)
{
  // set OpenGL options
  QGLFormat format;
  format.setAlpha(false);
  format.setStereo(false);
  format.setDepth(true);
  format.setDoubleBuffer(true);
  format.setDirectRendering(true);

  ToolView* w = new ToolView(doc, format, pWorkspace, 0, WDestructiveClose);
  w->installEventFilter(this);  
  w->resize(200, 200);
  QPixmap cgalIcon = QPixmap(cgal_xpm);
  w->setIcon(cgalIcon);
  doc->addView(w);
  connect(w, SIGNAL(invalidate_view(ToolView*)), this, SLOT(slotUpdateStates(ToolView*)));
  w->show();
  w->toggleMode(ToolView::FILL);
  w->toggleLighting();
  renderLighting->setOn(true);
  renderFill->setOn(true);
  renderingToolbar->repaint();
  w->show();
  slotUpdateStates(w);
}

void ToolApp::openDocumentFile(const char* file)
{
  statusBar()->message(tr("Opening file..."));
  ToolDoc* doc;
  // check, if document already open. If yes, set the focus to the first view
  for(doc=pDocList->first(); doc > 0; doc=pDocList->next())
  {
    if(doc->pathName()==file)
    {
      ToolView* view=doc->firstView();
      view->setFocus();
      return;
     }
  }
  doc = new ToolDoc();
  pDocList->append(doc);
  doc->newDocument();
  // Creates an untitled window if file is 0
  if(!file)
  {
    untitledCount+=1;
    QString fileName=QString(tr("Untitled%1")).arg(untitledCount);
    doc->setPathName(fileName);
    doc->setTitle(fileName);
  }
  // Open the file
  else
  {
    if(!doc->openDocument(file))
    {
      QMessageBox::critical(this, tr("Error !"),tr("Could not open document !"));
      delete doc;
      return;
    }
  }
  // create the window
  createClient(doc);  
  statusBar()->message(tr("Ready."));
}

bool ToolApp::queryExit()
{
  int exit=QMessageBox::information(this, tr("Quit..."), "Do your really want to quit?",
    QMessageBox::Ok, QMessageBox::Cancel);

  if (exit==1)
  {

  }
  else
  {

  };

  return (exit==1);
}

bool ToolApp::eventFilter(QObject* object, QEvent* event)
{
  if((event->type() == QEvent::Close)&&((ToolApp*)object!=this))
  {
    QCloseEvent* e=(QCloseEvent*)event;
    ToolView* pView=(ToolView*)object;
    ToolDoc* pDoc=pView->getDocument();
    if(pDoc->canCloseFrame(pView))
    {
      pDoc->removeView(pView);
      if(!pDoc->firstView())
        pDocList->remove(pDoc);

      e->accept();
    }
    else
      e->ignore();
  }
  return QWidget::eventFilter( object, event );    // standard event processing
}

/////////////////////////////////////////////////////////////////////
// SLOT IMPLEMENTATION
/////////////////////////////////////////////////////////////////////

/*
void ToolApp::slotFileNew()
{
  statusBar()->message(tr("Creating new file..."));
  openDocumentFile();
  statusBar()->message(tr("Ready."));
}
*/
void ToolApp::slotFileOpen()
{
  statusBar()->message(tr("Opening file..."));

  QString fileName = QFileDialog::getOpenFileName("../../../data",0,this);
  if (!fileName.isEmpty())
  {
     openDocumentFile(fileName);
  }
  renderingToolbar->show();
  viewRenderingToolbar->setEnabled(true);
  viewRenderingToolbar->setOn(true);
  renderActionGroup->setEnabled(true);
  subdivisionActionGroup->setEnabled(true);
  statusBar()->message(tr("Ready."));
}


void ToolApp::slotFileSave()
{
  statusBar()->message(tr("Saving file..."));

  ToolView* m = (ToolView*)pWorkspace->activeWindow();
  if(m)
  {
    //ToolDoc* doc = m->getDocument();
    //if(doc->title().contains(tr("Untitled")))
      slotFileSaveAs();
    //else
      //if(!doc->saveDocument(doc->pathName()))
        //QMessageBox::critical (this, tr("I/O Error !"),
	//tr("Could not save the current document !"));
  }

  statusBar()->message(tr("Ready."));
}

void ToolApp::slotFileSaveAs()
{
  statusBar()->message(tr("Saving file under new filename..."));
  QString fn = QFileDialog::getSaveFileName(0, 0, this);
  if (!fn.isEmpty())
  {
    ToolView* m = (ToolView*)pWorkspace->activeWindow();
    if( m )
    {
      ToolDoc* doc = m->getDocument();
      if(!doc->saveDocument(fn))
      {
         QMessageBox::critical (this, tr("I/O Error !"),
	 tr("Could not save the current document !"));
         return;
      }
      doc->changedViewList();
    }
  }
  statusBar()->message(tr("Ready."));
}

void ToolApp::slotFileClose()
{
  statusBar()->message(tr("Closing file..."));

  ToolView* m = (ToolView*)pWorkspace->activeWindow();
  if( m )
  {
    ToolDoc* doc=m->getDocument();
    doc->closeDocument();
  }

  statusBar()->message(tr("Ready."));
}

void ToolApp::slotFilePrint()
{
  statusBar()->message(tr("Printing..."));

  ToolView* m = (ToolView*) pWorkspace->activeWindow();
  if ( m )
    m->print( printer );

  statusBar()->message(tr("Ready."));
}

void ToolApp::slotFileQuit()
{
  statusBar()->message(tr("Exiting application..."));
  ///////////////////////////////////////////////////////////////////
  // exits the Application
//  if(doc->isModified())
//  {
//    if(queryExit())
//    {
//      qApp->quit();
//    }
//    else
//    {
//
//    };
//  }
//  else
//  {
    qApp->quit();
//  };

  statusBar()->message(tr("Ready."));
}

void ToolApp::slotEditUndo()
{
  statusBar()->message(tr("Reverting last action..."));

  ToolView* m = (ToolView*) pWorkspace->activeWindow();
  if ( m )
//   m->undo();

  statusBar()->message(tr("Ready."));
}

void ToolApp::slotEditCut()
{
  statusBar()->message(tr("Cutting selection..."));

  ToolView* m = (ToolView*) pWorkspace->activeWindow();
  if ( m )
//  m->cut();

  statusBar()->message(tr("Ready."));
}

void ToolApp::slotEditCopy()
{
  statusBar()->message(tr("Copying selection to clipboard..."));

  ToolView* m = (ToolView*) pWorkspace->activeWindow();
  if ( m )
//  m->copy();

  statusBar()->message(tr("Ready."));
}

void ToolApp::slotEditPaste()
{
  statusBar()->message(tr("Inserting clipboard contents..."));

  ToolView* m = (ToolView*) pWorkspace->activeWindow();
  if ( m )
//   m->paste();

  statusBar()->message(tr("Ready."));
}


void ToolApp::slotViewToolBar(bool toggle)
{
  statusBar()->message(tr("Toggle toolbar..."));
  ///////////////////////////////////////////////////////////////////
  // turn Toolbar on or off
   if (toggle== false)
  {
    fileToolbar->hide();
  }
  else
  {
    fileToolbar->show();
  };

 statusBar()->message(tr("Ready."));
}

void ToolApp::slotViewRenderingToolBar(bool toggle)
{
 statusBar()->message(tr("Toggle rendering toolbar..."));
  ///////////////////////////////////////////////////////////////////
  // turn Rendering Toolbar on or off
  if (toggle== false)
  {
    renderingToolbar->hide();
  }
  else
  {
    renderingToolbar->show();
  };

 statusBar()->message(tr("Ready."));
}

//**********************************************
// Copy current viewpoint
//**********************************************
void ToolApp::slotViewCopyViewpoint()
{
  statusBar()->message(tr("Copy viewpoint..."));
  ToolView* pView = (ToolView*)pWorkspace->activeWindow();
  if(pView)
    theArcball.Set(pView->theArcball);
  statusBar()->message(tr("Ready."));
}

//**********************************************
// Copy current viewpoint
//**********************************************
void ToolApp::slotViewPasteViewpoint()
{
  statusBar()->message(tr("Paste viewpoint..."));
  ToolView* pView = (ToolView*)pWorkspace->activeWindow();
  if(pView)
  {
    pView->theArcball.Set(theArcball);
    ToolDoc* pDoc = pView->getDocument();
    pDoc->updateAllViews(NULL);
  }
  statusBar()->message(tr("Ready."));
}

void ToolApp::slotViewStatusBar(bool toggle)
{
  statusBar()->message(tr("Toggle statusbar..."));
  if (toggle == false)
    statusBar()->hide();
  else
    statusBar()->show();
  statusBar()->message(tr("Ready."));
}

void ToolApp::slotWindowNewWindow()
{
  statusBar()->message(tr("Opening new document view..."));

  ToolView* m = (ToolView*) pWorkspace->activeWindow();
  if ( m ){
    ToolDoc* doc = m->getDocument();
    createClient(doc);
  }

  statusBar()->message(tr("Ready."));
}

void ToolApp::slotUpdateStates(ToolView* m){  
  if ( m ){
    if(m->m_Lighting)
      renderLighting->setOn(true);
    else
      renderLighting->setOn(false);
    if(m->m_Antialiasing)
      renderAntialiasing->setOn(true);
    else
      renderAntialiasing->setOn(false);
    if(m->m_Culling)
      renderCulling->setOn(true);
    else
      renderCulling->setOn(false);
    if(m->m_Smooth)
      renderSmooth->setOn(true);
    else
      renderSmooth->setOn(false);
    if(m->m_Superimpose_edges)
      renderSuperimpose->setOn(true);
    else
      renderSuperimpose->setOn(false);
    if(m->m_mode == ToolView::FILL)
      renderFill->setOn(true);
    else if (m->m_mode == ToolView::LINE)
      renderWireframe->setOn(true);
    else
      renderVertex->setOn(true);
    bgcolor_button->setPaletteBackgroundColor(
        QColor(m->m_clear_color[0]*255, m->m_clear_color[1]*255, m->m_clear_color[2]*255));
    bgcolor_button->setPaletteForegroundColor(
        QColor(m->m_clear_color[0]*255, m->m_clear_color[1]*255, m->m_clear_color[2]*255));

    fgcolor_button->setPaletteBackgroundColor(
        QColor(m->m_fore_color[0]*255, m->m_fore_color[1]*255, m->m_fore_color[2]*255));
    fgcolor_button->setPaletteForegroundColor(
        QColor(m->m_fore_color[0]*255, m->m_fore_color[1]*255, m->m_fore_color[2]*255));

    renderingToolbar->repaint();
  }
}

void ToolApp::slotHelpAbout()
{
}

void ToolApp::slotStatusHelpMsg(const QString &text)
{
  ///////////////////////////////////////////////////////////////////
  // change status message of whole statusbar temporary (text, msec)
  statusBar()->message(text, 2000);
}

void ToolApp::windowMenuAboutToShow()
{
  pWindowMenu->clear();
  windowNewWindow->addTo(pWindowMenu);
  windowCascade->addTo(pWindowMenu);
  windowTile->addTo(pWindowMenu);
  if ( pWorkspace->windowList().isEmpty() )
  {
    windowAction->setEnabled(false);
  }
  else
  {
    windowAction->setEnabled(true);
  }

  pWindowMenu->insertSeparator();

  QWidgetList windows = pWorkspace->windowList();
  for ( int i = 0; i < int(windows.count()); ++i )
  {
    int id = pWindowMenu->insertItem(QString("&%1").arg(i+1)+windows.at(i)->caption(),
    this, SLOT( windowMenuActivated( int ) ) );
    pWindowMenu->setItemParameter( id, i );
    pWindowMenu->setItemChecked( id, pWorkspace->activeWindow() == windows.at(i)
 );
  }
}

void ToolApp::windowMenuActivated( int id )
{
  QWidget* w = pWorkspace->windowList().at( id );
  if ( w )
    w->setFocus();
}

//////////////////////////////////////////////////////////
// CUSTOM MENUS
//////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////
// Rendering options
//////////////////////////////////////////////////////////

void ToolApp::slotSelectBgColor()
{
  ToolView* pView = (ToolView*)pWorkspace->activeWindow();
  if(pView)
  {    
    QColor bgcolor = QColorDialog::getColor();
    if(bgcolor.isValid()){
      pView->set_clear_color(float(bgcolor.red())/255, float(bgcolor.green())/255, float(bgcolor.blue())/255);
      bgcolor_button->setPaletteBackgroundColor(bgcolor);
      bgcolor_button->setPaletteForegroundColor(bgcolor);
    }
    ToolDoc* doc = pView->getDocument();
    doc->updateAllViews(NULL);
  }
}

void ToolApp::slotSelectFgColor()
{
  ToolView* pView = (ToolView*)pWorkspace->activeWindow();
  if(pView)
  {    
    QColor fgcolor = QColorDialog::getColor();
    if(fgcolor.isValid()){
      pView->set_fore_color(float(fgcolor.red())/255, float(fgcolor.green())/255, float(fgcolor.blue())/255);
      fgcolor_button->setPaletteBackgroundColor(fgcolor);
      fgcolor_button->setPaletteForegroundColor(fgcolor);
    }
    ToolDoc* doc = pView->getDocument();
    doc->updateAllViews(NULL);
  }
}


void ToolApp::slotRenderWireframe()
{
  ToolView* pView = (ToolView*)pWorkspace->activeWindow();
  if(pView)
  {
    pView->toggleMode(ToolView::LINE);
    ToolDoc* doc = pView->getDocument();
    doc->updateAllViews(NULL);
  }
}

void ToolApp::slotRenderFill()
{
  ToolView* pView = (ToolView*)pWorkspace->activeWindow();
  if(pView)
  {    
    pView->toggleMode(ToolView::FILL);
    ToolDoc* doc = pView->getDocument();
    doc->updateAllViews(NULL);
  }
}

void ToolApp::slotRenderSuperimposing()
{
  ToolView* pView = (ToolView*)pWorkspace->activeWindow();
  if(pView)
  {
    pView->toggleSuperimposing();
    ToolDoc* doc = pView->getDocument();
    doc->updateAllViews(NULL);
  }
}

void ToolApp::slotRenderSuperimposedV()
{
  ToolView* pView = (ToolView*)pWorkspace->activeWindow();
  if(pView)
  {
    pView->toggleSuperimposeV();
    ToolDoc* doc = pView->getDocument();
    doc->updateAllViews(NULL);
  }
}


void ToolApp::slotRenderVertex()
{
  ToolView* pView = (ToolView*)pWorkspace->activeWindow();
  if(pView)
  {
    pView->toggleMode(ToolView::VERTEX);
    ToolDoc* doc = pView->getDocument();
    doc->updateAllViews(NULL);    
  }
}



void ToolApp::slotRenderEdges()
{
}

void ToolApp::slotRenderCulling()
{
  ToolView* pView = (ToolView*)pWorkspace->activeWindow();
  if(pView)
  {
    pView->toggleCulling();
    ToolDoc* doc = pView->getDocument();
    doc->updateAllViews(NULL);
  }
}

void ToolApp::slotRenderSmooth()
{
  ToolView* pView = (ToolView*)pWorkspace->activeWindow();
  if(pView)
  {
    pView->toggleSmooth();
    ToolDoc* pDoc = pView->getDocument();
    pDoc->updateAllViews(NULL);
  }
}

void ToolApp::slotRenderAntialiasing()
{
  ToolView* pView = (ToolView*)pWorkspace->activeWindow();
  if(pView)
  {
    pView->toggleAntialiasing();
    ToolDoc* doc = pView->getDocument();
    doc->updateAllViews(NULL);
  }
}


void ToolApp::slotRenderLighting()
{
  ToolView* pView = (ToolView*)pWorkspace->activeWindow();
  if(pView)
  {
    pView->toggleLighting();
    ToolDoc* doc = pView->getDocument();
    doc->updateAllViews(NULL);
  }
}

// Stam-Loop subdivision
void ToolApp::slotSubdivisionStamLoop()
{
  ToolView* pView = (ToolView*)pWorkspace->activeWindow();
  if(pView)
  {
    ToolDoc* doc = pView->getDocument();
    doc->subdivisionStamLoop();
    doc->updateAllViews(NULL);
  }
}

// sqrt3 subdivision
void ToolApp::slotSubdivisionSqrt3()
{
  ToolView* pView = (ToolView*)pWorkspace->activeWindow();
  if(pView)
  {
    ToolDoc* doc = pView->getDocument();
    doc->subdivisionSqrt3();
    doc->updateAllViews(NULL);
  }
}

// sqrt3 subdivision
void ToolApp::slotSubdivisionSqrt3_twice()
{
  ToolView* pView = (ToolView*)pWorkspace->activeWindow();
  if(pView)
  {
    ToolDoc* doc = pView->getDocument();
    doc->subdivisionSqrt3Twice();
    doc->updateAllViews(NULL);
  }
}

// Doo-Sabin subdivision
void ToolApp::slotSubdivisionDooSabin()
{
  ToolView* pView = (ToolView*)pWorkspace->activeWindow();
  if(pView)
  {
    ToolDoc* doc = pView->getDocument();
    doc->subdivisionDooSabin();
    doc->updateAllViews(NULL);
  }
}

// Catmull-Clark subdivision
void ToolApp::slotSubdivisionCatmullClark()
{
  ToolView* pView = (ToolView*)pWorkspace->activeWindow();
  if(pView)
  {
    ToolDoc* doc = pView->getDocument();
    doc->subdivisionCatmullClark();
    doc->updateAllViews(NULL);
  }
}

// Loop subdivision
void ToolApp::slotSubdivisionLoop()
{
  ToolView* pView = (ToolView*)pWorkspace->activeWindow();
  if(pView)
  {
    ToolDoc* doc = pView->getDocument();
    doc->subdivisionLoop();
    doc->updateAllViews(NULL);
  }
}




#include "tool.moc"

