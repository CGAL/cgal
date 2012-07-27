// Copyright (c) 2002-2004  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Radu Ursu

#include <CGAL/basic.h>

#include <CGAL/IO/Qt_help_window.h>

/* XPM */
/* Drawn  by Mark Donohoe for the K Desktop Environment */
/* See http://www.kde.org */
static const char*backb[]={
"16 16 5 1",
"# c #000000",
"a c #ffffff",
"c c #808080",
"b c #c0c0c0",
". c None",
"................",
".......#........",
"......##........",
".....#a#........",
"....#aa########.",
"...#aabaaaaaaa#.",
"..#aabbbbbbbbb#.",
"...#abbbbbbbbb#.",
"...c#ab########.",
"....c#a#ccccccc.",
".....c##c.......",
"......c#c.......",
".......cc.......",
"........c.......",
"................",
"......................"};

/* XPM */
/* Drawn  by Mark Donohoe for the K Desktop Environment */
/* See http://www.kde.org */
static const char*forwardb[]={
"16 16 5 1",
"# c #000000",
"a c #ffffff",
"c c #808080",
"b c #c0c0c0",
". c None",
"................",
"................",
".........#......",
".........##.....",
".........#a#....",
"..########aa#...",
"..#aaaaaaabaa#..",
"..#bbbbbbbbbaa#.",
"..#bbbbbbbbba#..",
"..########ba#c..",
"..ccccccc#a#c...",
"........c##c....",
"........c#c.....",
"........cc......",
"........c.......",
"................",
"................"};

/* XPM */
/* Drawn  by Mark Donohoe for the K Desktop Environment */
/* See http://www.kde.org */
static const char*homeb[]={
"16 16 4 1",
"# c #000000",
"a c #ffffff",
"b c #c0c0c0",
". c None",
"........... ....",
"   ....##.......",
"..#...####......",
"..#..#aabb#.....",
"..#.#aaaabb#....",
"..##aaaaaabb#...",
"..#aaaaaaaabb#..",
".#aaaaaaaaabbb#.",
"###aaaaaaaabb###",
"..#aaaaaaaabb#..",
"..#aaa###aabb#..",
"..#aaa#.#aabb#..",
"..#aaa#.#aabb#..",
"..#aaa#.#aabb#..",
"..#aaa#.#aabb#..",
"..#####.######..",
"................"};

namespace CGAL{

Qt_help_window::Qt_help_window( const QString& home_, const QString& _path,
			QWidget* parent, const char *name )
  : QMainWindow( parent, name, WDestructiveClose ),
            pathCombo( 0 )
{
  readHistory();
  browser = new QTextBrowser( this );
  browser->mimeSourceFactory()->setFilePath( _path );
  //  browser->setFrameStyle( QFrame::Panel | QFrame::Sunken );
  setCentralWidget( browser );
  if ( !home_.isEmpty() )
       browser->setSource( home_ );

  connect( browser, SIGNAL( highlighted( const QString&) ),
           statusBar(), SLOT( message( const QString&)) );

  // The same three icons are used twice each.
  QIconSet icon_back( QPixmap((const char**)backb) );
  QIconSet icon_forward( QPixmap((const char**)forwardb) );
  QIconSet icon_home( QPixmap((const char**)homeb) );


  QPopupMenu* file = new QPopupMenu( this );
  file->insertItem( "&Print", this, SLOT( print() ), CTRL+Key_P );
  file->insertSeparator();
  file->insertItem( "&Close", this, SLOT( close() ), CTRL+Key_Q );

  QPopupMenu* go = new QPopupMenu( this );
  backwardId = go->insertItem( icon_back,
				 "&Backward", browser, SLOT( backward() ),
				 CTRL+Key_Left );
  forwardId = go->insertItem( icon_forward,
				"&Forward", browser, SLOT( forward() ),
				CTRL+Key_Right );
  go->insertItem( icon_home, "&Home", browser, SLOT( home() ) );

  menuBar()->insertItem("&File", file);
  menuBar()->insertItem("&Go", go); 

  QToolBar* toolbar = new QToolBar( "Toolbar", this, this);

  QToolButton* button;


  button = new QToolButton( icon_back, "Backward", "", 
                            browser, SLOT(backward()), toolbar );
  connect( browser, SIGNAL( backwardAvailable(bool) ), 
           button, SLOT( setEnabled(bool) ) );
  button->setEnabled( FALSE );
  button = new QToolButton( icon_forward, "Forward", "", 
                            browser, SLOT(forward()), toolbar );
  connect( browser, SIGNAL( forwardAvailable(bool) ), 
	   button, SLOT(setEnabled(bool) ) );
  button->setEnabled( FALSE );
  button = new QToolButton( icon_home, "Home", "", browser,
			    SLOT(home()), toolbar );
  toolbar->addSeparator();

  pathCombo = new QComboBox( TRUE, toolbar );
  connect( pathCombo, SIGNAL( activated( const QString & ) ),
	     this, SLOT( pathSelected( const QString & ) ) );
  toolbar->setStretchableWidget( pathCombo );
  setRightJustification( TRUE );
  pathCombo->insertItem( home_ );
}

Qt_help_window::~Qt_help_window(){}

void Qt_help_window::setBackwardAvailable( bool b)
{
    menuBar()->setItemEnabled( backwardId, b);
}

void Qt_help_window::setForwardAvailable( bool b)
{
    menuBar()->setItemEnabled( forwardId, b);
}



void Qt_help_window::print()
{
#ifndef QT_NO_PRINTER
    QPrinter printer;
    printer.setFullPage(TRUE);
    if ( printer.setup( this ) ) {
	QPainter p( &printer );
	QPaintDeviceMetrics metrics(p.device());
	int dpix = metrics.logicalDpiX();
	int dpiy = metrics.logicalDpiY();
	const int margin = 72; // pt
	QRect body(margin*dpix/72, margin*dpiy/72,
		   metrics.width()-margin*dpix/72*2,
		   metrics.height()-margin*dpiy/72*2 );
	QSimpleRichText richText( browser->text(), QFont(), browser->context(), browser->styleSheet(),
				  browser->mimeSourceFactory(), body.height() );
	richText.setWidth( &p, body.width() );
	QRect view( body );
	int page = 1;
	do {
	    richText.draw( &p, body.left(), body.top(), view, colorGroup() );
	    view.moveBy( 0, body.height() );
	    p.translate( 0 , -body.height() );
	    p.drawText( view.right() - p.fontMetrics().width( QString::number(page) ),
			view.bottom() + p.fontMetrics().ascent() + 5, QString::number(page) );
	    if ( view.top()  >= richText.height() )
		break;
	    printer.newPage();
	    page++;
	} while (TRUE);
    }
    #endif
}

void Qt_help_window::pathSelected( const QString &_path )
{
    browser->setSource( _path );
}
void Qt_help_window::histChosen( int i )
{
  if ( mHistory.contains( i ) )
    browser->setSource( mHistory[ i ] );
}

void Qt_help_window::readHistory()
{
    if ( QFile::exists( QDir::currentDirPath() + "/.history" ) ) {
	QFile f( QDir::currentDirPath() + "/.history" );
	f.open( IO_ReadOnly );
	QDataStream s( &f );
	s >> history;
	f.close();
	while ( history.count() > 20 )
	    history.remove( history.begin() );
    }
}
} //end CGAL namespace

#include "Qt_help_window.moc"
