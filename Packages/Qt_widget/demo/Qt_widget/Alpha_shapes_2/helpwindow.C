#include "helpwindow.h"


HelpWindow::HelpWindow( const QString& home_, const QString& _path,
			QWidget* parent, const char *name )
  : QMainWindow( parent, name, WDestructiveClose ),
            pathCombo( 0 )
{
  readHistory();
  readBookmarks();
  browser = new QTextBrowser( this );
  browser->mimeSourceFactory()->setFilePath( _path );
  //  browser->setFrameStyle( QFrame::Panel | QFrame::Sunken );
  setCentralWidget( browser );
  if ( !home_.isEmpty() )
       browser->setSource( home_ );

  connect( browser, SIGNAL( highlighted( const QString&) ),
           statusBar(), SLOT( message( const QString&)) );

  // The same three icons are used twice each.
  QIconSet icon_back( QPixmap("data/back.xpm") );
  QIconSet icon_forward( QPixmap("data/forward.xpm") );
  QIconSet icon_home( QPixmap("data/home.xpm") );


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



  QToolBar* toolbar = new QToolBar( this );
  addToolBar( toolbar, "Toolbar");
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

}

HelpWindow::~HelpWindow(){}

void HelpWindow::setBackwardAvailable( bool b)
{
    menuBar()->setItemEnabled( backwardId, b);
}

void HelpWindow::setForwardAvailable( bool b)
{
    menuBar()->setItemEnabled( forwardId, b);
}



void HelpWindow::print()
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

void HelpWindow::pathSelected( const QString &_path )
{
    browser->setSource( _path );
    //    if ( mHistory.values().contains(_path) )
    //	mHistory[ hist->insertItem( _path ) ] = _path;
}
void HelpWindow::histChosen( int i )
{
  if ( mHistory.contains( i ) )
    browser->setSource( mHistory[ i ] );
}

void HelpWindow::bookmChosen( int i )
{
  if ( mBookmarks.contains( i ) )
    browser->setSource( mBookmarks[ i ] );
}

void HelpWindow::addBookmark()
{
  mBookmarks[ bookm->insertItem( caption() ) ] = browser->context();
}

void HelpWindow::readHistory()
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

void HelpWindow::readBookmarks()
{
    if ( QFile::exists( QDir::currentDirPath() + "/.bookmarks" ) ) {
	QFile f( QDir::currentDirPath() + "/.bookmarks" );
	f.open( IO_ReadOnly );
	QDataStream s( &f );
	s >> bookmarks;
	f.close();
    }
}



#include "helpwindow.moc"
