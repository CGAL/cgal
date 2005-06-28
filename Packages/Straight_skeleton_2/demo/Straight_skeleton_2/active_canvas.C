// Copyright (c) 2002  Max Planck Institut fuer Informatik (Germany).
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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Radu Ursu

// if QT is not installed, a message will be issued in runtime.
#ifndef CGAL_USE_QT
#include <iostream>

int main(int, char*)
{

  std::cout << "Sorry, this demo needs QT...";
  std::cout << std::endl;

  return 0;
}

#else

#include <qplatinumstyle.h>
#include <qapplication.h>
#include <qmainwindow.h>
#include <qstatusbar.h>
#include <qfiledialog.h>
#include <qmessagebox.h>
#include <qpopupmenu.h>
#include <qmenubar.h>
#include <qtoolbutton.h>
#include <qtoolbar.h>
#include <qfiledialog.h>
#include <qtimer.h>
#include <qthread.h>
#include <qsocket.h>
#include <qserversocket.h>
#include <qstring.h>
#include <qtextstream.h>
#include <qstringlist.h>

#include "cgal_types.h"
#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_standard_toolbar.h>
#include <CGAL/IO/Qt_help_window.h>
#include <CGAL/IO/pixmaps/demoicon.xpm>

#include <fstream>
#include <iomanip>

std::ofstream sLog("/home/CGAL/demo/Active_canvas/log.txt");

class ActiveCanvasClientSocket : public QSocket
{
    Q_OBJECT
public:

    ActiveCanvasClientSocket( int sock, QObject *parent=0, const char *name=0 ) :
        QSocket( parent, name )
    {
        connect( this, SIGNAL(readyRead()), SLOT(readClient()) );
        connect( this, SIGNAL(connectionClosed()), SLOT(deleteLater()) );
        setSocket( sock );
    }

    ~ActiveCanvasClientSocket()
    {
    }

signals:
    void Command( const QString& );

private slots:
    void readClient()
    {
        QTextStream ts( this );
        while ( canReadLine() ) 
          emit Command(ts.readLine());
    }
    void deleteLater()
    {
sLog << "Connection closed" << std::endl << std::flush ;    
    }
};

class ActiveCanvasServer : public QServerSocket
{
    Q_OBJECT
public:
    ActiveCanvasServer( QObject* parent=0 ) :
        QServerSocket( 4242, 1, parent )
    {
        if ( !ok() ) {
            qWarning("Failed to bind to port 4242");
            exit(1);
        }
    }

    ~ActiveCanvasServer()
    {
    }

    void newConnection( int socket )
    {
        ActiveCanvasClientSocket *s = new ActiveCanvasClientSocket( socket, this );
sLog << "new connection" << std::endl << std::flush ;
        emit newConnect( s );
    }

signals:
    void newConnect( ActiveCanvasClientSocket* );
};



const QString my_title_string("Active Canvas");

int current_state;

#include <CGAL/IO/Qt_widget_layer.h>
#include "types.h"

class Layer : public CGAL::Qt_widget_layer
{
public:

  Layer(FigureList const& aFigures) : mFigures(aFigures) {};
  
  void draw()
  {
    widget->lock();
    
sLog << "drawing " << mFigures.size() << " figures" << std::endl ;
    
    for ( FigureList::const_iterator it = mFigures.begin(), eit = mFigures.end(); it != eit ; ++ it )
      (*it)->render(widget);
      
    widget->unlock();
  }
private:

  FigureList const& mFigures ;
};//end class

class MyWindow : public QMainWindow
{
  Q_OBJECT
public:
  MyWindow(int w, int h)
  {
    widget = new CGAL::Qt_widget(this);
    setCentralWidget(widget);

    //create a timer for checking if somthing changed
    QTimer *timer = new QTimer( this );
    connect( timer, SIGNAL(timeout()),
             this, SLOT(timerDone()) );
    timer->start( 200, FALSE );

    // file menu
    QPopupMenu * file = new QPopupMenu( this );
    menuBar()->insertItem( "&File", file );
    file->insertItem("New &Window", this, SLOT(new_window()), CTRL+Key_W);
    file->insertSeparator();
    file->insertItem( "&Close", this, SLOT(close()), CTRL+Key_X );
    file->insertItem( "&Quit", qApp, SLOT( closeAllWindows() ), CTRL+Key_Q );

    // help menu
    QPopupMenu * help = new QPopupMenu( this );
    menuBar()->insertItem( "&Help", help );
    help->insertItem("How To", this, SLOT(howto()), Key_F1);
    help->insertSeparator();
    help->insertItem("&About", this, SLOT(about()), CTRL+Key_A );
    help->insertItem("About &Qt", this, SLOT(aboutQt()) );

    //QPopupMenu * help = new QPopupMenu( this );
    menuBar()->insertItem( "&Clear", this, SLOT(clear_figures()), CTRL+Key_C );
    
    //the standard toolbar
    stoolbar = new CGAL::Qt_widget_standard_toolbar (widget, this, "ST");

    resize(w,h);
    widget->set_window(-1, 1, -1, 1);
    widget->setMouseTracking(true);

    //application flag stuff
    old_state = 0;

    mLayer.reset ( new Layer(mFigures) );
    
    widget->attach(boost::get_pointer(mLayer));
    
    ActiveCanvasServer *server = new ActiveCanvasServer( this );

    connect( server, SIGNAL(newConnect(ActiveCanvasClientSocket*)),
                SLOT(newConnect(ActiveCanvasClientSocket*)) );
 };

private:
  void something_changed(){current_state++;};

public slots:


private slots:
  
    
  void newConnect( ActiveCanvasClientSocket *s )
  {
sLog << "handling new connection" << std::endl << std::flush ;
    connect( s, SIGNAL(Command(const QString&)),
             this, SLOT(Command(const QString&)) );
    connect( s, SIGNAL(connectionClosed()),
            SLOT(connectionClosed()) );
    mClient = s ;
  }

  void connectionClosed()
  {
sLog << "connection closed." << std::endl << std::flush ;
    mClient = 0 ;
  }

  void Command(const QString& aCmd )
  {
sLog << "recieved:" << aCmd << std::endl << std::flush ;
    if ( aCmd.length() > 1 )
    {
      QString     lRawArgs = aCmd.right(aCmd.length()-1);
      QStringList lArgs    = QStringList::split(' ',lRawArgs);
      
      bool lOK = false  ;
      
      switch(aCmd[0])
      {
        case 'P' : lOK = AddPoint    (lArgs); break ;
        case 'S' : lOK = AddSegment  (lArgs); break ;
        case '~' : lOK = RemoveFigure(lArgs); break ;
      }
      
      if ( !lOK )
      {
        QString lMsg ;
        QTextOStream(&lMsg) << "Command syntax error: " << aCmd ;
        Error(lMsg);
      }
        
      widget->redraw();
    }
  }
   
  static CGAL::Color toColor ( int aC )
  {
    int r = (aC & 0xFF0000) >> 16 ;
    int g = (aC & 0x00FF00) >> 8 ;
    int b = (aC & 0x0000FF) ;
sLog << "color=" << std::hex << aC << std::dec << " r:" << r << " g:" << g << " b:" << b << std::endl ;    
    return CGAL::Color(r,g,b);
  }
  
  bool AddPoint( const QStringList& aDesc )
  {
    bool rOK = false ;
    if ( aDesc.size() == 4 )
    {
      std::size_t id = aDesc[0].toUInt  (&rOK);
      int         c  = aDesc[1].toInt   (&rOK);
      double      x  = aDesc[2].toDouble(&rOK);
      double      y  = aDesc[3].toDouble(&rOK);
      if ( rOK )
      {
        Point_2 lP(x,y);
        AddFigure(id, FigureBasePtr(new Figure<Point_2>(toColor(c),lP)) );
      }
    }
    return rOK ;
  }
  
  bool AddSegment( const QStringList& aDesc )
  {
    bool rOK = false ;
    if ( aDesc.size() == 6 )
    { 
      std::size_t id = aDesc[0].toUInt  (&rOK);
      int         c  = aDesc[1].toInt   (&rOK);
      double      sx = aDesc[2].toDouble(&rOK);
      double      sy = aDesc[3].toDouble(&rOK);
      double      tx = aDesc[4].toDouble(&rOK);
      double      ty = aDesc[5].toDouble(&rOK);
      if ( rOK )
      {
        Point_2   lSrc(sx,sy);
        Point_2   lTgt(tx,ty);
        Segment_2 lSeg(lSrc,lTgt);
        AddFigure(id, FigureBasePtr(new Figure<Segment_2>(toColor(c),lSeg)) );
      }
    }
    return rOK ;
  }
  
  bool RemoveFigure( const QStringList& aDesc )
  {
    bool rOK = false ;
    if ( aDesc.size() == 1 )
    {
      std::size_t lID = aDesc[0].toUInt(&rOK);
      if ( rOK )
      {
        if ( lID < mFigures.size() )
        {
          mFigures[lID]->set_is_visible(false);
        }
        else
        {
          QString lMsg;
          QTextOStream(&lMsg) << "There is no figure with ID " << lID << ". Cannot remove it." ;   
          Error(lMsg);
        }  
      }
    }
    return rOK ;
  }

  void AddFigure( std::size_t aID, FigureBasePtr aFigure )
  {
    if ( aID == mFigures.size() )
    {
      mFigures.push_back(aFigure);
    }
    else if ( aID < mFigures.size() )
    {
      QString lMsg;
      QTextOStream(&lMsg) << "Figure IDs must be unique, and " << aID << " was used already.";   
      Error(lMsg);
    }
    else
    {
      QString lMsg;
      QTextOStream(&lMsg) << "Sorry, in this implementation Figure IDs must be consecutive, and "
                         << aID << " is not since the last ID recieved was " << mFigures.size() << ".";   
      Error(lMsg);
    }
  }  
  void Error( const QString& aA0 )
  {
    sLog << "ERROR:" << aA0 << std::endl << std::flush ;
    QString lMsg ;
    QTextOStream(&lMsg) << "E" << aA0 ;
    SendClient(lMsg);
  }
  
  void SendClient ( const QString& aMessage )
  {
    if ( mClient != 0 )
      QTextStream(mClient) << aMessage << "\n";
  }   
  void about()
  {
    QMessageBox::about( this, my_title_string,
                        "Active Canvas\n"
                        "Copyright CGAL @2005");
  };

  void aboutQt()
  {
    QMessageBox::aboutQt( this, my_title_string );
  }

  void howto()
  {
    QString home;
    home = "help/index.html";
    CGAL::Qt_help_window *help = new
                                 CGAL::Qt_help_window(home, ".", 0, "help viewer");
    help->resize(400, 400);
    help->setCaption("Demo HowTo");
    help->show();
  }

  void clear_figures()
  {
    mFigures.clear();
    widget->redraw();
  }
  
  void new_window()
  {
    MyWindow *ed = new MyWindow(500, 500);
    ed->setCaption("View");
    ed->stoolbar->clear_history();
    ed->widget->set_window(-1.1, 1.1, -1.1, 1.1);
    ed->show();
    something_changed();
  }

  void timerDone()
  {
    if(old_state!=current_state)
    {
      widget->redraw();
      old_state = current_state;
    }
  }


private:
  CGAL::Qt_widget        *widget;
  CGAL::Qt_widget_standard_toolbar *stoolbar;
  int                    old_state;
  boost::shared_ptr<Layer>  mLayer ;
  FigureList                mFigures ;
  ActiveCanvasClientSocket* mClient ;
};


#include "active_canvas.moc"


int
main(int argc, char **argv)
{
  QApplication app( argc, argv );
  current_state = -1;

  MyWindow widget(500,500); // physical window size
  app.setMainWidget(&widget);
  widget.setCaption(my_title_string);
  widget.setMouseTracking(TRUE);
  QPixmap cgal_icon = QPixmap((const char**)demoicon_xpm);
  widget.setIcon(cgal_icon);
  widget.show();
  return app.exec();
  return 1;
}

#endif // CGAL_USE_QT
