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

#include<string>
#include<iostream>
#include<sstream>
#include<fstream>
#include<iomanip>
#include<list>
#include<map>

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
#include <qtextstream.h>
#include <qprogressbar.h>

#include "cgal_types.h"
#include <CGAL/Unique_hash_map.h>
#include <CGAL/IO/Color.h>
#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_standard_toolbar.h>
#include <CGAL/IO/Qt_help_window.h>
#include <CGAL/IO/pixmaps/demoicon.xpm>


class ActiveCanvasClient : QObject
{
    Q_OBJECT
public:
    ActiveCanvasClient() : mID(0)
    {
        // create the socket and connect various of its signals
        socket = new QSocket( this );
        connect( socket, SIGNAL(connected()), SLOT(socketConnected()) );
        connect( socket, SIGNAL(readyRead()), SLOT(socketReadyRead()) );
        connect( socket, SIGNAL(connectionClosed()), SLOT(socketConnectionClosed()) );
        connect( socket, SIGNAL(error(int)), SLOT(socketError(int)) );
                
        Connect();        
    }

    ~ActiveCanvasClient()
    {
    }

    void Connect()
    {
      socket->connectToHost( "localhost", 4242 );
    } 
    
    void undraw_object ( int n )
    {
      if ( !is_connected() )
        Connect();
        
      if ( is_connected() )
      {
        QString lCmd;
        QTextOStream(&lCmd) << '~' << n << '\n' ;
        sendToServer(lCmd);  
      }
    }

    int toInt( CGAL::Color color )
    {
      return ( color.red() << 16 ) + (color.green() << 8 ) + color.blue() ;
    }
    
    int draw_point ( double x, double y, CGAL::Color color, char const* layer )
    {
      int rID = -1 ;
      
      if ( !is_connected() )
        Connect();
        
      if ( is_connected() )
      {
        QString lCmd;
        QTextOStream(&lCmd) << 'P' << mID << ' ' << toInt(color) << ' ' << x << ' ' << y << '\n'  ;
        sendToServer(lCmd);  
        rID = mID++;
      }
      
      return rID ;
    }

    int draw_segment ( double sx, double sy, double tx, double ty, CGAL::Color color, char const* layer )
    {
      int rID = -1 ;
      
      if ( !is_connected() )
        Connect();
        
      if ( is_connected() )
      {
        QString lCmd;
        QTextOStream(&lCmd) << 'S' << mID << ' ' << toInt(color) << ' ' << sx << ' ' << sy << ' ' << tx << ' ' << ty << '\n' ;
        sendToServer(lCmd);  
        rID = mID++;
      }
      
      return rID ;
    }
    
private slots:
    void closeConnection()
    {
        socket->close();
        if ( socket->state() == QSocket::Closing ) 
        {
            // We have a delayed close.
            connect( socket, SIGNAL(delayedCloseFinished()),
                    SLOT(socketClosed()) );
        }
        else 
        {
            // The socket is closed.
            socketClosed();
        }
    }

    void sendToServer( const QString& aMessage )
    {
      socket->writeBlock(aMessage,aMessage.length());
      socket->flush();
    }
        
    void socketConnected()
    {
    }

    void socketConnectionClosed()
    {
    }

    void socketClosed()
    {
    }

    void socketError( int e )
    {
      QTextStream ts(socket);
      while ( socket->canReadLine() ) 
        std::cerr << "Active Canvas Server socket error: " << e << std::endl ;
    }

    void socketReadyRead()
    {
      QTextStream ts(socket);
      while ( socket->canReadLine() ) 
        std::cerr << "Active Canvas Server Response: " << ((char const*)ts.readLine()) << std::endl ;
    }
    
    bool is_connected() { return socket->state() == QSocket::Connected ; } 
    
private:

    QSocket *socket;
    int mID ;
    
};

ActiveCanvasClient sAC_Client ;

//#define CGAL_STRAIGHT_SKELETON_ENABLE_TRACE
//#define CGAL_STRAIGHT_SKELETON_ENABLE_SHOW
//#define CGAL_STRAIGHT_SKELETON_ENABLE_SHOW_AUX

#ifdef CGAL_STRAIGHT_SKELETON_ENABLE_TRACE
void Straight_skeleton_external_trace ( std::string s )
{
  static std::ofstream lout("ss_builder_log");
  lout << s << std::flush << std::endl ;
  std::printf("%s\n",s.c_str());
}
#endif

#ifdef CGAL_STRAIGHT_SKELETON_ENABLE_SHOW
void Straight_skeleton_external_undraw_object ( int n )
{
  sAC_Client.undraw_object(n);
}

int Straight_skeleton_external_draw_point ( double x, double y, CGAL::Color color, char const* layer )
{
  return sAC_Client.draw_point(x,y,color,layer);
}

int Straight_skeleton_external_draw_segment ( double sx
                                            , double sy
      				            , double tx
				            , double ty
					    , CGAL::Color color
					    , char const* layer
					    )
{
  return sAC_Client.draw_segment(sx,sy,tx,ty,color,layer);
}
#endif

#include "ss_types.h"
#include "straight_skeleton_2_toolbar.h"
#include "straight_skeleton_2_toolbar_layers.h"

const QString my_title_string("Straight_skeleton_2 Demo");

Ssds ssds;
PolygonalRegion region ;
int current_state;

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
    file->insertItem("&New", this, SLOT(new_instance()), CTRL+Key_N);
    file->insertItem("New &Window", this, SLOT(new_window()), CTRL+Key_W);
    file->insertSeparator();
    file->insertItem("&Load Polygon", this, SLOT(load_polygon()), CTRL+Key_L);
    file->insertItem("&Save Polygon", this, SLOT(save_polygon()), CTRL+Key_S);
    file->insertItem("&Save Edges", this, SLOT(save_edges()), CTRL+Key_S);
    file->insertSeparator();
    file->insertItem("Print", widget, SLOT(print_to_ps()), CTRL+Key_P);
    file->insertSeparator();
    file->insertItem( "&Close", this, SLOT(close()), CTRL+Key_X );
    file->insertItem( "&Quit", qApp, SLOT( closeAllWindows() ), CTRL+Key_Q );

    // drawing menu
    QPopupMenu * draw = new QPopupMenu( this );
    menuBar()->insertItem( "&Draw", draw );
    draw->insertItem("Generate Skeleton", this, SLOT(create_ss()), CTRL+Key_G );
    
    // help menu
    QPopupMenu * help = new QPopupMenu( this );
    menuBar()->insertItem( "&Help", help );
    help->insertItem("How To", this, SLOT(howto()), Key_F1);
    help->insertSeparator();
    help->insertItem("&About", this, SLOT(about()), CTRL+Key_A );
    help->insertItem("About &Qt", this, SLOT(aboutQt()) );

    //the standard toolbar
    stoolbar = new CGAL::Qt_widget_standard_toolbar (widget, this, "ST");
    //the new tools toolbar
    newtoolbar = new Tools_toolbar(widget, this);

    //the new scenes toolbar
    vtoolbar = new Layers_toolbar(widget, this, region, ssds);

    resize(w,h);
    widget->set_window(-1, 1, -1, 1);
    widget->setMouseTracking(true);

    //connect the widget to the main function that receives the objects
    connect(widget, SIGNAL(new_cgal_object(CGAL::Object)),
            this, SLOT(get_new_object(CGAL::Object)));

    //application flag stuff
    old_state = 0;

 };

private:
  void something_changed(){current_state++;};

public slots:
  void new_instance()
  {
    widget->lock();
    widget->clear();
    ssds.clear();
    region.clear();
    // set the Visible Area to the Interval
    widget->set_window(-1.1, 1.1, -1.1, 1.1);
    widget->unlock();
  }


private slots:

  void get_new_object(CGAL::Object obj)
  {
    PolygonPtr poly(new Polygon());
    if (CGAL::assign(*poly, obj))
    {
      CGAL::Orientation expected = ( region.size() == 0 ? CGAL::COUNTERCLOCKWISE : CGAL::CLOCKWISE ) ;
      if ( poly->orientation() != expected )
        poly->reverse_orientation();
      region.push_back(poly);  
    }
    widget->redraw();
  };


  void create_ss()
  {
    Builder builder ;
    
    for( PolygonalRegion::const_iterator bit = region.begin(), ebit = region.end() ; bit != ebit ; ++ bit )
    {
      builder.insert_CCB((*bit)->vertices_begin(),(*bit)->vertices_end());
    }  
    std::cout << "Proceesing..." << std::endl ;
    ssds = builder.proceed() ;
    std::cout << "Done." << std::endl ;
    widget->redraw();
    something_changed();
  }
    
  void about()
  {
    QMessageBox::about( this, my_title_string,
                        "Polygon partition demo\n"
                        "Copyright CGAL @2003");
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


  void save_polygon()
  {
    QString fileName = QFileDialog::getSaveFileName(
                         "sample.poly", "Polygonal PolygonalRegion files (*.poly)", this );
    if ( !fileName.isNull() && region.size() > 0 )
    {
      std::ofstream out(fileName);
      
      CGAL::set_ascii_mode(out);
     
      out << region.size() << std::endl ;
      
      for ( PolygonalRegion::const_iterator bit = region.begin(), ebit = region.end() ; bit != ebit ; ++ bit )
        out << **bit ;
        
    }
  }

  void save_edges()
  {
    QString fileName = QFileDialog::getSaveFileName(
                         "sample.edg", "CDT edges file (*.edg)", this );
    if ( !fileName.isNull() && region.size() > 0 )
    {
      std::ofstream out(fileName);
      
      CGAL::set_ascii_mode(out);
     
      std::vector<Segment> lEdges ;
      
      for ( PolygonalRegion::const_iterator bit = region.begin(), ebit = region.end() ; bit != ebit ; ++ bit )
      {
        Polygon::const_iterator first = (*bit)->vertices_begin();
        Polygon::const_iterator end   = (*bit)->vertices_end  ();
        Polygon::const_iterator last  = end - 1 ;
        for ( Polygon::const_iterator it = first ; it != end ; ++ it )
        {
          Polygon::const_iterator nx = ( it != last ? it + 1 : first ) ;
          lEdges.push_back( Segment(*it,*nx) ) ;  
        }  
      }
      
      out << lEdges.size() << '\n' ;
      for ( std::vector<Segment>::const_iterator sit = lEdges.begin(), esit = lEdges.end() ; sit != esit ; ++ sit )
        out << sit->source() << ' ' << sit->target() << '\n' ;
    }
  }
  
  void load_polygon()
  {
    QString s( QFileDialog::getOpenFileName(
                 QString::null, "Polygonal PolygonalRegion Files (*.poly)", this ) );
    if ( s.isEmpty() )
      return;
    std::ifstream in(s);
    CGAL::set_ascii_mode(in);
    region.clear();
    int ccb_count ;
    in >> ccb_count ;
    
    for ( int i = 0 ; i < ccb_count ; ++ i )
    {
      PolygonPtr poly( new Polygon() );
      in >> *poly;
      if ( i == 0 )
      {
        CGAL::Bbox_2 lBbox = poly->bbox();
        widget->set_window(lBbox.xmin(), lBbox.xmax(), lBbox.ymin(), lBbox.ymax());
      }
      CGAL::Orientation expected = ( region.size() == 0 ? CGAL::COUNTERCLOCKWISE : CGAL::CLOCKWISE ) ;
      if ( poly->orientation() != expected )
        poly->reverse_orientation();
      region.push_back(poly);
    } 
         
    ssds.clear();
    widget->redraw();
    something_changed();
  }

private:
  CGAL::Qt_widget        *widget;
  CGAL::Qt_widget_standard_toolbar *stoolbar;
  Tools_toolbar          *newtoolbar;
  Layers_toolbar         *vtoolbar;
  int                    old_state;
};

#include "straight_skeleton_2.moc"

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
