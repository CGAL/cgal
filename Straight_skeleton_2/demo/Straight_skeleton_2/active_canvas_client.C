// Copyright (c) 2005, 2006 Fernando Luis Cacciola Carballal. All rights reserved.
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
// $URL: svn+ssh://fcacciola@scm.gforge.inria.fr/svn/cgal/trunk/Straight_skeleton_2/demo/Straight_skeleton_2/straight_skeleton_2.C $
// $Id: straight_skeleton_2.C 29640 2006-03-20 21:25:42Z fcacciola $
//
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>
//

#include <qsocket.h>

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
