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

#include <qapplication.h>
#include <qfont.h>
#include <qstring.h>
#include <qtextcodec.h>
#include <qtranslator.h>

#include "tool.h"
#include "icons/cgal.xpm"

int main(int argc, char *argv[])
{
  QApplication a(argc, argv);
  //a.setFont(QFont("helvetica", 12));
  QTranslator tor( 0 );
  // set the location where your .qm files are in load() below as the last parameter instead of "."
  // for development, use "/" to use the english original as
  // .qm files are stored in the base project directory.
  tor.load( QString("tool.") + QTextCodec::locale(), "." );
  a.installTranslator( &tor );

  ToolApp *tool=new ToolApp();
  QPixmap cgalIcon = QPixmap(cgal_xpm);
  tool->setIcon(cgalIcon);
  a.setMainWidget(tool);

  tool->show();

  if(argc>1)
    tool->openDocumentFile(argv[1]);
	//else
	  //tool->openDocumentFile();
	
  return a.exec();
}
