// ======================================================================
//
// Copyright (c) 1997-2002 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// file          : include/CGAL/IO/Qt_widget_helpwindow.h
// package       : Qt_widget (1.2.46)
// author(s)     : Radu Ursu
// release       : $CGAL_Revision: CGAL-2.5-I-60 $
// release_date  : $CGAL_Date: 2002/10/18 $
//
// coordinator   : Laurent Rineau <rineau@clipper.ens.fr>
//
// ======================================================================

#ifndef CGAL_QT_WIDGET_HELPWINDOW_H
#define CGAL_QT_WIDGET_HELPWINDOW_H

#include <qmainwindow.h>
#include <qtextbrowser.h>
#include <qstringlist.h>
#include <qmap.h>
#include <qdir.h>
#include <qstatusbar.h>
#include <qpixmap.h>
#include <qpopupmenu.h>
#include <qmenubar.h>
#include <qtoolbar.h>
#include <qtoolbutton.h>
#include <qiconset.h>
#include <qfile.h>
#include <qstylesheet.h>
#include <qmessagebox.h>
#include <qfiledialog.h>
#include <qapplication.h>
#include <qcombobox.h>
#include <qevent.h>
#include <qlineedit.h>
#include <qobjectlist.h>
#include <qfileinfo.h>
#include <qfile.h>
#include <qdatastream.h>
#include <qprinter.h>
#include <qsimplerichtext.h>
#include <qpainter.h>
#include <qpaintdevicemetrics.h>




class QComboBox;
class QPopupMenu;

class HelpWindow : public QMainWindow
{
    Q_OBJECT
public:
    HelpWindow( const QString& home_,  const QString& path, 
                QWidget* parent = 0, const char *name=0 );
    ~HelpWindow();
private slots:
    void setBackwardAvailable( bool );
    void setForwardAvailable( bool );

    void print();

    void pathSelected( const QString & );
    void histChosen( int );

private:
    void readHistory();

    QTextBrowser* browser;
    QComboBox *pathCombo;
    int backwardId, forwardId;
    QStringList history;
    QPopupMenu *hist;
    QMap<int, QString> mHistory;
};

static char* back[100];
static char* forward[100];
static char* home[100];

#endif
