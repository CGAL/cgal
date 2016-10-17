// ============================================================================
//
// Copyright (c) 2001-2006 Max-Planck-Institut Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of EXACUS (http://www.mpi-inf.mpg.de/projects/EXACUS/).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// ----------------------------------------------------------------------------
//
// Library       : AlciX
// File          : demos/xalci/include/misc.h
// SoX_release   : $Name:  $
// Revision      : $Revision: 1.2 $
// Revision_date : $Date: 2009-06-30 13:14:59 $
//
// Author(s)     : Pavel Emeliyanenko <asm@mpi-sb.mpg.de>
//
// ============================================================================

#ifndef MISC_H
#define MISC_H

#include <qlistbox.h>
#include <qpushbutton.h>
#include <qlayout.h>
#include <qvbox.h>

#include <sstream>

#include <CGAL/IO/Qt_widget_standard_toolbar.h>
#include <CGAL/IO/Qt_help_window.h>
#include <CGAL/IO/Qt_widget_layer.h>

/*!\brief
 * declares miscellaneous routines in a separate file (to speed-up
 * compilation)
 */

class Curve_selection_dialog : public QDialog {

public:

    QPushButton* ok_btn,*cancel_btn;
    QListBox* curve_list;

    template<typename InputIterator>
    Curve_selection_dialog(QWidget* widget, const char* name,
            bool modal, InputIterator begin, InputIterator end)
    {
        QHBoxLayout *hbox = new QHBoxLayout(this,10,10);
        curve_list = new QListBox(this);
        curve_list->setSelectionMode(QListBox::Multi);
        hbox->addWidget(curve_list,8);

        int i = 0;
        for(InputIterator it = begin;it != end; it++, i++) {
            std::stringstream curr_item;
            ::CGAL::set_pretty_mode(curr_item);
            curr_item << i << "th: ";
            curr_item << (*it);
            curve_list->insertItem(curr_item.str());
        }

        //this->addWidget(curve_list);
        QVBox *vbox = new QVBox(this);
        hbox->addWidget(vbox);
        ok_btn = new QPushButton("Analyse",vbox);
        cancel_btn = new QPushButton("Cancel",vbox);

        connect(ok_btn, SIGNAL(clicked()), SLOT(accept()));
        connect(cancel_btn, SIGNAL(clicked()), SLOT(reject()));
    }

};

class Graphic_layer : public CGAL::Qt_widget_layer
{
Q_OBJECT
 
public:
    Graphic_layer(CGAL::Qt_widget *attach_to, int index_,
        int color_index_, int fst_ind = -1,
        int snd_ind = -1, QObject* parent = NULL, const char* name = NULL) :
         Qt_widget_layer(parent, name),
          erase(false), index(index_), color_index(color_index_),
          first_curve_index(fst_ind), second_curve_index(snd_ind)
    { 
        attach_to->attach(this);
    }

    void draw();

    int get_first_index() { return this->first_curve_index; }
    int get_second_index() { return this->second_curve_index; }

    bool erase;
protected:

    int index, color_index;
    int first_curve_index, second_curve_index;
        
};

typedef std::vector<Graphic_layer *> Layers;

#endif // MISC_H
