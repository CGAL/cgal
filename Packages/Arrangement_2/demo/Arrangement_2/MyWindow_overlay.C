// Copyright (c) 2005  Tel-Aviv University (Israel).
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
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>

#include <CGAL/basic.h>

#ifdef CGAL_USE_QT

#include "MyWindow.h"
#include "forms.h"
#include "qt_layer.h"
#include "demo_tab.h"

/*! open the overlay dialog form and read its output */
void MyWindow::overlay_pm()
{
  OverlayForm *form = new OverlayForm( myBar , this ,tab_number - 1);
  if ( form->exec() ) 
  {    
    unsigned int i = 2;
    if (form->listBox2->count() < i)
    {
      QMessageBox::information( this, "Overlay",
          "Please!!! you need more than one planar map to make an overlay...");
      return;
    }
    i = 12;
    if (form->listBox2->count() > i)
    {
      QMessageBox::information( this, "Overlay",
                              "Max number of Planar Maps to overlay is 12!!!");
      return;
    }
	CheckForm *check_form = new CheckForm( form , this );
	if ( ! check_form->exec() )
	  return;
    
	std::list<int> indexes;
	std::list<int> paint_flags;
    TraitsType t;
    Qt_widget_base_tab *w_demo_p;
    int index,real_index;
    
    for (unsigned int i = 0; i < form->listBox2->count(); i++)
    {
      form->listBox2->setCurrentItem(i);
      char s[100];
      strcpy(s, form->listBox2->currentText()); 
      char * pch;
      pch = strtok(s," ");
      pch = strtok(NULL, " ");
      index = atoi(pch);
      real_index = realIndex(index);
      indexes.push_back(real_index);
	  QCheckBox *b = 
		  static_cast<QCheckBox *> (check_form->button_group->find(i));
	  if ( b->isChecked() )
        paint_flags.push_back(1);
	  else
	    paint_flags.push_back(0);
    }
	delete check_form;
    
    w_demo_p = static_cast<Qt_widget_base_tab *> (myBar->page( real_index ));
    t = w_demo_p->traits_type;

    FileOpenOptionsForm * form = new FileOpenOptionsForm(false);
    if ( form->exec() ) 
    {
		int id = form->buttonGroup->id(form->buttonGroup->selected());
		switch ( id ) 
		{
			case 0: // open file in a new tab
			make_overlay( indexes , paint_flags , t , true);    
			break;
			case 1: // open file in current tab (delete current Pm)
			make_overlay( indexes , paint_flags , t , false);
			break;        
		}// switch
    }// if
	
  }
  delete form;
}

/*! make the overlay
 * \param indexes - list of the planar maps indexes in the overlay
 * \param t - the traits type of the overlay
 */
void MyWindow::make_overlay( std::list<int> indexes ,
               std::list<int> paint_flags ,TraitsType t , bool new_tab)
{}

/*! real index - finds the tab bar index of a tab
 * \param index - the tab index (the same one it was 
 *                initialized with)
 *\ return the tab bar index of this tab
 */
int MyWindow::realIndex(int index)
{
  Qt_widget_base_tab * w_demo_p;
  for (int i = 0; i < tab_number-1; i++)
  {
    if ( myBar->isTabEnabled( myBar->page(i) ) )
    {
      // We peform downcasting from QWigdet* to 
      // Qt_widget_base_tab*, as we know that only
      // Qt_widget_base_tab objects are stored in the tab pages.
      w_demo_p = static_cast<Qt_widget_base_tab *> (myBar->page(i));
      if (w_demo_p->index == index)
        return i;
    }
  }
  return -1;
}


#endif // CGAL_USE_QT
