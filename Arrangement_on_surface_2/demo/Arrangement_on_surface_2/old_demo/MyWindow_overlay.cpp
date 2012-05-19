// Copyright (c) 2005  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
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
// $URL: svn+ssh://scm.gforge.inria.fr/svn/cgal/branches/features/gsoc2012-Arrangement_on_surface_2-demo-atsui/Arrangement_on_surface_2/demo/Arrangement_on_surface_2/MyWindow_overlay.cpp $
// $Id: MyWindow_overlay.cpp 67117 2012-01-13 18:14:48Z lrineau $
//
//
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>

#include <CGAL/basic.h>


#include "arrangement_2.h"
#include "forms.h"
#include "qt_layer.h"
#include "demo_tab.h"
#include "overlay_functor.h"

#include <CGAL/Arr_overlay_2.h>

/*! open the overlay dialog form and read its output */
void MyWindow::overlay_pm()
{
  OverlayForm *form = new OverlayForm( myBar , this ,tab_number);
  if ( form->exec() )
  {
    unsigned int i = 2;
    if (form->listBox2->count() < i)
    {
      QMessageBox::information( this, "Overlay",
          "You need more than one arrangement to perform overlay...");
      return;
    }
    i = 12;
    if (form->listBox2->count() > i)
    {
      QMessageBox::information
        (this, "Overlay",
         "Maximal number of arrangements to overlay is 12 !");
      return;
    }
    CheckForm *check_form = new CheckForm( form , this );
    if ( ! check_form->exec() )
      return;

    std::list<int>      indexes;
    std::list<int>      paint_flags;
    TraitsType          t;
    Qt_widget_base_tab *w_demo_p;
    int                 index,real_index = 0;

    for (unsigned int i = 0; i < form->listBox2->count(); i++)
    {
      form->listBox2->setCurrentItem(i);
      char s[100];
      strcpy(s, form->listBox2->currentText());
      char * pch;
      pch = strtok(s," ");
      pch = strtok(NULL, " ");
      index = atoi(pch);
      real_index = realIndex(index-1);
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
void MyWindow::make_overlay(std::list<int> indexes,
                            std::list<int> /* paint_flags */, TraitsType t,
                            bool new_tab)
{
  switch ( t )
  {
    case SEGMENT_TRAITS:
    {
       if(new_tab)
         add_segment_tab();
       else
         updateTraitsType( setSegmentTraits );
        Qt_widget_demo_tab<Segment_tab_traits> *w_demo_p_new =
         static_cast<Qt_widget_demo_tab<Segment_tab_traits> *>
       (myBar->currentPage());

   QCursor old = w_demo_p_new->cursor();
     w_demo_p_new->setCursor(Qt::WaitCursor);

   Qt_widget_demo_tab<Segment_tab_traits> *w_demo_p1, *w_demo_p2;

     int ind1, ind2;

     ind1 = indexes.front();
     ind2 = indexes.back();

     w_demo_p1 =
       static_cast<Qt_widget_demo_tab<Segment_tab_traits> *>
       (myBar->page( ind1 ));

     w_demo_p2 =
       static_cast<Qt_widget_demo_tab<Segment_tab_traits> *>
       (myBar->page( ind2 ));

     w_demo_p_new->bbox = w_demo_p1->bbox + w_demo_p2->bbox;

     // update the vector of colors of unbounded faces

     typedef  Segment_tab_traits::Arrangement_2    Arrangement_2;
     Overlay_functor<Arrangement_2> func;
     CGAL::overlay(*(w_demo_p1->m_curves_arr),
                   *(w_demo_p2->m_curves_arr),
                   *(w_demo_p_new->m_curves_arr),
                   func);

   w_demo_p_new->set_window(w_demo_p_new->bbox.xmin(),
                            w_demo_p_new->bbox.xmax(),
                            w_demo_p_new->bbox.ymin() ,
                            w_demo_p_new->bbox.ymax());

     w_demo_p_new->setCursor(old);
     break;
     }

    case POLYLINE_TRAITS:
    {
      if(new_tab)
        add_polyline_tab();
      else
        updateTraitsType( setPolylineTraits );


       Qt_widget_demo_tab<Polyline_tab_traits> *w_demo_p_new =
         static_cast<Qt_widget_demo_tab<Polyline_tab_traits> *>
       (myBar->currentPage());

   QCursor old = w_demo_p_new->cursor();
     w_demo_p_new->setCursor(Qt::WaitCursor);

   Qt_widget_demo_tab<Polyline_tab_traits> *w_demo_p1, *w_demo_p2;

     int ind1, ind2;

     ind1 = indexes.front();
     ind2 = indexes.back();

     w_demo_p1 =
       static_cast<Qt_widget_demo_tab<Polyline_tab_traits> *>
       (myBar->page( ind1 ));

     w_demo_p2 =
       static_cast<Qt_widget_demo_tab<Polyline_tab_traits> *>
       (myBar->page( ind2 ));

     w_demo_p_new->bbox = w_demo_p1->bbox + w_demo_p2->bbox;

     // update the vector of colors of unbounded faces

     typedef  Polyline_tab_traits::Arrangement_2    Arrangement_2;
     Overlay_functor<Arrangement_2> func;
     CGAL::overlay(*(w_demo_p1->m_curves_arr),
                   *(w_demo_p2->m_curves_arr),
                   *(w_demo_p_new->m_curves_arr),
                   func);

   w_demo_p_new->set_window(w_demo_p_new->bbox.xmin(),
                            w_demo_p_new->bbox.xmax(),
                            w_demo_p_new->bbox.ymin() ,
                            w_demo_p_new->bbox.ymax());

     w_demo_p_new->setCursor(old);
     break;
    }

    case CONIC_TRAITS:
    {
#ifdef CGAL_USE_CORE
      if(new_tab)
        add_conic_tab();
      else
        updateTraitsType( setConicTraits );

       Qt_widget_demo_tab<Conic_tab_traits> *w_demo_p_new =
         static_cast<Qt_widget_demo_tab<Conic_tab_traits> *>
       (myBar->currentPage());

   QCursor old = w_demo_p_new->cursor();
     w_demo_p_new->setCursor(Qt::WaitCursor);

   Qt_widget_demo_tab<Conic_tab_traits> *w_demo_p1, *w_demo_p2;

     int ind1, ind2;

     ind1 = indexes.front();
     ind2 = indexes.back();

     w_demo_p1 =
       static_cast<Qt_widget_demo_tab<Conic_tab_traits> *>
       (myBar->page( ind1 ));

     w_demo_p2 =
       static_cast<Qt_widget_demo_tab<Conic_tab_traits> *>
       (myBar->page( ind2 ));

     w_demo_p_new->bbox = w_demo_p1->bbox + w_demo_p2->bbox;

     // update the vector of colors of unbounded faces

     typedef  Conic_tab_traits::Arrangement_2    Arrangement_2;
     Overlay_functor<Arrangement_2> func;
     CGAL::overlay(*(w_demo_p1->m_curves_arr),
                   *(w_demo_p2->m_curves_arr),
                   *(w_demo_p_new->m_curves_arr),
                   func);



   w_demo_p_new->set_window(w_demo_p_new->bbox.xmin(),
                            w_demo_p_new->bbox.xmax(),
                            w_demo_p_new->bbox.ymin() ,
                            w_demo_p_new->bbox.ymax());

     w_demo_p_new->setCursor(old);
#endif

     break;

    }
  }
}



/*! real index - finds the tab bar index of a tab
 * \param index - the tab index (the same one it was
 *                initialized with)
 *\ return the tab bar index of this tab
 */
int MyWindow::realIndex(int index)
{
  Qt_widget_base_tab * w_demo_p;
  for (int i = 0; i < tab_number; i++)
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


