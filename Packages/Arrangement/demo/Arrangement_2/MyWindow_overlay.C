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
     
	 std::vector<QColor> ubf_colors(20) ; // vector of colors of the unbounded faces od the planar maps //..//
	 Qt_widget_demo_tab<Segment_tab_traits> *w_demo_p;
     
     std::list<Pm_seg_2> seg_list;
     std::list<Pm_seg_2> antenna_list;
     Pm_seg_const_iter itp;
     int current, flag;
     bool all_pm_are_empty = true; //..//
     while (! indexes.empty())
     {
       current = indexes.front();
       indexes.pop_front();
	     flag = paint_flags.front();
       paint_flags.pop_front();
       w_demo_p = 
         static_cast<Qt_widget_demo_tab<Segment_tab_traits> *> 
         (myBar->page( current ));
       
	   if (! w_demo_p->empty)//..//
	   {
         Seg_arr::Edge_iterator hei;
         for (hei = w_demo_p->m_curves_arr->edges_begin(); 
              hei != w_demo_p->m_curves_arr->edges_end(); ++hei) 
         {
           Pm_xseg_2 xcurve = hei->curve();
		       Curve_data cd = xcurve.get_data();
           if ( flag == 0 )
		         cd.m_index = w_demo_p_new->index;
		       xcurve.set_data( cd );
           Pm_seg_2 curve(xcurve, cd);
           if ( w_demo_p->antenna(*hei))
             antenna_list.push_back(curve);
           else
             seg_list.push_back(curve);
         }// for
		     all_pm_are_empty = false;
	   }// if //..//
     w_demo_p_new->bbox = w_demo_p_new->bbox + w_demo_p->bbox;
       
	   // update the vector of colors of unbounded faces
	   ubf_colors[current] = w_demo_p->unbounded_face_color();//..//
	 }// while

	 if (! all_pm_are_empty)//..//
	 {
     if ( ! seg_list.empty())
       w_demo_p_new->m_curves_arr->insert(seg_list.begin(),seg_list.end());
     w_demo_p_new->update_colors(ubf_colors); 
     Seg_notification  seg_notif;
     std::list<Pm_seg_2>::iterator it;
     for (it = antenna_list.begin(); it != antenna_list.end(); it++)
       w_demo_p_new->m_curves_arr->insert( *it, &seg_notif );

     w_demo_p_new->empty = false;

     if (!colors_flag)
     {
       // update new planner map indexes
       Seg_arr::Halfedge_iterator hei;
       for (hei = w_demo_p_new->m_curves_arr->halfedges_begin(); 
           hei != w_demo_p_new->m_curves_arr->halfedges_end(); ++hei) 
       {
         Curve_data cd = hei->curve().get_data();
         cd.m_type = Curve_data::INTERNAL;
         cd.m_index = w_demo_p_new->index;
         hei->curve().set_data( cd );
       }
	   }

	   
	 w_demo_p_new->set_window(w_demo_p_new->bbox.xmin() , 
                 w_demo_p_new->bbox.xmax() ,w_demo_p_new->bbox.ymin() , 
	                                         w_demo_p_new->bbox.ymax());
	 
	 }
	 else//..//
	 {
		
    w_demo_p_new->update_colors(ubf_colors);   
		//w_demo_p_new->set_window(-10, 10, -10, 10);
	 }

	 // update the colors of the faces of the new PM //..//
	
	 
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

  	 std::vector<QColor> ubf_colors(20) ; // vector of colors of the unbounded faces od the planar maps //..//

     Qt_widget_demo_tab<Polyline_tab_traits> *w_demo_p;
     
     std::list<Pm_pol_2> pol_list;
     std::list<Pm_pol_2> antenna_list;
     Pm_pol_const_iter itp;
     int current,flag;
     bool all_pm_are_empty = true; //..//
     while (! indexes.empty())
     {
       current = indexes.front();
       indexes.pop_front();
   	   flag = paint_flags.front();
       paint_flags.pop_front();

       w_demo_p = 
         static_cast<Qt_widget_demo_tab<Polyline_tab_traits> *> 
         (myBar->page( current ));
       
   	   if (! w_demo_p->empty)//..//
       {
         Pol_arr::Edge_iterator hei;
         for (hei = w_demo_p->m_curves_arr->edges_begin(); 
            hei != w_demo_p->m_curves_arr->edges_end(); ++hei) 
         {
           Pm_xpol_2 xcurve = hei->curve();
 		       Curve_pol_data cd = xcurve.get_data();
           if ( flag == 0 )
		         cd.m_index = w_demo_p_new->index;
           xcurve.set_data( cd );
           Pm_pol_2 curve(xcurve, cd);
            if ( w_demo_p->antenna(*hei))
             antenna_list.push_back(curve);
           else
             pol_list.push_back(curve);
         }// for
         all_pm_are_empty = false; //..//
       }// if
       w_demo_p_new->bbox = w_demo_p_new->bbox + w_demo_p->bbox;
   	   
       // update the vector of colors of unbounded faces
  	   ubf_colors[current] = w_demo_p->unbounded_face_color();//..//

     }// while
     
     if (! all_pm_are_empty)//..//
	   {
       if ( ! pol_list.empty())
       w_demo_p_new->m_curves_arr->insert(
                                       pol_list.begin(),pol_list.end());
       w_demo_p_new->update_colors(ubf_colors);   
       Pol_notification pol_notif;
       std::list<Pm_pol_2>::iterator it;
       for (it = antenna_list.begin(); it != antenna_list.end(); it++)
         w_demo_p_new->m_curves_arr->insert( *it , &pol_notif );

          
       if (!colors_flag)
       {
         // update new planner map indexes
         Pol_arr::Halfedge_iterator hei;
         for (hei = w_demo_p_new->m_curves_arr->halfedges_begin(); 
             hei != w_demo_p_new->m_curves_arr->halfedges_end(); ++hei) 
         {
           Curve_pol_data cd = hei->curve().get_data();
           cd.m_type = Curve_pol_data::INTERNAL;
           cd.m_index = w_demo_p_new->index;
           hei->curve().set_data( cd );
         }// for 
       }// if
	     w_demo_p_new->set_window(w_demo_p_new->bbox.xmin() , 
                 w_demo_p_new->bbox.xmax() ,w_demo_p_new->bbox.ymin() , 
	                                         w_demo_p_new->bbox.ymax());
     }
	   else//..//
	   {
		   w_demo_p_new->update_colors(ubf_colors);   
       //w_demo_p_new->set_window(-10, 10, -10, 10);
	   }

     w_demo_p_new->setCursor(old);

     break;
    }
   case CONIC_TRAITS:
    {
	 if(new_tab)
	   add_conic_tab();
	 else
	   updateTraitsType( setConicTraits );

     Qt_widget_demo_tab<Conic_tab_traits> *w_demo_p_new = 
       static_cast<Qt_widget_demo_tab<Conic_tab_traits> *> 
       (myBar->currentPage());

	 QCursor old = w_demo_p_new->cursor();
     w_demo_p_new->setCursor(Qt::WaitCursor);
     w_demo_p_new->setCursor(Qt::WaitCursor);
    
   std::vector<QColor> ubf_colors(20) ; // vector of colors of the unbounded faces od the planar maps //..//
	 Qt_widget_demo_tab<Conic_tab_traits> *w_demo_p;
     std::list<Pm_conic_2> antenna_list;
     Pm_xconic_const_iter itp;
     int current,flag;
     bool all_pm_are_empty = true; //..//
     while (! indexes.empty())
     {
       current = indexes.front();
       indexes.pop_front();
   	   flag = paint_flags.front();
       paint_flags.pop_front();

       w_demo_p = static_cast<Qt_widget_demo_tab<Conic_tab_traits> *>
         (myBar->page( current ));
     
       if (! w_demo_p->empty)//..//
  	   {

          Conic_arr::Edge_iterator hei;
          for (hei = w_demo_p->m_curves_arr->edges_begin(); 
                hei != w_demo_p->m_curves_arr->edges_end(); ++hei) 
          {
            Pm_xconic_2 xcurve = hei->curve();
  		      Curve_conic_data cd = xcurve.get_data();
            if ( flag == 0 )
		        cd.m_index = w_demo_p_new->index;
            xcurve.set_data( cd );
            Pm_conic_2 curve(xcurve, cd);
            if ( w_demo_p->antenna(*hei))
              antenna_list.push_back(curve);
            else
              w_demo_p_new->m_curves_arr->insert(curve);
          }// for
		      all_pm_are_empty = false;
       }// if
	   w_demo_p_new->bbox = w_demo_p_new->bbox + w_demo_p->bbox;

 	   // update the vector of colors of unbounded faces
	   ubf_colors[current] = w_demo_p->unbounded_face_color();//..//

     }// while

 	    if (! all_pm_are_empty)//..//
	    {
        //w_demo_p_new->m_curves_arr->insert(seg_list.begin(),seg_list.end());
        // update the colors of the faces of the new PM //..//
  	    w_demo_p_new->update_colors(ubf_colors);  
        Conic_notification conic_notif;
        std::list<Pm_conic_2>::iterator it;
        for (it = antenna_list.begin(); it != antenna_list.end(); it++)
          w_demo_p_new->m_curves_arr->insert( *it , &conic_notif );
        w_demo_p_new->empty = false;

	      if (!colors_flag)
        {
          // update new planner map indexes
          Conic_arr::Halfedge_iterator hei;
          for (hei = w_demo_p_new->m_curves_arr->halfedges_begin(); 
                hei != w_demo_p_new->m_curves_arr->halfedges_end(); ++hei) 
          {
            Curve_conic_data cd = hei->curve().get_data();
            cd.m_type = Curve_conic_data::INTERNAL;
            cd.m_index = w_demo_p_new->index;
            hei->curve().set_data( cd );
          }
        }
        w_demo_p_new->set_window(w_demo_p_new->bbox.xmin() , 
                    w_demo_p_new->bbox.xmax() ,w_demo_p_new->bbox.ymin() , 
	                                            w_demo_p_new->bbox.ymax());
      }
      else
      {
     	    w_demo_p_new->update_colors(ubf_colors);   
	  	    //w_demo_p_new->set_window(-10, 10, -10, 10);
	    }

	   	 
      w_demo_p_new->setCursor(old);
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
