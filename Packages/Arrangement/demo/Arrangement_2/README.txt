README - Planar Map Demo

This file contains a general description of the demo program structure and files.

structure:
the main function creates a MyWindow object which is the main window. a MyWindow
object has QTabWidget data member that holds the tab widgets.
the widget are Qt_widget_demo_tab, which inherit from CGAL's Qt_widget and gets
a template parameter of the tab's traits (segment, polyline or conic).

-----------------------------------------------------------------
|                           MyWindow                            |
|											                    |
|				  ----------------------------------------      |
| QTabWidget ->  1|	Qt_widget_demo_tab< Segment_Traits > |      |  
|				  |--------------------------------------|      |  
|                2|	Qt_widget_demo_tab< Polyline_Traits >|      |  
|                 |--------------------------------------|      |
|                3|	Qt_widget_demo_tab< Conic_Traits >   |      |  
|                 |--------------------------------------|      |
|                4|	Qt_widget_demo_tab< Polyline_Traits >|      |  
|                 |--------------------------------------|      |
|                5|	Qt_widget_demo_tab< Segment_Traits > |      |  
|                 |--------------------------------------|      |
|                 |	...................................  |      |
|                 |	...................................  |      |
|                 |	...................................  |      |
|                 |	...................................  |      |
|                 |	                                     |      |  
|                 |--------------------------------------|      |
|                                                               |
|                                                               |
|---------------------------------------------------------------|


          CGAL::Qt_widget
                 |
                 v
        Qt_widget_base_tab
                 |
                 v
Qt_widget_demo_tab < class Traits >


we can divide the program operations into 5 categories:
1. operations of the main window: 
 - MyWindow.h: declerations and definitions.
 - MyWindow.C: constructor destructor and main.
 - MyWindow_files.C: file handling procedures.
 - MyWindow_overlay.C: overlay procedures.
 - MyWindow_operations.C: other operations.

2. forms and dialog boxes (forms.h, forms.C): implementation of the program
forms like properties and overlay forms.

3. Qt_widget_base_tab operations (demo_tab.h - class Qt_widget_base_tab): procedures 
with no concern of the actual traits.

4. Qt_widget_demo_tab < class Traits > operations (demo_tab.h - traits classes): 
special treatment in any trait class. 

5. notification operations (seg_notif.h, pol_notif.h, conic_notif.h): inheritance
of class Pmwx_change_notification in each traits type.


more files:
- cgal_types: type definitions
- qt_layer.h/qt_layer.C: inheritance of class Qt_widget_layer
- makefile

known bugs:
1. read polyline Planar map file fail because 
   operator >> in Arr_polyline_traits doesn't work.
2. writing to postscript.
3. insert range of conic curves (sweepline).
4. overlay of overlaped curves.

things to do:
1. insert Parabula and Hyperbula conics.
2. allow choosing the point location strategy.
