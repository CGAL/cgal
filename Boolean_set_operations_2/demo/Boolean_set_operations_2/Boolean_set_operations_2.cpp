// Copyright (c) 2012  Tel-Aviv University (Israel).
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
// WARRANTY OF describingIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s) : Apurva Bhatt <response2apurva@gmail.com>
//             Efi Fogel <efifogel@gmain.com>

//The demo contains no error handling

#include <fstream>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <list>

#include <boost/shared_ptr.hpp>

#include <QApplication>
#include <qmessagebox.h>

#include <QMainWindow>
#include <QGraphicsScene>
#include <QActionGroup>
#include <QPen>
#include <QtGui>
#include <QString>
#include <QFileDialog>
#include <QInputDialog>
#include <QDragEnterEvent>
#include <QDropEvent>
#include <QPainter>
#include <QSlider>
#include <QProgressBar>
#include <QMessageBox>

#include <CGAL/basic.h>
#include <CGAL/Point_2.h>
#include <CGAL/number_utils.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/Timer.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/iterator.h>
#include <CGAL/assertions_behaviour.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Lazy_exact_nt.h>
#include <CGAL/General_polygon_set_2.h>
#include <CGAL/CORE_algebraic_number_traits.h>
#include <CGAL/Arr_Bezier_curve_traits_2.h>
#include <CGAL/Arr_conic_traits_2.h>
#include <CGAL/Iso_rectangle_2.h>
#include <CGAL/Gps_traits_2.h>
#include <CGAL/minkowski_sum_2.h>
#include <CGAL/approximated_offset_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Boolean_set_operations_2.h>

#ifdef CGAL_USE_GMP
  #include <CGAL/Gmpq.h>
#else
  #include <CGAL/MP_Float.h>
  #include <CGAL/Quotient.h>
#endif

#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Qt/Converter.h>
#include <CGAL/Qt/DemosMainWindow.h>
#include <CGAL/Qt/utility.h>
#include <CGAL/IO/Dxf_bsop_reader.h>
#include <CGAL/Qt/GraphicsViewNavigation.h>

#include "QT5/Circular_polygons.h"
#include "QT5/Linear_polygons.h"
#include "QT5/Graphics_view_circular_polygon_input.h"
#include "QT5/Graphics_view_linear_polygon_input.h"
#include "QT5/Graphics_view_linear_polygon_input.h"
#include "QT5/General_polygon_2.h"
#include "QT5/General_polygon_set_2.h"
#include "QT5/General_polygon_set_on_surface_2.h"
#include "QT5/Gps_circle_segment_traits_2.h"
#include "QT5/Gps_segment_traits_2.h"
#include "QT5/Gps_traits_2.h"
#include "QT5/connect_holes.h"
#include "QT5/Polygon_set_2.h"
#include "QT5/BezierCurves.h"
#include "QT5/GraphicsViewBezierPolygonInput.h"



#include "ui_Boolean_set_operations_2.h"

#include "Typedefs.h"

// Forward Declarations
void show_warning(std::string aS);
void show_error(std::string aS);
void error(std::string aS);
void error_handler(char const* what, char const* expr, char const* file,
                   int line, char const* msg);


//Circular_polygon linearPart_2_circ(Circular_Linear_polygon const& pgn);
//Circular_polygon_with_holes linearPart_2_circ(Circular_Linear_polygon_with_holes const& pwh);
//bool read_linear(QString aFileName, Circular_polygon_set& rSet,
                   //Circular_region_source_container& rSources);
//bool read_dxf ( QString aFileName, Circular_polygon_set& rSet, 
                //    Circular_region_source_container& rSources );
//bool read_bezier ( QString aFileName, Bezier_polygon_set& rSet,
     //               Bezier_region_source_container& rSources  );
//Bezier_curve read_bezier_curve ( std::istream& is, bool aDoubleFormat );

// Typedefs
typedef CGAL::Qt::Circular_set_graphics_item<Circular_polygon_set,
                                             Circular_traits>   Circular_GI;
typedef CGAL::Qt::Linear_set_graphics_item<Linear_polygon_set, Linear_traits>
                                                                Linear_GI;

// Functions to show errors
void show_warning(std::string aS)
{ QMessageBox::warning(NULL, "Warning", QString(aS.c_str())); }

void show_error(std::string aS)
{ QMessageBox::critical(NULL, "Critical Error", QString(aS.c_str())); }

void error(std::string aS)
{
  show_error(aS);
  throw std::runtime_error(aS);
}

void error_handler(char const* what, char const* expr, char const* file,
                   int line, char const* msg)
{
  std::ostringstream ss;

  ss << "CGAL error: " << what << " violation!" << std::endl
     << "Expr: " << expr << std::endl
     << "File: " << file << std::endl
     << "Line: " << line << std::endl;
  if (msg != 0)
    ss << "Explanation:" << msg << std::endl;

  error(ss.str());
}

//****************************************************

//A way to maintain 3 set of polygons namely red,blue and result for all
// boolean operations

enum {
  BLUE_GROUP, RED_GROUP, BLACK_GROUP, BROWN_GROUP, YELLOW_GROUP,
  MAGENTA_GROUP, AQUA_GROUP, RESULT_GROUP , UNIVERSAL_GROUP
};

//A way to maintain 3 category of polygons namely linear,circular
//enum genrates errors so, we wil use LINEAR_TYPE=1, CIRCULAR_TYPE=2and BEZIER_TPYE = 3 ConicCircle = 4 and ConicEllipse = 5
//enum { LINEAR_TYPE, CIRCULAR_TYPE, BEZIER_TPYE,ConicCircle,ConicEllipse};

//dawing tools
QPen sPens[] = {
  QPen(QColor(0,0,255),0,Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin),   //blue
  QPen(QColor(255,0,0),0,Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin),   //red
  QPen(QColor(0,0,0),0,Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin),     //black
  QPen(QColor(210,105,30),0,Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin),//brown
  QPen(QColor(255,255,0),0,Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin), //yellow
  QPen(QColor(255,0,255),0,Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin), //magenta
  QPen(QColor(0,255,255),0,Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin), //aqua
  QPen(QColor(0,255,0),0,Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin),    //green
  QPen(Qt::gray,0,Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin)    //gray
};


QBrush sBrushes[] = {
  QBrush(QColor(0,0,255,32)),           //blue
  QBrush(QColor(255,0,0,32)),           //red
  QBrush(QColor(0,0,0,32)),             //black
  QBrush(QColor(210,105,30,32)),        //brown
  QBrush(QColor(255,255,0,32)),         //yellow
  QBrush(QColor(255,0,255,32)),         //magenta
  QBrush(QColor(0,255,255,32)),         //aqua
  QBrush(QColor(0,255,0,220)),           //green(reserved for result)
  QBrush(QColor(255,255,255,32))		//gray univerasal
};
//**************************************

//A base call for rep class
struct Rep_base {
  virtual ~Rep_base() {}

  virtual int type () const = 0;

  virtual CGAL::Qt::GraphicsItem* gi() const = 0;
  virtual CGAL::Qt::GraphicsItem* gi() = 0;

  virtual void set_pen(QPen const& aPen) = 0;
  virtual void set_brush(QBrush const& aBrush) = 0;

  virtual QRectF bounding_rect() const { return gi()->boundingRect(); }

  virtual bool is_empty() const = 0;

  virtual void clear() = 0;
  virtual void complement() = 0;
  virtual void assign(Rep_base const& aOther) = 0;
  virtual void intersect(Rep_base const& aOther) = 0;
  virtual void join(Rep_base const& aOther) = 0;
  virtual void difference(Rep_base const& aOther) = 0;
  virtual void symmetric_difference(Rep_base const& aOther) = 0;
};

//Class for initializing
template <typename GI_, typename Set_, typename Gps_traits>
class Rep : public Rep_base {
public:
  typedef GI_  GI;
  typedef Set_ Set;
  typedef Rep<GI,Set,Gps_traits> Self;
  typedef Gps_traits m_tratis;

  Rep() { m_GI = new GI(&m_set,m_traits); }

  Set const& set() const { return m_set; }
  Set & set() { return m_set; }

  virtual CGAL::Qt::GraphicsItem* gi() const { return m_GI; }
  virtual CGAL::Qt::GraphicsItem* gi()       { return m_GI; }

  virtual void set_pen(QPen const& aPen) { m_GI->setPen  (aPen);   }
  virtual void set_brush(QBrush const& aBrush) { m_GI->setBrush(aBrush); }

  virtual bool is_empty() const { return m_set.is_empty(); }

  virtual void clear()
  {
    try
    {
      m_set.clear();
    }
    catch(...)
    {
      show_error("Exception thrown during boolean operation clear");
    }
  }

  virtual void complement()
  {
    try
    {
      m_set.complement();
    }
    catch(...)
    {
      show_error("Exception thrown during boolean operation complement");
    }
  }

  virtual void assign(Rep_base const& aOther)
  {
    try
    {
      m_set = cast(aOther).m_set;
    }
    catch(...)
    {
      show_error("Exception thrown during boolean operation assign");
    }
  }

  virtual void intersect(Rep_base const& aOther)
  {
    try
    {
      m_set.intersection(cast(aOther).m_set);
    }
    catch(...)
    {
      show_error("Exception thrown during boolean operation intersect");
    }
  }

  virtual void join(Rep_base const& aOther)
  {
    try
    {
      m_set.join(cast(aOther).m_set);
    }
    catch(...)
    {
      show_error("Exception thrown during boolean operation union");
    }
  }

  virtual void difference(Rep_base const& aOther)
  {
    try
    {
      m_set.difference(cast(aOther).m_set);
    }
    catch(...)
    {
      show_error("Exception thrown during boolean operation difference");
    }
  }

  virtual void symmetric_difference(Rep_base const& aOther)
  {
    try
    {
      m_set.symmetric_difference(cast(aOther).m_set);
    }
    catch(...)
    {
      show_error("Exception thrown during boolean operation symmetric difference");
    }
  }

  static Self const& cast(Rep_base const& aOther)
  { return dynamic_cast<Self const&>(aOther); }

  static Self& cast(Rep_base& aOther) { return dynamic_cast<Self&>(aOther); }

private:
  //For maintaining all drawing operations
  GI* m_GI;
  //Storage for all polygons of one type. It is used as a base to perform all boolean operations
  Set m_set;
protected:
  //pass it
  Gps_traits m_traits;
};



//Implementing Bezier's rep class
template<class GI_, class Set_>
class Rep_o : public Rep_base
{
public:

  typedef GI_  GI  ;
  typedef Set_ Set ;
  
  typedef Rep_o<GI,Set> Self ;
  
  Rep_o() { mGI = new GI(&mSet) ; }
  
  Set const& set() const { return mSet ; }
  Set      & set()       { return mSet ; }
  
  virtual CGAL::Qt::GraphicsItem* gi() const { return mGI; }
  virtual CGAL::Qt::GraphicsItem* gi()       { return mGI; }
  
  virtual void set_pen  ( QPen   const& aPen   ) { mGI->setPen  (aPen);   } 
  virtual void set_brush( QBrush const& aBrush ) { mGI->setBrush(aBrush); }
      
  virtual bool is_empty() const { return mSet.is_empty() ; }
  
  virtual void clear()                         
  { 
    try
    {
      mSet.clear() ; 
    } 
    catch(...)
    {
      show_error("Exception thrown during boolean operation");
    } 
  }
  
  virtual void complement()                         
  { 
    try
    {
      mSet.complement(); 
    } 
    catch(...)
    {
      show_error("Exception thrown during boolean operation");
    } 
  }
  
  virtual void assign( Rep_base const& aOther ) 
  { 
    try
    {
      mSet = cast(aOther).mSet; 
    } 
    catch(...)
    {
      show_error("Exception thrown during boolean operation");
    } 
  }
  
  virtual void intersect( Rep_base const& aOther ) 
  { 
    try
    {
      mSet.intersection( cast(aOther).mSet); 
    } 
    catch(...)
    {
      show_error("Exception thrown during boolean operation");
    } 
  }
  
  virtual void join( Rep_base const& aOther ) 
  { 
    try
    {
      mSet.join( cast(aOther).mSet); 
    } 
    catch(...)
    {
      show_error("Exception thrown during boolean operation");
    } 
  }
  
  virtual void difference( Rep_base const& aOther ) 
  { 
    try
    {
      mSet.difference( cast(aOther).mSet); 
    } 
    catch(...)
    {
      show_error("Exception thrown during boolean operation");
    } 
  }
  
  virtual void symmetric_difference( Rep_base const& aOther ) 
  { 
    try
    {
      mSet.symmetric_difference( cast(aOther).mSet); 
    } 
    catch(...)
    {
      show_error("Exception thrown during boolean operation");
    } 
  }
  
  static Self const& cast( Rep_base const& aOther ) { return dynamic_cast<Self const&>(aOther); }
  static Self      & cast( Rep_base      & aOther ) { return dynamic_cast<Self      &>(aOther); }
  
private:

  GI* mGI;
  Set mSet ;
} ;


//A class for connecting GUI and this file
class Circular_rep : public Rep<Circular_GI, Circular_polygon_set,
                                Circular_traits>
{
  typedef Rep<Circular_GI, Circular_polygon_set,Circular_traits> Base;

public:
  Circular_rep () : Base() {}

  virtual int type() const { return 2; }
};

class Bezier_rep : public Rep_o<Bezier_GI, Bezier_polygon_set>
{
  typedef Rep_o<Bezier_GI, Bezier_polygon_set> Base ;
  
public:
  
  
  Bezier_rep () : Base() {} 
  
  virtual int type() const { return 3; }
} ;


//A class for connecting GUI and this file
class Linear_rep : public Rep<Linear_GI, Linear_polygon_set, Linear_traits>
{

typedef Rep<Linear_GI, Linear_polygon_set, Linear_traits> Base;
public:
  Linear_rep () : Base() {}

  virtual int type() const { return 1; }
};

class Curve_set {
  // a conatiner which deletes an object when last shared_ptr gets deleted or
  // re-initiated
  typedef boost::shared_ptr<Rep_base> Rep_ptr;

public:
  //constructor
  Curve_set(int aType, QPen aPen, QBrush aBrush) : m_pen(aPen), m_brush(aBrush)
  { reset_type(aType); }

  void reset_type(int aType)
  {
    //setting shared_ptr for repective polygon
    if (aType == 1) m_rep = Rep_ptr(new Linear_rep());
    else if(aType == 2) m_rep=Rep_ptr(new Circular_rep());
    else if(aType == 3) m_rep=Rep_ptr(new Bezier_rep());
    //setting pen and brush
    m_rep->set_pen  (m_pen);
    m_rep->set_brush(m_brush);
  }

  
  CGAL::Qt::GraphicsItem const* gi() const { return m_rep->gi(); }
  CGAL::Qt::GraphicsItem* gi() { return m_rep->gi(); }

  QRectF bounding_rect() const { return m_rep->bounding_rect(); }

  bool is_empty() const { return !m_rep || m_rep->is_empty(); }

  void clear() { m_rep->clear(); }
  //boolean operations
  void complement () { m_rep->complement(); }

  void assign (Curve_set const& aOther)
  {
    if (is_circular() && aOther.is_circular())
      get_circular_rep()->assign(*aOther.get_circular_rep());
    else if (is_bezier() && aOther.is_bezier()) 
      get_bezier_rep()->assign( *aOther.get_bezier_rep()) ;
    else if (is_linear() && aOther.is_linear())
      get_linear_rep()->assign(*aOther.get_linear_rep());
  }

  void intersect(Curve_set const& aOther)
  {
    if (is_circular() && aOther.is_circular())
      get_circular_rep()->intersect(*aOther.get_circular_rep());
    else if ( is_bezier() && aOther.is_bezier() )
      get_bezier_rep()->intersect( *aOther.get_bezier_rep() ) ;
    else if (is_linear() && aOther.is_linear())
      get_linear_rep()->intersect(*aOther.get_linear_rep());
  }

  void join(Curve_set const& aOther)
  {
    if (is_circular() && aOther.is_circular())
      get_circular_rep()->join(*aOther.get_circular_rep());
    else if ( is_bezier() && aOther.is_bezier() )
      get_bezier_rep()->join( *aOther.get_bezier_rep() ) ;
    else if (is_linear() && aOther.is_linear())
      get_linear_rep()->join(*aOther.get_linear_rep());
  }

  void difference(Curve_set const& aOther)
  {
    if (is_circular() && aOther.is_circular())
      get_circular_rep()->difference(*aOther.get_circular_rep());
    else if ( is_bezier() && aOther.is_bezier() )
      get_bezier_rep()->difference( *aOther.get_bezier_rep() ) ;
    else if (is_linear() && aOther.is_linear())
      get_linear_rep()->difference(*aOther.get_linear_rep());
  }

  void symmetric_difference(Curve_set const& aOther)
  {
    if (is_circular() && aOther.is_circular())
      get_circular_rep()->symmetric_difference(*aOther.get_circular_rep());
    else if ( is_bezier() && aOther.is_bezier() )
      get_bezier_rep()->symmetric_difference( *aOther.get_bezier_rep() ) ;
    else if (is_linear() && aOther.is_linear())
      get_linear_rep()->symmetric_difference(*aOther.get_linear_rep());
  }

  //see its need keep it for now
  const Rep_base& rep() const { return *m_rep; }
  Rep_base& rep() { return *m_rep; }

  bool is_circular() const { return m_rep->type() == 2; }
  bool is_bezier  () const { return m_rep->type() == 3; }
  bool is_linear() const { return m_rep->type() == 1; }

  //to get rep for circualr polygons
  const Circular_rep* get_circular_rep() const
  { return dynamic_cast<Circular_rep const*>(boost::get_pointer(m_rep)); }

  Circular_rep* get_circular_rep()
  { return dynamic_cast<Circular_rep*  >(boost::get_pointer(m_rep)); }

  //to get Circular_polygon_set
  const Circular_polygon_set& circular() const
  { return get_circular_rep()->set(); }

  Circular_polygon_set& circular() { return get_circular_rep()->set(); }


  const Bezier_rep* get_bezier_rep() const 
  { return dynamic_cast<Bezier_rep   const*>( boost::get_pointer(m_rep) ); }

  Bezier_rep* get_bezier_rep()       
  { return dynamic_cast<Bezier_rep*  >( boost::get_pointer(m_rep) ); }
  
  const Bezier_polygon_set& bezier() const 
  { return get_bezier_rep()->set(); }

  Bezier_polygon_set& bezier() { return get_bezier_rep ()->set(); }


   //to get rep for linear polygons
  const Linear_rep* get_linear_rep() const
  { return dynamic_cast<Linear_rep const*>(boost::get_pointer(m_rep)); }

  Linear_rep* get_linear_rep()
  { return dynamic_cast<Linear_rep*>(boost::get_pointer(m_rep)); }

  //to get Linear_polygon_set
  const Linear_polygon_set& linear() const { return get_linear_rep()->set(); }
  Linear_polygon_set& linear() { return get_linear_rep()->set(); }

private:
  //drawing tools
  QPen m_pen;
  QBrush m_brush;
  // a conatiner which deletes an object when last shared_ptr gets deleted or
  // re-initiated
  boost::shared_ptr<Rep_base> m_rep;
};

typedef std::vector<Curve_set>                  Curve_set_container;

typedef Curve_set_container::const_iterator     Curve_set_const_iterator;
typedef Curve_set_container::iterator           Curve_set_iterator;


class MainWindow : public CGAL::Qt::DemosMainWindow ,
                   public Ui::Boolean_set_operations_2
{
  Q_OBJECT// removing it gives error ui not declared

private:
  QGraphicsScene m_scene;
  //keep it intact for now check it out
  bool m_circular_active;
  bool m_bezier_active;
  //which type is currently active now
  //bool m_blue_active;
  size_t m_color_active;
  size_t m_color_complement; //allows only one color {0:blue,1:red,2:black,3:brown,4:yellow,5:magenta,6:aqua}
  bool m_blue_int;
  bool m_red_int;
  bool m_black_int;
  bool m_brown_int;
  bool m_yellow_int;
  bool m_magenta_int;
  bool m_aqua_int;

  bool m_blue_union;
  bool m_red_union;
  bool m_black_union;
  bool m_brown_union;
  bool m_yellow_union;
  bool m_magenta_union;
  bool m_aqua_union;

  size_t m_color_diff_A;
  size_t m_color_diff_B;

  bool m_blue_sym_diff;
  bool m_red_sym_diff;
  bool m_black_sym_diff;
  bool m_brown_sym_diff;
  bool m_yellow_sym_diff;
  bool m_magenta_sym_diff;
  bool m_aqua_sym_diff;
  

  bool m_visible_black;
  bool m_visible_brown;
  bool m_visible_yellow;
  bool m_visible_magenta;
  bool m_visible_aqua;

  Curve_set_container m_curve_sets;
  //container for curves
  Circular_region_source_container m_blue_circular_sources;
  Circular_region_source_container m_red_circular_sources;
  Circular_region_source_container m_black_circular_sources;
  Circular_region_source_container m_brown_circular_sources;
  Circular_region_source_container m_yellow_circular_sources;
  Circular_region_source_container m_magenta_circular_sources;
  Circular_region_source_container m_aqua_circular_sources;
  Circular_region_source_container m_univ_circular_sources;

  Linear_region_source_container m_blue_linear_sources;
  Linear_region_source_container m_red_linear_sources;
  Linear_region_source_container m_black_linear_sources;
  Linear_region_source_container m_brown_linear_sources;
  Linear_region_source_container m_yellow_linear_sources;
  Linear_region_source_container m_magenta_linear_sources;
  Linear_region_source_container m_aqua_linear_sources;
  Linear_region_source_container m_univ_linear_sources;

  Bezier_region_source_container m_blue_bezier_sources;
  Bezier_region_source_container m_red_bezier_sources;
  Bezier_region_source_container m_black_bezier_sources;
  Bezier_region_source_container m_brown_bezier_sources;
  Bezier_region_source_container m_yellow_bezier_sources;
  Bezier_region_source_container m_magenta_bezier_sources;
  Bezier_region_source_container m_aqua_bezier_sources;
  Bezier_region_source_container m_univ_bezier_sources;

  //typedefs of classes used to draw circular and linear polygon
  CGAL::Qt::Graphics_view_linear_polygon_input<Kernel>* m_linear_input;
  CGAL::Qt::Graphics_view_circular_polygon_input<Kernel>* m_circular_input;
  CGAL::Qt::GraphicsViewBezierPolygonInput<Bezier_traits>* m_bezier_input ;

public:
  MainWindow();

private:
  void dragEnterEvent(QDragEnterEvent* event);
  void dropEvent(QDropEvent* event);
  void zoomToFit();
  void ToogleView(size_t aGROUP, bool a_check);

protected slots:
  void open(QString filename);//for file handling


public slots:
  void processInput(CGAL::Object o);
  void on_actionNew_triggered();
  void on_actionRecenter_triggered();
  void on_actionComplement_triggered();
  void on_actionUnion_triggered();
  void on_actionIntersection_triggered();
  void on_actionDifference_triggered();
  void on_actionSymmetric_Difference_triggered();
  void on_actionMinkowski_Sum_triggered();
  void on_actionInsertLinear_toggled(bool aChecked);
  void on_actionInsertCircular_toggled(bool aChecked);
  void on_actionInsertBezier_toggled(bool aChecked);
  void on_actionInsertConicCircle_triggered();
  void on_actionInsertConicEclipse_triggered();
  void on_actionOpenLinear_triggered();
  void on_actionOpenDXF_triggered();
  void on_actionOpenBezier_triggered();
  void on_actionSaveResult_triggered();
  void on_actionAddColor_triggered();
  void on_actionMinusColor_triggered();

  void on_showBlue_toggled  (bool a_check);
  void on_showRed_toggled   (bool a_check);
  void on_showBlack_toggled   (bool a_check);
  void on_showBrown_toggled   (bool a_check);
  void on_showYellow_toggled   (bool a_check);
  void on_showMagenta_toggled   (bool a_check);
  void on_showAqua_toggled   (bool a_check);
  void on_showResult_toggled(bool a_check);

  void on_showBlueComp_toggled(bool aCheck);
  void on_showRedComp_toggled(bool aCheck);
  void on_showBlackComp_toggled(bool aCheck);
  void on_showBrownComp_toggled(bool aCheck);
  void on_showYellowComp_toggled(bool aCheck);
  void on_showMagentaComp_toggled(bool aCheck);
  void on_showAquaComp_toggled(bool aCheck);

  void on_showBlueInt_toggled(bool aCheck);
  void on_showRedInt_toggled(bool aCheck);
  void on_showBlackInt_toggled(bool aCheck);
  void on_showBrownInt_toggled(bool aCheck);
  void on_showYellowInt_toggled(bool aCheck);
  void on_showMagentaInt_toggled(bool aCheck);
  void on_showAquaInt_toggled(bool aCheck);

  void on_showBlueUnion_toggled(bool aCheck);
  void on_showRedUnion_toggled(bool aCheck);
  void on_showBlackUnion_toggled(bool aCheck);
  void on_showBrownUnion_toggled(bool aCheck);
  void on_showYellowUnion_toggled(bool aCheck);
  void on_showMagentaUnion_toggled(bool aCheck);
  void on_showAquaUnion_toggled(bool aCheck);

  void on_showBlueDiff_toggled(bool aCheck);
  void on_showRedDiff_toggled(bool aCheck);
  void on_showBlackDiff_toggled(bool aCheck);
  void on_showBrownDiff_toggled(bool aCheck);
  void on_showYellowDiff_toggled(bool aCheck);
  void on_showMagentaDiff_toggled(bool aCheck);
  void on_showAquaDiff_toggled(bool aCheck);


  void on_showBlueSym_Diff_toggled(bool aCheck);
  void on_showRedSym_Diff_toggled(bool aCheck);
  void on_showBlackSym_Diff_toggled(bool aCheck);
  void on_showBrownSym_Diff_toggled(bool aCheck);
  void on_showYellowSym_Diff_toggled(bool aCheck);
  void on_showMagentaSym_Diff_toggled(bool aCheck);
  void on_showAquaSym_Diff_toggled(bool aCheck);

  void on_showBlueMink_Sum_toggled(bool aCheck);
  void on_showRedMink_Sum_toggled(bool aCheck);
  void on_showBlackMink_Sum_toggled(bool aCheck);
  void on_showBrownMink_Sum_toggled(bool aCheck);
  void on_showYellowMink_Sum_toggled(bool aCheck);
  void on_showMagentaMink_Sum_toggled(bool aCheck);
  void on_showAquaMink_Sum_toggled(bool aCheck);

  void on_drawBlue_toggled(bool a_check);
  void on_drawRed_toggled (bool a_check);
  void on_drawBlack_toggled (bool a_check);
  void on_drawBrown_toggled (bool a_check);
  void on_drawYellow_toggled (bool a_check);
  void on_drawMagenta_toggled (bool a_check);
  void on_drawAqua_toggled (bool a_check);
  void on_drawUniversalPolygon_toggled(bool a_check);


  //void on_actionAdd_new_polygon_triggered();
  void on_actionDeleteResult_triggered();
  void on_actionDelete_triggered();
  void on_actionDeleteAll_triggered();
  void on_actionPAN_toggled(bool aChecked);

signals:
  void changed();

private:
  void modelChanged() { emit(changed()); }

  //warning message for user
  bool ask_user_yesno(const char* aTitle, const char* aQuestion)
  {
    return QMessageBox::warning(this, aTitle, QString(aQuestion), "&Yes", "&No",
                                QString::null, 1, 1) == 0;
  }

  // for setting Curve_set of aGroup type an int representing a set of polygon
  // of a specific type
  Curve_set& set(size_t aGroup) { return m_curve_sets[aGroup]; }

  //setting curve
  Curve_set& blue_set() { return set(BLUE_GROUP); }
  Curve_set& red_set() { return set(RED_GROUP); }
  Curve_set& black_set() { return set(BLACK_GROUP); }
  Curve_set& brown_set() { return set(BROWN_GROUP); }
  Curve_set& yellow_set() { return set(YELLOW_GROUP); }
  Curve_set& magenta_set() { return set(MAGENTA_GROUP); }
  Curve_set& aqua_set() { return set(AQUA_GROUP); }
  Curve_set& result_set() { return set(RESULT_GROUP); }
  Curve_set& universal_set() {return set(UNIVERSAL_GROUP);}

  //gets which group is currently active now
  size_t active_group() const { return m_color_active; }
  // size_t complement_group() const {return m_color_complement; } //see if needed

  //sets the current active group
  Curve_set& active_set()  { return set(active_group()); }

  //returns circular containers
  // Circular_region_source_container const& blue_circular_sources() const
  // { return m_blue_circular_sources; }

  Circular_region_source_container& blue_circular_sources()
  { return m_blue_circular_sources; }

  // Circular_region_source_container const& red_circular_sources () const
  // { return m_red_circular_sources; }

  Circular_region_source_container& red_circular_sources()
  { return m_red_circular_sources; }

  // Circular_region_source_container const& red_circular_sources () const
  // { return m_red_circular_sources; }

  Circular_region_source_container& black_circular_sources()
  { return m_black_circular_sources; }

  Circular_region_source_container& brown_circular_sources()
  { return m_brown_circular_sources; }

  Circular_region_source_container& yellow_circular_sources()
  { return m_yellow_circular_sources; }

  Circular_region_source_container& magenta_circular_sources()
  { return m_magenta_circular_sources; }

  Circular_region_source_container& aqua_circular_sources()
  { return m_aqua_circular_sources; }

  Circular_region_source_container& univ_circular_sources()
  { return m_univ_circular_sources; }

  //returns linear containers
  // Linear_region_source_container const& blue_linear_sources() const
  // { return m_blue_linear_sources; }
  Linear_region_source_container& blue_linear_sources()
  { return m_blue_linear_sources; }

  // Linear_region_source_container const& red_linear_sources () const
  // { return m_red_linear_sources; }

  Linear_region_source_container& red_linear_sources()
  { return m_red_linear_sources; }

  Linear_region_source_container& black_linear_sources()
  { return m_black_linear_sources; }

  Linear_region_source_container& brown_linear_sources()
  { return m_brown_linear_sources; }

  Linear_region_source_container& yellow_linear_sources()
  { return m_yellow_linear_sources; }

  Linear_region_source_container& magenta_linear_sources()
  { return m_magenta_linear_sources; }

  Linear_region_source_container& aqua_linear_sources()
  { return m_aqua_linear_sources; }

  Linear_region_source_container& univ_linear_sources()
  { return m_univ_linear_sources; }

//Same for Bezier 

  Bezier_region_source_container& blue_bezier_sources()
  { return m_blue_bezier_sources; }

  Bezier_region_source_container& red_bezier_sources()
  { return m_red_bezier_sources; }

  Bezier_region_source_container& black_bezier_sources()
  { return m_black_bezier_sources; }

  Bezier_region_source_container& brown_bezier_sources()
  { return m_brown_bezier_sources; }

  Bezier_region_source_container& yellow_bezier_sources()
  { return m_yellow_bezier_sources; }

  Bezier_region_source_container& magenta_bezier_sources()
  { return m_magenta_bezier_sources; }

  Bezier_region_source_container& aqua_bezier_sources()
  { return m_aqua_bezier_sources; }

  Bezier_region_source_container& univ_bezier_sources()
  { return m_univ_bezier_sources; }


  //returns active blue container
  // Circular_region_source_container const& active_circular_sources() const
  // { return m_blue_active ? m_blue_circular_sources : m_red_circular_sources; }



  Circular_region_source_container& active_circular_sources()
  {
    //return m_blue_active ? m_blue_circular_sources : m_red_circular_sources;
    switch(m_color_active) {
     case 0: return m_blue_circular_sources;
     case 1: return m_red_circular_sources;
     case 2: return m_black_circular_sources;
     case 3: return m_brown_circular_sources;
     case 4: return m_yellow_circular_sources;
     case 5: return m_magenta_circular_sources;
     case 6: return m_aqua_circular_sources;
     case 8: return m_univ_circular_sources;

     default: break;
    }

    CGAL_warning_msg(true, "Should not reach here!");
    return m_blue_circular_sources;
  }

  //returns active linear container
  //  Linear_region_source_container const& active_linear_sources() const
  // { return m_blue_active ? m_blue_linear_sources : m_red_linear_sources; }

  Linear_region_source_container& active_linear_sources()
  {
    //return m_blue_active ? m_blue_linear_sources : m_red_linear_sources;
    switch(m_color_active) {
     case 0: return m_blue_linear_sources;
     case 1: return m_red_linear_sources;
     case 2: return m_black_linear_sources;
     case 3: return m_brown_linear_sources;
     case 4: return m_yellow_linear_sources;
     case 5: return m_magenta_linear_sources;
     case 6: return m_aqua_linear_sources;
     case 8: return m_univ_linear_sources;

     default: break;
    }

    CGAL_warning_msg(true, "Should not reach here!");
    return m_blue_linear_sources;
  }

  Bezier_region_source_container& active_bezier_sources()
  {
    //return m_blue_active ? m_blue_linear_sources : m_red_linear_sources;
    switch(m_color_active) {
     case 0: return m_blue_bezier_sources;
     case 1: return m_red_bezier_sources;
     case 2: return m_black_bezier_sources;
     case 3: return m_brown_bezier_sources;
     case 4: return m_yellow_bezier_sources;
     case 5: return m_magenta_bezier_sources;
     case 6: return m_aqua_bezier_sources;
     case 8: return m_univ_bezier_sources;

     default: break;
    }

    CGAL_warning_msg(true, "Should not reach here!");
    return m_blue_bezier_sources;
  }

  void SetViewBlue(bool a_check) { showBlue->setChecked(a_check); }
  void SetViewRed(bool a_check) { showRed->setChecked(a_check); }
  void SetViewBlack(bool a_check) { showBlack->setChecked(a_check); }
  void SetViewBrown(bool a_check) { showBrown->setChecked(a_check); }
  void SetViewYellow(bool a_check) { showYellow->setChecked(a_check); }
  void SetViewMagenta(bool a_check) { showMagenta->setChecked(a_check); }
  void SetViewAqua(bool a_check) { showAqua->setChecked(a_check); }
  void SetViewResult(bool a_check) { showResult->setChecked(a_check); }

  //changes the set of polygons of a specific type
  //void ToogleView(size_t aGROUP, bool a_check);

  void link_GI(CGAL::Qt::GraphicsItem* aGI)
  {
    QObject::connect(this, SIGNAL(changed()), aGI, SLOT(modelChanged()));
    m_scene.addItem(aGI);
  }

  void unlink_GI(CGAL::Qt::GraphicsItem* aGI)
  {
    m_scene.removeItem(aGI);
    QObject::disconnect(this, SIGNAL(changed()), aGI, SLOT(modelChanged()));
  }

  void switch_set_type(Curve_set& aSet, int aType);

  void switch_sets_type(int aType);

  bool ensure_circular_mode();

  bool ensure_bezier_mode();

  bool ensure_linear_mode();//see if it is need

};

MainWindow::MainWindow() :
  DemosMainWindow(),
  m_bezier_active(false), //default
  m_circular_active(false), //default
  m_color_active(0),    //default
  m_color_complement(0), //default
  m_blue_int(true), //default
  m_red_int(false), //default
  m_black_int(false), //default
  m_brown_int(false), //default
  m_yellow_int(false), //default
  m_magenta_int(false), //default
  m_aqua_int(false), //default
  m_blue_union(true), //default
  m_red_union(false), //default
  m_black_union(false), //default
  m_brown_union(false), //default
  m_yellow_union(false), //default
  m_magenta_union(false), //default
  m_aqua_union(false), //default
  m_color_diff_A(0), //default
  m_color_diff_B(0), //default
  m_blue_sym_diff(true), //default
  m_red_sym_diff(false), //default
  m_black_sym_diff(false), //default
  m_brown_sym_diff(false), //default
  m_yellow_sym_diff(false), //default
  m_magenta_sym_diff(false), //default
  m_aqua_sym_diff(false), //default
  m_visible_black(false), //default 
  m_visible_brown(false), //default 
  m_visible_yellow(false), //default 
  m_visible_magenta(false), //default 
  m_visible_aqua(false) //default 

  // setting color rows to setVisible-> False 
  /*showBlack->setVisible(false), //default
  drawBlack->setVisible(false),  //default
  showBlackComp->setVisible(false), //default
  showBlackDiff->setVisible(false), //default
  showBlackUnion->setVisible(false), //default
  showBlackInt->setVisible(false), //default
  showBlackSym_Diff->setVisible(false), //default
  showBlackMink_Sum->setVisible(false), //default
  showBrown->setVisible(false), //default
  drawBrown->setVisible(false),  //default
  showBrownComp->setVisible(false), //default
  showBrownDiff->setVisible(false), //default
  showBrownUnion->setVisible(false), //default
  showBrownInt->setVisible(false), //default
  showBrownSym_Diff->setVisible(false), //default
  showBrownMink_Sum->setVisible(false), //default
  showYellow->setVisible(false), //default
  drawYellow->setVisible(false),  //default
  showYellowComp->setVisible(false), //default
  showYellowDiff->setVisible(false), //default
  showYellowUnion->setVisible(false), //default
  showYellowInt->setVisible(false), //default
  showYellowSym_Diff->setVisible(false), //default
  showYellowMink_Sum->setVisible(false), //default
  showMagenta->setVisible(false), //default
  drawMagenta->setVisible(false),  //default
  showMagentaComp->setVisible(false), //default
  showMagentaDiff->setVisible(false), //default
  showMagentaUnion->setVisible(false), //default
  showMagentaInt->setVisible(false), //default
  showMagentaSym_Diff->setVisible(false), //default
  showMagentaMink_Sum->setVisible(false), //default
  showAqua->setVisible(false), //default
  drawAqua->setVisible(false),  //default
  showAquaComp->setVisible(false), //default
  showAquaDiff->setVisible(false), //default
  showAquaUnion->setVisible(false), //default
  showAquaInt->setVisible(false), //default
  showAquaSym_Diff->setVisible(false), //default
  showAquaMink_Sum->setVisible(false), //default*/
{
  CGAL::set_error_handler  (error_handler);
  CGAL::set_warning_handler(error_handler);

  setupUi(this);	

  setAcceptDrops(true);
  //default setups with Linear
  m_curve_sets.push_back(Curve_set(1, sPens[BLUE_GROUP], sBrushes[BLUE_GROUP]));
  m_curve_sets.push_back(Curve_set(1, sPens[RED_GROUP], sBrushes[RED_GROUP]));
  m_curve_sets.push_back(Curve_set(1, sPens[BLACK_GROUP],
                                   sBrushes[BLACK_GROUP]));
  m_curve_sets.push_back(Curve_set(1, sPens[BROWN_GROUP],
                                   sBrushes[BROWN_GROUP]));
  m_curve_sets.push_back(Curve_set(1, sPens[YELLOW_GROUP],
                                   sBrushes[YELLOW_GROUP]));
  m_curve_sets.push_back(Curve_set(1, sPens[MAGENTA_GROUP],
                                   sBrushes[MAGENTA_GROUP]));
  m_curve_sets.push_back(Curve_set(1, sPens[AQUA_GROUP], sBrushes[AQUA_GROUP]));
  m_curve_sets.push_back(Curve_set(1, sPens[RESULT_GROUP],
                                   sBrushes[RESULT_GROUP]));
  m_curve_sets.push_back(Curve_set(1, sPens[UNIVERSAL_GROUP],
                                   sBrushes[UNIVERSAL_GROUP]));

  for (auto si = m_curve_sets.begin(); si != m_curve_sets.end(); ++ si)
    link_GI(si->gi());
  //
  // Setup the m_scene and the view
  //
  m_scene.setItemIndexMethod(QGraphicsScene::NoIndex);
  m_scene.setSceneRect(-100, -100, 100, 100);
  this->graphicsView->setScene(&m_scene);
  this->graphicsView->setMouseTracking(true);
  //this->on_actionInsertLinear_triggered();

  // Turn the vertical axis upside down
  this->graphicsView->scale(1, -1);

  //adding basic setups

  //need to be given finishing touch
  // The navigation adds zooming and translation functionality to the
  // QGraphicsView
  this->addNavigation(this->graphicsView);
  //setting the menus
  this->setupStatusBar();
  this->setupOptionsMenu();
  //link to a page describing
  //link for about page of CGAL
  this->addAboutCGAL();

  this->addRecentFiles(this->menuFile, this->actionClose);

  //initializing classes to draw respective polygons using mouse
  m_bezier_input = new CGAL::Qt::GraphicsViewBezierPolygonInput<Bezier_traits>(this, &m_scene);
  m_linear_input = new CGAL::Qt::Graphics_view_linear_polygon_input<Kernel>(this, &m_scene);
  m_circular_input = new CGAL::Qt::Graphics_view_circular_polygon_input<Kernel>(this, &m_scene);

  //connecting GUI and the code base
  QObject::connect(m_linear_input, SIGNAL(generate(CGAL::Object)), this,
                   SLOT(processInput(CGAL::Object)));
  QObject::connect(m_circular_input, SIGNAL(generate(CGAL::Object)), this,
                   SLOT(processInput(CGAL::Object)));
  QObject::connect(m_bezier_input, SIGNAL(generate(CGAL::Object)), this,
                   SLOT(processInput(CGAL::Object)));
  
  m_scene.installEventFilter(m_linear_input);

  QObject::connect(this->actionQuit, SIGNAL(triggered()), this,
                    SLOT(close()));
  QObject::connect(this, SIGNAL(openRecentFile(QString)), this,
                    SLOT(open(QString)));//for file handling
  QObject::connect(drawBlue, SIGNAL(toggled(bool)), this,
                   SLOT(on_drawBlue_toggled (bool)));
  QObject::connect(drawRed , SIGNAL(toggled(bool)), this,
                   SLOT(on_drawRed_toggled(bool)));
  QObject::connect(drawBlack , SIGNAL(toggled(bool)), this,
                   SLOT(on_drawBlack_toggled(bool)));
  QObject::connect(drawBrown , SIGNAL(toggled(bool)), this,
                   SLOT(on_drawBrown_toggled(bool)));
  QObject::connect(drawYellow , SIGNAL(toggled(bool)), this,
                   SLOT(on_drawYellow_toggled(bool)));
  QObject::connect(drawMagenta , SIGNAL(toggled(bool)), this,
                   SLOT(on_drawMagenta_toggled(bool)));
  QObject::connect(drawAqua , SIGNAL(toggled(bool)), this,
                   SLOT(on_drawAqua_toggled(bool)));
  QObject::connect(drawUniversalPolygon , SIGNAL(toggled(bool)), this,
                   SLOT(on_drawUniversalPolygon_toggled(bool)));
  QObject::connect(actionPAN,SIGNAL(triggered()), this,
  				   SLOT(on_actionPAN_toggled(bool)));
  QObject::connect(showBlue, SIGNAL(toggled(bool)), this,
                   SLOT(on_showBlue_toggled(bool)));
  QObject::connect(showRed, SIGNAL(toggled(bool)), this,
                   SLOT(on_showRed_toggled(bool)));
  QObject::connect(showBlack, SIGNAL(toggled(bool)), this,
                   SLOT(on_showBlack_toggled(bool)));
  QObject::connect(showBrown, SIGNAL(toggled(bool)), this,
                   SLOT(on_showBrown_toggled(bool)));
  QObject::connect(showYellow, SIGNAL(toggled(bool)), this,
                   SLOT(on_showYellow_toggled(bool)));
  QObject::connect(showMagenta, SIGNAL(toggled(bool)), this,
                   SLOT(on_showMagenta_toggled(bool)));
  QObject::connect(showAqua, SIGNAL(toggled(bool)), this,
                   SLOT(on_showAqua_toggled(bool)));
  QObject::connect(showResult, SIGNAL(toggled(bool)), this,
                   SLOT(on_showResult_toggled(bool)));


  //for the colour buckets

  //complement
  QObject::connect(showBlueComp, SIGNAL(toggled(bool)), this,
                   SLOT(on_showBlueComp_toggled(bool)));
  QObject::connect(showRedComp, SIGNAL(toggled(bool)), this,
                   SLOT(on_showRedComp_toggled(bool)));
  QObject::connect(showBlackComp, SIGNAL(toggled(bool)), this,
                   SLOT(on_showBlackComp_toggled(bool)));
  QObject::connect(showBrownComp, SIGNAL(toggled(bool)), this,
                   SLOT(on_showBrownComp_toggled(bool)));
  QObject::connect(showYellowComp, SIGNAL(toggled(bool)), this,
                   SLOT(on_showYellowComp_toggled(bool)));
  QObject::connect(showMagentaComp, SIGNAL(toggled(bool)), this,
                   SLOT(on_showMagentaComp_toggled(bool)));
  QObject::connect(showAquaComp, SIGNAL(toggled(bool)), this,
                   SLOT(on_showAquaComp_toggled(bool)));

  //intersection

  QObject::connect(showBlueInt, SIGNAL(toggled(bool)), this,
                   SLOT(on_showBlueInt_toggled(bool)));
  QObject::connect(showRedInt, SIGNAL(toggled(bool)), this,
                   SLOT(on_showRedInt_toggled(bool)));
  QObject::connect(showBlackInt, SIGNAL(toggled(bool)), this,
                   SLOT(on_showBlackInt_toggled(bool)));
  QObject::connect(showBrownInt, SIGNAL(toggled(bool)), this,
                   SLOT(on_showBrownInt_toggled(bool)));
  QObject::connect(showYellowInt, SIGNAL(toggled(bool)), this,
                   SLOT(on_showYellowInt_toggled(bool)));
  QObject::connect(showMagentaInt, SIGNAL(toggled(bool)), this,
                   SLOT(on_showMagentaInt_toggled(bool)));
  QObject::connect(showAquaInt, SIGNAL(toggled(bool)), this,
                   SLOT(on_showAquaInt_toggled(bool)));

  //union

  QObject::connect(showBlueUnion, SIGNAL(toggled(bool)), this,
                   SLOT(on_showBlueUnion_toggled(bool)));
  QObject::connect(showRedUnion, SIGNAL(toggled(bool)), this,
                   SLOT(on_showRedUnion_toggled(bool)));
  QObject::connect(showBlackUnion, SIGNAL(toggled(bool)), this,
                   SLOT(on_showBlackUnion_toggled(bool)));
  QObject::connect(showBrownUnion, SIGNAL(toggled(bool)), this,
                   SLOT(on_showBrownUnion_toggled(bool)));
  QObject::connect(showYellowUnion, SIGNAL(toggled(bool)), this,
                   SLOT(on_showYellowUnion_toggled(bool)));
  QObject::connect(showMagentaUnion, SIGNAL(toggled(bool)), this,
                   SLOT(on_showMagentaUnion_toggled(bool)));
  QObject::connect(showAquaUnion, SIGNAL(toggled(bool)), this,
                   SLOT(on_showAquaUnion_toggled(bool)));

  //difference

  QObject::connect(showBlueDiff, SIGNAL(toggled(bool)), this,
                   SLOT(on_showBlueDiff_toggled(bool)));
  QObject::connect(showRedDiff, SIGNAL(toggled(bool)), this,
                   SLOT(on_showRedDiff_toggled(bool)));
  QObject::connect(showBlackDiff, SIGNAL(toggled(bool)), this,
                   SLOT(on_showBlackDiff_toggled(bool)));
  QObject::connect(showBrownDiff, SIGNAL(toggled(bool)), this,
                   SLOT(on_showBrownDiff_toggled(bool)));
  QObject::connect(showYellowDiff, SIGNAL(toggled(bool)), this,
                   SLOT(on_showYellowDiff_toggled(bool)));
  QObject::connect(showMagentaDiff, SIGNAL(toggled(bool)), this,
                   SLOT(on_showMagentaDiff_toggled(bool)));
  QObject::connect(showAquaDiff, SIGNAL(toggled(bool)), this,
                   SLOT(on_showAquaDiff_toggled(bool)));

  //sym_difference

  QObject::connect(showBlueSym_Diff, SIGNAL(toggled(bool)), this,
                   SLOT(on_showBlueSym_Diff_toggled(bool)));
  QObject::connect(showRedSym_Diff, SIGNAL(toggled(bool)), this,
                   SLOT(on_showRedSym_Diff_toggled(bool)));
  QObject::connect(showBlackSym_Diff, SIGNAL(toggled(bool)), this,
                   SLOT(on_showBlackSym_Diff_toggled(bool)));
  QObject::connect(showBrownSym_Diff, SIGNAL(toggled(bool)), this,
                   SLOT(on_showBrownSym_Diff_toggled(bool)));
  QObject::connect(showYellowSym_Diff, SIGNAL(toggled(bool)), this,
                   SLOT(on_showYellowSym_Diff_toggled(bool)));
  QObject::connect(showMagentaSym_Diff, SIGNAL(toggled(bool)), this,
                   SLOT(on_showMagentaSym_Diff_toggled(bool)));
  QObject::connect(showAquaSym_Diff, SIGNAL(toggled(bool)), this,
                   SLOT(on_showAquaSym_Diff_toggled(bool)));

  //minkowski_sum

  QObject::connect(showBlueMink_Sum, SIGNAL(toggled(bool)), this,
                   SLOT(on_showBlueMink_Sum_toggled(bool)));
  QObject::connect(showRedMink_Sum, SIGNAL(toggled(bool)), this,
                   SLOT(on_showRedMink_Sum_toggled(bool)));
  QObject::connect(showBlackMink_Sum, SIGNAL(toggled(bool)), this,
                   SLOT(on_showBlackMink_Sum_toggled(bool)));
  QObject::connect(showBrownMink_Sum, SIGNAL(toggled(bool)), this,
                   SLOT(on_showBrownMink_Sum_toggled(bool)));
  QObject::connect(showYellowMink_Sum, SIGNAL(toggled(bool)), this,
                   SLOT(on_showYellowMink_Sum_toggled(bool)));
  QObject::connect(showMagentaMink_Sum, SIGNAL(toggled(bool)), this,
                   SLOT(on_showMagentaMink_Sum_toggled(bool)));
  QObject::connect(showAquaMink_Sum, SIGNAL(toggled(bool)), this,
                   SLOT(on_showAquaMink_Sum_toggled(bool)));

}

void MainWindow::on_showBlue_toggled(bool a_check)
{ ToogleView(BLUE_GROUP, a_check); }

void MainWindow::on_showRed_toggled(bool a_check)
{ ToogleView(RED_GROUP, a_check); }

void MainWindow::on_showBlack_toggled(bool a_check)
{ ToogleView(BLACK_GROUP, a_check); }

void MainWindow::on_showBrown_toggled(bool a_check)
{ ToogleView(BROWN_GROUP, a_check); }

void MainWindow::on_showYellow_toggled(bool a_check)
{ ToogleView(YELLOW_GROUP, a_check); }

void MainWindow::on_showMagenta_toggled(bool a_check)
{ ToogleView(MAGENTA_GROUP, a_check); }

void MainWindow::on_showAqua_toggled(bool a_check)
{ ToogleView(AQUA_GROUP, a_check); }


void MainWindow::on_showResult_toggled(bool a_check)
{ ToogleView(RESULT_GROUP, a_check); }

void MainWindow::on_showBlueComp_toggled(bool aCheck)
{
	if(aCheck)
	{
		m_color_complement = 0;
		showRedComp->setChecked(false);
		showBlackComp->setChecked(false);
		showBrownComp->setChecked(false);
		showYellowComp->setChecked(false);
		showMagentaComp->setChecked(false);
		showAquaComp->setChecked(false);
	} 
}

void MainWindow::on_showRedComp_toggled(bool aCheck)
{
	if(aCheck)
	{
		m_color_complement = 1;
		showBlueComp->setChecked(false);
		showBlackComp->setChecked(false);
		showBrownComp->setChecked(false);
		showYellowComp->setChecked(false);
		showMagentaComp->setChecked(false);
		showAquaComp->setChecked(false);
	} 
}

void MainWindow::on_showBlackComp_toggled(bool aCheck)
{
	if(aCheck)
	{
		m_color_complement = 2;
		showRedComp->setChecked(false);
		showBlueComp->setChecked(false);
		showBrownComp->setChecked(false);
		showYellowComp->setChecked(false);
		showMagentaComp->setChecked(false);
		showAquaComp->setChecked(false);
	} 
}

void MainWindow::on_showBrownComp_toggled(bool aCheck)
{
	if(aCheck)
	{
		m_color_complement = 3;
		showRedComp->setChecked(false);
		showBlackComp->setChecked(false);
		showBlueComp->setChecked(false);
		showYellowComp->setChecked(false);
		showMagentaComp->setChecked(false);
		showAquaComp->setChecked(false);
	} 
}

void MainWindow::on_showYellowComp_toggled(bool aCheck)
{
	if(aCheck)
	{
		m_color_complement = 4;
		showRedComp->setChecked(false);
		showBlackComp->setChecked(false);
		showBrownComp->setChecked(false);
		showBlueComp->setChecked(false);
		showMagentaComp->setChecked(false);
		showAquaComp->setChecked(false);
	} 
}

void MainWindow::on_showMagentaComp_toggled(bool aCheck)
{
	if(aCheck)
	{
		m_color_complement = 5;
		showRedComp->setChecked(false);
		showBlackComp->setChecked(false);
		showBrownComp->setChecked(false);
		showYellowComp->setChecked(false);
		showBlueComp->setChecked(false);
		showAquaComp->setChecked(false);
	} 
}

void MainWindow::on_showAquaComp_toggled(bool aCheck)
{
	if(aCheck)
	{
		m_color_complement = 6;
		showRedComp->setChecked(false);
		showBlackComp->setChecked(false);
		showBrownComp->setChecked(false);
		showYellowComp->setChecked(false);
		showMagentaComp->setChecked(false);
		showBlueComp->setChecked(false);
	} 
}



void MainWindow::on_showBlueInt_toggled(bool aCheck)
{
	if(aCheck)
	{
		m_blue_int = true;
	}
	else
	{
		m_blue_int = false;
	}
}
void MainWindow::on_showRedInt_toggled(bool aCheck)
{
	if(aCheck)
	{
		m_red_int = true;
	}

	else
	{
		m_red_int = false;
	}

}
void MainWindow::on_showBlackInt_toggled(bool aCheck)
{
	if(aCheck)
	{
		m_black_int = true;
	}
	else
	{
		m_black_int = false;
	}

}
void MainWindow::on_showBrownInt_toggled(bool aCheck)
{
	if(aCheck)
	{
		m_brown_int = true;
	}
	else
	{
		m_brown_int = false;
	}

}
void MainWindow::on_showYellowInt_toggled(bool aCheck)
{

	if(aCheck)
	{
		m_yellow_int = true;
	}
	else
	{
		m_yellow_int = false;
	}
}
void MainWindow::on_showMagentaInt_toggled(bool aCheck)
{
	if(aCheck)
	{
		m_magenta_int = true;
	}
	else
	{
		m_magenta_int = false;
	}

}
void MainWindow::on_showAquaInt_toggled(bool aCheck)
{
	if(aCheck)
	{
		m_aqua_int = true;
	}
	else
	{
		m_aqua_int = false;
	}

}



void MainWindow::on_showBlueUnion_toggled(bool aCheck)
{
	if(aCheck)
	{
		m_blue_union = true;
	}
	else
	{
		m_blue_union = false;
	}
}
void MainWindow::on_showRedUnion_toggled(bool aCheck)
{
	if(aCheck)
	{
		m_red_union = true;
	}
	else
	{
		m_red_union = false;
	}

}
void MainWindow::on_showBlackUnion_toggled(bool aCheck)
{

	if(aCheck)
	{
		m_black_union = true;
	}
	else
	{
		m_black_union = false;
	}
}
void MainWindow::on_showBrownUnion_toggled(bool aCheck)
{

	if(aCheck)
	{
		m_brown_union = true;
	}
	else
	{
		m_brown_union = false;
	}
}
void MainWindow::on_showYellowUnion_toggled(bool aCheck)
{

	if(aCheck)
	{
		m_yellow_union = true;
	}
	else
	{
		m_yellow_union = false;
	}
}
void MainWindow::on_showMagentaUnion_toggled(bool aCheck)
{

	if(aCheck)
	{
		m_magenta_union = true;
	}
	else
	{
		m_magenta_union = false;
	}
}
void MainWindow::on_showAquaUnion_toggled(bool aCheck)
{

	if(aCheck)
	{
		m_aqua_union = true;
	}
	else
	{
		m_aqua_union = false;
	}
}



void MainWindow::on_showBlueDiff_toggled(bool aCheck)
{
	if(aCheck)
	{
		size_t checkCount = 0;

		if(showAquaDiff->isChecked()) checkCount++;
		if(showRedDiff->isChecked()) checkCount++;
		if(showBlackDiff->isChecked()) checkCount++;
		if(showBrownDiff->isChecked()) checkCount++;
		if(showYellowDiff->isChecked()) checkCount++;
		if(showMagentaDiff->isChecked()) checkCount++;
		
		bool lProceed;
		if(checkCount >1)
		{
			lProceed = ask_user_yesno("Difference Operation Error", "More than 2 polygons are not supported for this operation. If you want to proceed further, please deselect one coloured polygon set.\nContinue?\n");
			if(lProceed)
			{
				if (m_color_diff_B == 2) showBlackDiff->setChecked(false);
				if (m_color_diff_B == 3) showBrownDiff->setChecked(false);
				if (m_color_diff_B == 4) showYellowDiff->setChecked(false);
				if (m_color_diff_B == 5) showMagentaDiff->setChecked(false);
				if (m_color_diff_B == 6) showAquaDiff->setChecked(false);
				m_color_diff_B = m_color_diff_A;
				m_color_diff_A = 0;
			}
			else
			{
				showBlueDiff->setChecked(false);
			}
		}
		else
		{
			m_color_diff_B = m_color_diff_A;
			m_color_diff_A = 0;
		}
	}
	else
	{
		showBlueDiff->setChecked(false);
		if(m_color_diff_B == 1) m_color_diff_A = 1;
		if(m_color_diff_B ==2) m_color_diff_A =2;
		if(m_color_diff_B ==3) m_color_diff_A =3;
		if(m_color_diff_B ==4) m_color_diff_A =4;
		if(m_color_diff_B ==5) m_color_diff_A =5;
		if(m_color_diff_B ==6) m_color_diff_A =6;
	}
}
void MainWindow::on_showRedDiff_toggled(bool aCheck)
{
	if(aCheck)
	{
		size_t checkCount = 0;

		if(showBlueDiff->isChecked()) checkCount++;
		if(showAquaDiff->isChecked()) checkCount++;
		if(showBlackDiff->isChecked()) checkCount++;
		if(showBrownDiff->isChecked()) checkCount++;
		if(showYellowDiff->isChecked()) checkCount++;
		if(showMagentaDiff->isChecked()) checkCount++;
		

		bool lProceed;
		if(checkCount >1)
		{
			lProceed = ask_user_yesno("Difference Operation Error", "More than 2 polygons are not supported for this operation. If you want to proceed further, please deselect one coloured polygon set.\nContinue?\n");
			if(lProceed)
			{
				if (m_color_diff_A == 0)
				{
					if (m_color_diff_B == 2) showBlackDiff->setChecked(false);
					if (m_color_diff_B == 3) showBrownDiff->setChecked(false);
					if (m_color_diff_B == 4) showYellowDiff->setChecked(false);
					if (m_color_diff_B == 5) showMagentaDiff->setChecked(false);
					if (m_color_diff_B == 6) showAquaDiff->setChecked(false);
					m_color_diff_B = 1;
				}
				else
				{
					if (m_color_diff_B == 3) showBrownDiff->setChecked(false);
					if (m_color_diff_B == 4) showYellowDiff->setChecked(false);
					if (m_color_diff_B == 5) showMagentaDiff->setChecked(false);
					if (m_color_diff_B == 6) showAquaDiff->setChecked(false);
					m_color_diff_B = m_color_diff_A;
					m_color_diff_A = 1;
				}
			}
			else
			{
				showRedDiff ->setChecked(false);
			}
		}
		else
		{
			if(m_color_diff_A <= 1)
			{
				m_color_diff_B = 1;
			}
			else
			{
				m_color_diff_B = m_color_diff_A;
				m_color_diff_A = 1;
			}
		}
	}
	else
	{
		showRedDiff->setChecked(false);
		if(m_color_diff_B == 1) m_color_diff_B = 0;
		if(m_color_diff_B == 2) m_color_diff_A = 2;
		if(m_color_diff_B == 3) m_color_diff_A = 3;
		if(m_color_diff_B == 4) m_color_diff_A = 4;
		if(m_color_diff_B == 5) m_color_diff_A = 5;
		if(m_color_diff_B == 6) m_color_diff_A = 6;
	}
}
void MainWindow::on_showBlackDiff_toggled(bool aCheck)
{
	if(aCheck)
	{
		size_t checkCount = 0;

		if(showBlueDiff->isChecked()) checkCount++;
		if(showRedDiff->isChecked()) checkCount++;
		if(showAquaDiff->isChecked()) checkCount++;
		if(showBrownDiff->isChecked()) checkCount++;
		if(showYellowDiff->isChecked()) checkCount++;
		if(showMagentaDiff->isChecked()) checkCount++;
		

		bool lProceed;
		if(checkCount >1)
		{
			lProceed = ask_user_yesno("Difference Operation Error", "More than 2 polygons are not supported for this operation. If you want to proceed further, please deselect one coloured polygon set.\nContinue?\n");
			if(lProceed)
			{
				if(m_color_diff_A == 0)
				{
					if(m_color_diff_B == 1)
					{
						lProceed = ask_user_yesno("Difference Operation Error"," You possibly selected a lower priority color polygon when 2 higher ones are selected.\nWould you like to replace the lower one with Black?\n");						
						if(lProceed)
						{
							showRedDiff->setChecked(false);
							m_color_diff_B = 2;
						}
						else
						{
							showBlackDiff ->setChecked(false);
						}							
					}
					if(m_color_diff_B >2)
					{
						if (m_color_diff_B == 3) showBrownDiff->setChecked(false);
						if (m_color_diff_B == 4) showYellowDiff->setChecked(false);
						if (m_color_diff_B == 5) showMagentaDiff->setChecked(false);
						if (m_color_diff_B == 6) showAquaDiff->setChecked(false);
						m_color_diff_B = 2;
					}
				}

				if (m_color_diff_A == 1)
				{
					if (m_color_diff_B == 3) showBrownDiff->setChecked(false);
					if (m_color_diff_B == 4) showYellowDiff->setChecked(false);
					if (m_color_diff_B == 5) showMagentaDiff->setChecked(false);
					if (m_color_diff_B == 6) showAquaDiff->setChecked(false);
					m_color_diff_B = 2;
				}
				if(m_color_diff_A > 2)
				{
					if (m_color_diff_B == 4) showYellowDiff->setChecked(false);
					if (m_color_diff_B == 5) showMagentaDiff->setChecked(false);
					if (m_color_diff_B == 6) showAquaDiff->setChecked(false);
					m_color_diff_B = m_color_diff_A;
					m_color_diff_A = 2;
				}
			}
			else
			{
				showBlackDiff->setChecked(false);
			}
		}
		else
		{
			if(m_color_diff_A <= 2)
			{
				m_color_diff_B = 2;
			}
			else
			{
				m_color_diff_B = m_color_diff_A;
				m_color_diff_A = 2;
			}
		}
	}
	else
	{
		showBlackDiff->setChecked(false);
		if(m_color_diff_A == 0) m_color_diff_B = 0;
		if(m_color_diff_A == 1) m_color_diff_B = 1;
		if(m_color_diff_B == 3) m_color_diff_A =3;
		if(m_color_diff_B == 4) m_color_diff_A =4;
		if(m_color_diff_B == 5) m_color_diff_A =5;
		if(m_color_diff_B == 6) m_color_diff_A =6;
	}	
}
void MainWindow::on_showBrownDiff_toggled(bool aCheck)
{
	if(aCheck)
	{
		size_t checkCount = 0;

		if(showBlueDiff->isChecked()) checkCount++;
		if(showRedDiff->isChecked()) checkCount++;
		if(showBlackDiff->isChecked()) checkCount++;
		if(showAquaDiff->isChecked()) checkCount++;
		if(showYellowDiff->isChecked()) checkCount++;
		if(showMagentaDiff->isChecked()) checkCount++;
		

		bool lProceed;
		if(checkCount >1)
		{
			lProceed = ask_user_yesno("Difference Operation Error", "More than 2 polygons are not supported for this operation. If you want to proceed further, please deselect one coloured polygon set.\nContinue?\n");
			if(lProceed)
			{
				if(m_color_diff_A == 0)
				{
					if(m_color_diff_B <= 2)
					{
						lProceed = ask_user_yesno("Difference Operation Error"," You possibly selected a lower priority color polygon when 2 higher ones are selected.\nWould you like to replace the lower one with Brown?\n");						
						if(lProceed)
						{
							if(m_color_diff_B == 1) showRedDiff->setChecked(false);
							if(m_color_diff_B == 2) showBlackDiff->setChecked(false);
							m_color_diff_B = 3;
						}
						else
						{
							showBrownDiff->setChecked(false);
						}
					}
					if(m_color_diff_B > 3)
					{
						if (m_color_diff_B == 4) showYellowDiff->setChecked(false);
						if (m_color_diff_B == 5) showMagentaDiff->setChecked(false);
						if (m_color_diff_B == 6) showAquaDiff->setChecked(false);
						m_color_diff_B = 3; 
					}
				}
				if (m_color_diff_A == 1)
				{
					if(m_color_diff_B == 2)					
					{
						lProceed = ask_user_yesno("Difference Operation Error"," You possibly selected a lower priority color polygon when 2 higher ones are selected.\nWould you like to replace the lower one with Brown ?\n");						
						if(lProceed)
						{
							if(m_color_diff_B == 2) showBlackDiff->setChecked(false);
							m_color_diff_B = 3;
						}
						else
						{
							showBrownDiff->setChecked(false);
						}
					}
					if(m_color_diff_B > 3)
					{
						if (m_color_diff_B == 4) showYellowDiff->setChecked(false);
						if (m_color_diff_B == 5) showMagentaDiff->setChecked(false);
						if (m_color_diff_B == 6) showAquaDiff->setChecked(false);
						m_color_diff_B = 3; 
					}
				}
				if(m_color_diff_A == 2)
				{
					if (m_color_diff_B == 4) showYellowDiff->setChecked(false);
					if (m_color_diff_B == 5) showMagentaDiff->setChecked(false);
					if (m_color_diff_B == 6) showAquaDiff->setChecked(false);
					m_color_diff_B = 3;
				}
				if(m_color_diff_A > 3)
				{
					if (m_color_diff_B == 5) showMagentaDiff->setChecked(false);
					if (m_color_diff_B == 6) showAquaDiff->setChecked(false);
					m_color_diff_B = m_color_diff_A;
					m_color_diff_A = 3;
				}
			}
			else
			{
				showBrownDiff->setChecked(false);
			}
		}
		else
		{
			if(m_color_diff_A <= 3)
			{
				m_color_diff_B = 3;
			}
			else
			{
				m_color_diff_B = m_color_diff_A;
				m_color_diff_A = 3;
			}
		}
	}
	else
	{
		showBrownDiff->setChecked(false);
		if(m_color_diff_A == 0) m_color_diff_B = 0;
		if(m_color_diff_A == 1) m_color_diff_B = 1;
		if(m_color_diff_A == 2) m_color_diff_B =2;
		if(m_color_diff_B == 4) m_color_diff_A =4;
		if(m_color_diff_B == 5) m_color_diff_A =5;
		if(m_color_diff_B == 6) m_color_diff_A =6;
	}	
}
void MainWindow::on_showYellowDiff_toggled(bool aCheck)
{
	if(aCheck)
	{
		size_t checkCount = 0;

		if(showBlueDiff->isChecked()) checkCount++;
		if(showRedDiff->isChecked()) checkCount++;
		if(showBlackDiff->isChecked()) checkCount++;
		if(showBrownDiff->isChecked()) checkCount++;
		if(showAquaDiff->isChecked()) checkCount++;
		if(showMagentaDiff->isChecked()) checkCount++;
		

		bool lProceed;
		if(checkCount >1)
		{
			lProceed = ask_user_yesno("Difference Operation Error", "More than 2 polygons are not supported for this operation. If you want to proceed further, please deselect one coloured polygon set.\nContinue?\n");
			if(lProceed)
			{
				if(m_color_diff_A == 0)
				{
					if(m_color_diff_B < 4)
					{
						lProceed = ask_user_yesno("Difference Operation Error"," You possibly selected a lower priority color polygon when 2 higher ones are selected.\nWould you like to replace the lower one with Yellow?\n");						
						if(lProceed)
						{
							if(m_color_diff_B == 1) showRedDiff->setChecked(false);
							if(m_color_diff_B == 2) showBlackDiff->setChecked(false);
							if(m_color_diff_B == 3) showBrownDiff->setChecked(false);
							m_color_diff_B = 4;
						}
						else
						{
							showYellowDiff->setChecked(false);
						}
					}
					if(m_color_diff_B >4)
					{
						if (m_color_diff_B == 5) showMagentaDiff->setChecked(false);
						if (m_color_diff_B == 6) showAquaDiff->setChecked(false);
						m_color_diff_B = 4;
					}
				}
				if(m_color_diff_A == 1)
				{
					if(m_color_diff_B < 4)
					{
						lProceed = ask_user_yesno("Difference Operation Error"," You possibly selected a lower priority color polygon when 2 higher ones are selected.\nWould you like to replace the lower one with Yellow?\n");						
						if(lProceed)
						{
							if(m_color_diff_B == 2) showBlackDiff->setChecked(false);
							if(m_color_diff_B == 3) showBrownDiff->setChecked(false);
							m_color_diff_B = 4;
						}
						else
						{
							showYellowDiff->setChecked(false);
						}
					}
					if(m_color_diff_B >4)
					{
						if (m_color_diff_B == 5) showMagentaDiff->setChecked(false);
						if (m_color_diff_B == 6) showAquaDiff->setChecked(false);
						m_color_diff_B = 4;
					}
				}
				if(m_color_diff_A == 2)
				{
					if(m_color_diff_B ==3 )
					{
						lProceed = ask_user_yesno("Difference Operation Error"," You possibly selected a lower priority color polygon when 2 higher ones are selected.\nWould you like to replace the lower one with Yellow?\n");						
						if(lProceed)
						{
							if(m_color_diff_B == 3) showBrownDiff->setChecked(false);
							m_color_diff_B = 4;
						}
						else
						{
							showYellowDiff->setChecked(false);
						}
					}
					if(m_color_diff_B >4)
					{
						if (m_color_diff_B == 5) showMagentaDiff->setChecked(false);
						if (m_color_diff_B == 6) showAquaDiff->setChecked(false);
						m_color_diff_B = 4;
					}
				}
				if(m_color_diff_A == 3)
				{
					if (m_color_diff_B == 5) showMagentaDiff->setChecked(false);
					if (m_color_diff_B == 6) showAquaDiff->setChecked(false);
					m_color_diff_B = 4;
				}
				if(m_color_diff_A > 4)
				{
					showAquaDiff->setChecked(false);
					m_color_diff_B = m_color_diff_A;
					m_color_diff_A = 4;
				}
			}
			else
			{
				showYellowDiff->setChecked(false);
			}
		}

		else
		{
			if(m_color_diff_A <= 4)
			{
				m_color_diff_B = 4;
			}
			else
			{
				m_color_diff_B = m_color_diff_A;
				m_color_diff_A = 4;
			}
		}
		
	}
	else
	{
		showYellowDiff->setChecked(false);
		if(m_color_diff_A == 0) m_color_diff_B = 0;
		if(m_color_diff_A == 1) m_color_diff_B = 1;
		if(m_color_diff_A == 2) m_color_diff_B =2;
		if(m_color_diff_A == 3) m_color_diff_B =3;
		if(m_color_diff_B == 5) m_color_diff_A =5;
		if(m_color_diff_B == 6) m_color_diff_A =6;
	}
}
void MainWindow::on_showMagentaDiff_toggled(bool aCheck)
{
	if(aCheck)
	{
		size_t checkCount = 0;

		if(showBlueDiff->isChecked()) checkCount++;
		if(showRedDiff->isChecked()) checkCount++;
		if(showBlackDiff->isChecked()) checkCount++;
		if(showBrownDiff->isChecked()) checkCount++;
		if(showYellowDiff->isChecked()) checkCount++;
		if(showAquaDiff->isChecked()) checkCount++;
		

		bool lProceed;
		if(checkCount >1)
		{
			lProceed = ask_user_yesno("Difference Operation Error", "More than 2 polygons are not supported for this operation. If you want to proceed further, please deselect one coloured polygon set.\nContinue?\n");
			if(lProceed)
			{
				if(m_color_diff_A == 0)
				{
					if(m_color_diff_B < 5)
					{
						lProceed = ask_user_yesno("Difference Operation Error"," You possibly selected a lower priority color polygon when 2 higher ones are selected.\nWould you like to replace the lower one with Magenta?\n");						
						if(lProceed)
						{
							if(m_color_diff_B == 1) showRedDiff->setChecked(false);
							if(m_color_diff_B == 2) showBlackDiff->setChecked(false);
							if(m_color_diff_B == 3) showBrownDiff->setChecked(false);
							if(m_color_diff_B == 4) showYellowDiff->setChecked(false); 
							m_color_diff_B = 5;
						}
						else
						{
							showMagentaDiff->setChecked(false);
						}
					}
					if(m_color_diff_B >5 )
					{
						showAquaDiff->setChecked(false);
						m_color_diff_B = 5;
					}
				}
				if(m_color_diff_A == 1)
				{
					if(m_color_diff_B < 5)
					{
						lProceed = ask_user_yesno("Difference Operation Error"," You possibly selected a lower priority color polygon when 2 higher ones are selected.\nWould you like to replace the lower one with Magenta?\n");						
						if(lProceed)
						{
							if(m_color_diff_B == 2) showBlackDiff->setChecked(false);
							if(m_color_diff_B == 3) showBrownDiff->setChecked(false);
							if(m_color_diff_B == 4) showYellowDiff->setChecked(false); 
							m_color_diff_B = 5;
						}
						else
						{
							showMagentaDiff->setChecked(false);
						}
					}
					if(m_color_diff_B >5 )
					{
						showAquaDiff->setChecked(false);
						m_color_diff_B = 5;
					}
				}
				if(m_color_diff_A == 2)
				{
					if(m_color_diff_B<5)
					{
						lProceed = ask_user_yesno("Difference Operation Error"," You possibly selected a lower priority color polygon when 2 higher ones are selected.\nWould you like to replace the lower one with Magenta?\n");						
						if(lProceed)
						{
							if(m_color_diff_B == 3) showBrownDiff->setChecked(false);
							if(m_color_diff_B == 4) showYellowDiff->setChecked(false); 
							m_color_diff_B = 5;
						}
						else
						{
							showMagentaDiff->setChecked(false);
						}
					}
					if(m_color_diff_B >5 )
					{
						showAquaDiff->setChecked(false);
						m_color_diff_B = 5;
					}
				}
				if(m_color_diff_A == 3)
				{
					if(m_color_diff_B < 5)
					{
						lProceed = ask_user_yesno("Difference Operation Error"," You possibly selected a lower priority color polygon when 2 higher ones are selected.\nWould you like to replace the lower one with Magenta?\n");						
						if(lProceed)
						{
							showYellowDiff->setChecked(false); 
							m_color_diff_B = 5;
						}
						else
						{
							showMagentaDiff->setChecked(false);
						}
					}
					if(m_color_diff_B >5 )
					{
						showAquaDiff->setChecked(false);
						m_color_diff_B = 5;
					}
				}
				if(m_color_diff_A == 4)
				{
					showAquaDiff->setChecked(false); 
					m_color_diff_B = 5;
				}
			}
		}

		else
		{
			if(m_color_diff_A <= 5)
			{
				m_color_diff_B = 5;
			}
			else
			{
				m_color_diff_A = 5;
			}
		}
	}
	else
	{
		showMagentaDiff->setChecked(false);
		if(m_color_diff_A == 0) m_color_diff_B = 0;
		if(m_color_diff_A == 1) m_color_diff_B = 1;
		if(m_color_diff_A == 2) m_color_diff_B =2;
		if(m_color_diff_A == 3) m_color_diff_B =3;
		if(m_color_diff_A == 4) m_color_diff_B =4;
		if(m_color_diff_B == 6) m_color_diff_A =6;
	}
}
void MainWindow::on_showAquaDiff_toggled(bool aCheck)
{
	if(aCheck)
	{
		size_t checkCount = 0;

		if(showBlueDiff->isChecked()) checkCount++;
		if(showRedDiff->isChecked()) checkCount++;
		if(showBlackDiff->isChecked()) checkCount++;
		if(showBrownDiff->isChecked()) checkCount++;
		if(showYellowDiff->isChecked()) checkCount++;
		if(showMagentaDiff->isChecked()) checkCount++;
		

		bool lProceed;
		if(checkCount >1)
		{
			lProceed = ask_user_yesno("Difference Operation Error","You possibly selected a lower priority color polygon when 2 higher ones are selected.\nWould you like to replace the lower one with Aqua?\n");						
			if(lProceed)
			{
				if(m_color_diff_B == 1) showRedDiff->setChecked(false);
				if(m_color_diff_B == 2) showBlackDiff->setChecked(false);
				if(m_color_diff_B == 3) showBrownDiff->setChecked(false);
				if(m_color_diff_B == 4) showYellowDiff->setChecked(false); 
				if(m_color_diff_B == 5) showMagentaDiff->setChecked(false);
				m_color_diff_B = 6;
			}
			else
			{
				showAquaDiff->setChecked(false);
			}
		}
		else
		{
			m_color_diff_B = 6;
		}
	}
	else
	{
		showAquaDiff->setChecked(false);
		if(m_color_diff_A == 0) m_color_diff_B = 0;
		if(m_color_diff_A == 1) m_color_diff_B = 1;
		if(m_color_diff_A == 2) m_color_diff_B =2;
		if(m_color_diff_A == 3) m_color_diff_B =3;
		if(m_color_diff_A == 4) m_color_diff_B =4;
		if(m_color_diff_A == 5) m_color_diff_B =5;
	}

}




void MainWindow::on_showBlueSym_Diff_toggled(bool aCheck)
{
	if(aCheck)
	{
		m_blue_sym_diff = true;
	}
	else
	{
		m_blue_sym_diff = false;
	}
}
void MainWindow::on_showRedSym_Diff_toggled(bool aCheck)
{

	if(aCheck)
	{
		m_red_sym_diff = true;
	}
	else
	{
		m_red_sym_diff = false;
	}
}
void MainWindow::on_showBlackSym_Diff_toggled(bool aCheck)
{
	if(aCheck)
	{
		m_black_sym_diff = true;
	}
	else
	{
		m_black_sym_diff = false;
	}
}
void MainWindow::on_showBrownSym_Diff_toggled(bool aCheck)
{
	if(aCheck)
	{
		m_brown_sym_diff = true;
	}
	else
	{
		m_brown_sym_diff = false;
	}
}
void MainWindow::on_showYellowSym_Diff_toggled(bool aCheck)
{

	if(aCheck)
	{
		m_yellow_sym_diff = true;
	}
	else
	{
		m_yellow_sym_diff = false;
	}
}
void MainWindow::on_showMagentaSym_Diff_toggled(bool aCheck)
{
	if(aCheck)
	{
		m_magenta_sym_diff = true;
	}
	else
	{
		m_magenta_sym_diff = false;
	}

}
void MainWindow::on_showAquaSym_Diff_toggled(bool aCheck)
{
	if(aCheck)
	{
		m_aqua_sym_diff = true;
	}
	else
	{
		m_aqua_sym_diff = false;
	}

}




void MainWindow::on_showBlueMink_Sum_toggled(bool aCheck)
{

}
void MainWindow::on_showRedMink_Sum_toggled(bool aCheck)
{

}
void MainWindow::on_showBlackMink_Sum_toggled(bool aCheck)
{

}
void MainWindow::on_showBrownMink_Sum_toggled(bool aCheck)
{

}
void MainWindow::on_showYellowMink_Sum_toggled(bool aCheck)
{

}
void MainWindow::on_showMagentaMink_Sum_toggled(bool aCheck)
{

}
void MainWindow::on_showAquaMink_Sum_toggled(bool aCheck)
{

}



//////////###################### Add Color Plus Button ###########################////////////
void MainWindow::on_actionAddColor_triggered()
{
	if(!m_visible_black)
	{
		m_visible_black = true;
		showBlack ->setVisible(true);
		drawBlack -> setVisible(true);
		showBlackComp -> setVisible(true);
		showBlackDiff -> setVisible(true);
		showBlackUnion -> setVisible(true);
		showBlackInt -> setVisible(true);
		showBlackSym_Diff -> setVisible(true);
		showBlackMink_Sum -> setVisible(true);


		line11 -> setVisible(true);

		line1 -> setGeometry(QRect(111,0,7,95));
  		line2 -> setGeometry(QRect(200,0,7,95));
  		line3 -> setGeometry(QRect(290,0,7,95));
  		line4 -> setGeometry(QRect(345,0,7,95));
  		line5 -> setGeometry(QRect(420,0,7,95));
  		line6 -> setGeometry(QRect(530,0,7,95));

  		actionMinusColor -> setText("Remove Black");
	}

	else if(!m_visible_brown)
	{
		m_visible_brown = true;
		showBrown ->setVisible(true);
		drawBrown -> setVisible(true);
		showBrownComp -> setVisible(true);
		showBrownDiff -> setVisible(true);
		showBrownUnion -> setVisible(true);
		showBrownInt -> setVisible(true);
		showBrownSym_Diff -> setVisible(true);
		showBrownMink_Sum -> setVisible(true);

		line10 -> setVisible(true);

		line1 -> setGeometry(QRect(111,0,7,125));
  		line2 -> setGeometry(QRect(200,0,7,125));
  		line3 -> setGeometry(QRect(290,0,7,125));
  		line4 -> setGeometry(QRect(345,0,7,125));
  		line5 -> setGeometry(QRect(420,0,7,125));
  		line6 -> setGeometry(QRect(530,0,7,125));

  		actionMinusColor -> setText("Remove Brown");
	}

	else if(!m_visible_yellow)
	{
		m_visible_yellow = true;
		showYellow -> setVisible(true);
		drawYellow -> setVisible(true);
		showYellowComp -> setVisible(true);
		showYellowDiff -> setVisible(true);
		showYellowUnion -> setVisible(true);
		showYellowInt -> setVisible(true);
		showYellowSym_Diff -> setVisible(true);
		showYellowMink_Sum -> setVisible(true);

		line9 -> setVisible(true);

		line1 -> setGeometry(QRect(111,0,7,155));
  		line2 -> setGeometry(QRect(200,0,7,155));
  		line3 -> setGeometry(QRect(290,0,7,155));
  		line4 -> setGeometry(QRect(345,0,7,155));
  		line5 -> setGeometry(QRect(420,0,7,155));
  		line6 -> setGeometry(QRect(530,0,7,155));

  		actionMinusColor -> setText("Remove Yellow");
	}

	else if(!m_visible_magenta)
	{
		m_visible_magenta = true;
		showMagenta -> setVisible(true);
		drawMagenta -> setVisible(true);
		showMagentaComp -> setVisible(true);
		showMagentaDiff -> setVisible(true);
		showMagentaUnion -> setVisible(true);
		showMagentaInt -> setVisible(true);
		showMagentaSym_Diff -> setVisible(true);
		showMagentaMink_Sum -> setVisible(true);

		line8 -> setVisible(true);

		line1 -> setGeometry(QRect(111,0,7,185));
  		line2 -> setGeometry(QRect(200,0,7,185));
  		line3 -> setGeometry(QRect(290,0,7,185));
  		line4 -> setGeometry(QRect(345,0,7,185));
  		line5 -> setGeometry(QRect(420,0,7,185));
  		line6 -> setGeometry(QRect(530,0,7,185));

  		actionMinusColor -> setText("Remove Magenta");
	}

	else if(!m_visible_aqua)
	{
		m_visible_aqua = true;
		showAqua ->setVisible(true);
		drawAqua -> setVisible(true);
		showAquaComp -> setVisible(true);
		showAquaDiff -> setVisible(true);
		showAquaUnion -> setVisible(true);
		showAquaInt -> setVisible(true);
		showAquaSym_Diff -> setVisible(true);
		showAquaMink_Sum -> setVisible(true);	

		line7 -> setVisible(true);

		line1 -> setGeometry(QRect(111,0,7,220));
  		line2 -> setGeometry(QRect(200,0,7,220));
  		line3 -> setGeometry(QRect(290,0,7,220));
  		line4 -> setGeometry(QRect(345,0,7,220));
  		line5 -> setGeometry(QRect(420,0,7,220));
  		line6 -> setGeometry(QRect(530,0,7,220));

  		actionMinusColor -> setText("Remove Aqua");
	}

	else
	{
		ask_user_yesno("Adding Color Error","Maximum Limit of Color has reached.\n");						
	}

}


//////#########Remove Color Minus

void MainWindow::on_actionMinusColor_triggered()
{
	if (!(m_visible_brown || m_visible_yellow || m_visible_magenta || m_visible_aqua))
	{	
		m_curve_sets[2].clear();
		m_curve_sets[7].clear();
		m_color_active = 1;

		black_circular_sources().clear();
		black_linear_sources().clear();
		black_bezier_sources().clear();

		m_visible_black = false;
		showBlack ->setVisible(false);
		drawBlack -> setVisible(false);
		drawBlue -> setChecked(true);
		showBlackComp -> setVisible(false);
		showBlackDiff -> setVisible(false);
		showBlackUnion -> setVisible(false);
		showBlackInt -> setVisible(false);
		showBlackSym_Diff -> setVisible(false);
		showBlackMink_Sum -> setVisible(false);

		


		line11 -> setVisible(false);

		line1 -> setGeometry(QRect(111,0,7,65));
	  	line2 -> setGeometry(QRect(200,0,7,65));
	  	line3 -> setGeometry(QRect(290,0,7,65));
	  	line4 -> setGeometry(QRect(345,0,7,65));
	  	line5 -> setGeometry(QRect(420,0,7,65));
	  	line6 -> setGeometry(QRect(530,0,7,65));	

	  	actionMinusColor -> setText("Color Removal Not Allowed");

	  	modelChanged();

	}

	else if (!(m_visible_yellow || m_visible_magenta || m_visible_aqua))
	{	
		m_curve_sets[3].clear();
		m_curve_sets[7].clear();
		m_color_active = 2;

		brown_circular_sources().clear();
		brown_linear_sources().clear();
		brown_bezier_sources().clear();

		m_visible_brown = false;
		showBrown ->setVisible(false);
		drawBrown -> setVisible(false);
		showBrownComp -> setVisible(false);
		showBrownDiff -> setVisible(false);
		showBrownUnion -> setVisible(false);
		showBrownInt -> setVisible(false);
		showBrownSym_Diff -> setVisible(false);
		showBrownMink_Sum -> setVisible(false);
		

		drawBlue -> setChecked(true);


		line10 -> setVisible(false);

		line1 -> setGeometry(QRect(111,0,7,95));
	  	line2 -> setGeometry(QRect(200,0,7,95));
	  	line3 -> setGeometry(QRect(290,0,7,95));
	  	line4 -> setGeometry(QRect(345,0,7,95));
	  	line5 -> setGeometry(QRect(420,0,7,95));
	  	line6 -> setGeometry(QRect(530,0,7,95));

	  	actionMinusColor -> setText("Remove Black");

	  	modelChanged();

	 }


	else if (!(m_visible_magenta || m_visible_aqua))
	{	
		m_curve_sets[4].clear();
		m_curve_sets[7].clear();
		m_color_active = 3;

		yellow_circular_sources().clear();
		yellow_linear_sources().clear();
		yellow_bezier_sources().clear();

		m_visible_yellow = false;
		showYellow ->setVisible(false);
		drawYellow -> setVisible(false);
		showYellowComp -> setVisible(false);
		showYellowDiff -> setVisible(false);
		showYellowUnion -> setVisible(false);
		showYellowInt -> setVisible(false);
		showYellowSym_Diff -> setVisible(false);
		showYellowMink_Sum -> setVisible(false);

		


		drawBlue -> setChecked(true);


		line9 -> setVisible(false);

		line1 -> setGeometry(QRect(111,0,7,125));
	  	line2 -> setGeometry(QRect(200,0,7,125));
	  	line3 -> setGeometry(QRect(290,0,7,125));
	  	line4 -> setGeometry(QRect(345,0,7,125));
	  	line5 -> setGeometry(QRect(420,0,7,125));
	  	line6 -> setGeometry(QRect(530,0,7,125));

	  	actionMinusColor -> setText("Remove Brown");	

	  	modelChanged();
	}

	else if (! m_visible_aqua)
	{	
		m_curve_sets[5].clear();
		m_curve_sets[7].clear();
		m_color_active = 4;

		magenta_circular_sources().clear();
		magenta_linear_sources().clear();
		magenta_bezier_sources().clear();

		m_visible_magenta = false;
		showMagenta ->setVisible(false);
		drawMagenta -> setVisible(false);
		showMagentaComp -> setVisible(false);
		showMagentaDiff -> setVisible(false);
		showMagentaUnion -> setVisible(false);
		showMagentaInt -> setVisible(false);
		showMagentaSym_Diff -> setVisible(false);
		showMagentaMink_Sum -> setVisible(false);

		

		drawBlue -> setChecked(true);

		line8 -> setVisible(false);

		line1 -> setGeometry(QRect(111,0,7,155));
	  	line2 -> setGeometry(QRect(200,0,7,155));
	  	line3 -> setGeometry(QRect(290,0,7,155));
	  	line4 -> setGeometry(QRect(345,0,7,155));
	  	line5 -> setGeometry(QRect(420,0,7,155));
	  	line6 -> setGeometry(QRect(530,0,7,155));	

	  	actionMinusColor -> setText("Remove Yellow");

	  	modelChanged();

	}

	else
	{

		m_curve_sets[6].clear();
		m_curve_sets[7].clear();
		m_color_active = 5;

		aqua_circular_sources().clear();
		aqua_linear_sources().clear();
		aqua_bezier_sources().clear();


		m_visible_aqua = false;
		showAqua ->setVisible(false);
		drawAqua -> setVisible(false);
		showAquaComp -> setVisible(false);
		showAquaDiff -> setVisible(false);
		showAquaUnion -> setVisible(false);
		showAquaInt -> setVisible(false);
		showAquaSym_Diff -> setVisible(false);
		showAquaMink_Sum -> setVisible(false);

		

		drawBlue -> setChecked(true);

		line7 -> setVisible(false);

		line1 -> setGeometry(QRect(111,0,7,185));
	  	line2 -> setGeometry(QRect(200,0,7,185));
	  	line3 -> setGeometry(QRect(290,0,7,185));
	  	line4 -> setGeometry(QRect(345,0,7,185));
	  	line5 -> setGeometry(QRect(420,0,7,185));
	  	line6 -> setGeometry(QRect(530,0,7,185));	

	  	actionMinusColor -> setText("Remove Magenta");

	  	modelChanged();
	}
}


//////////////#################################
void MainWindow::on_actionNew_triggered() 
{
  for( Curve_set_iterator si = m_curve_sets.begin(); si != m_curve_sets.end() ; ++ si )
    si->clear();
    
  blue_circular_sources().clear();
  red_circular_sources().clear();
  black_circular_sources().clear();
  brown_circular_sources().clear();
  yellow_circular_sources().clear();
  magenta_circular_sources().clear();
  aqua_circular_sources().clear();
  univ_circular_sources().clear();

  blue_linear_sources().clear();
  red_linear_sources().clear();
  black_linear_sources().clear();
  brown_linear_sources().clear();
  yellow_linear_sources().clear();
  magenta_linear_sources().clear();
  aqua_linear_sources().clear();
  univ_linear_sources().clear();

  blue_bezier_sources().clear();
  red_bezier_sources().clear();
  black_bezier_sources().clear();
  brown_bezier_sources().clear();
  yellow_bezier_sources().clear();
  magenta_bezier_sources().clear();
  aqua_bezier_sources().clear();
  univ_bezier_sources().clear();
    
  SetViewBlue  (true);
  SetViewRed   (true);
  SetViewBlack (true);
  SetViewBrown (true);
  SetViewYellow (true);
  SetViewMagenta (true);
  SetViewAqua (true);
  SetViewResult(true);
  
  m_circular_active = false ;
  m_bezier_active = false;
  
  m_color_active = 0;

  actionComplement -> setChecked(false);
  actionUnion -> setChecked(false);
  actionIntersection -> setChecked(false);
  actionDifference -> setChecked(false);
  actionSymmetric_Difference -> setChecked(false);
  actionMinkowski_Sum -> setChecked(false);

  actionInsertLinear -> setChecked(true);

  m_color_complement = 0; //default
  m_blue_int = true; //default
  m_red_int = false ; //default
  m_black_int = false; //default
  m_brown_int = false; //default
  m_yellow_int = false; //default
  m_magenta_int = false; //default
  m_aqua_int = false; //default
  m_blue_union = true; //default
  m_red_union = false; //default
  m_black_union = false; //default
  m_brown_union = false; //default
  m_yellow_union = false; //default
  m_magenta_union = false; //default
  m_aqua_union = false; //default
  m_color_diff_A = 0; //default
  m_color_diff_B = 0; //default
  m_blue_sym_diff = true; //default
  m_red_sym_diff = false; //default
  m_black_sym_diff = false; //default
  m_brown_sym_diff = false; //default
  m_yellow_sym_diff = false; //default
  m_magenta_sym_diff = false; //default
  m_aqua_sym_diff = false; //default

  drawBlue -> setChecked(true);

  showBlue->setChecked(true);
  showRed->setChecked(true);
  showBlack->setChecked(true);
  showBrown->setChecked(true);
  showYellow->setChecked(true);
  showBrown->setChecked(true);
  showMagenta->setChecked(true);
  showAqua->setChecked(true);
  showResult->setChecked(true);

  showBlueComp->setChecked(true);
  showRedComp->setChecked(false);
  showBlackComp->setChecked(false);
  showBrownComp->setChecked(false);
  showYellowComp->setChecked(false);
  showMagentaComp->setChecked(false);
  showAquaComp->setChecked(false);

  showBlueInt->setChecked(true);
  showRedInt->setChecked(false);
  showBlackInt->setChecked(false);
  showBrownInt->setChecked(false);
  showYellowInt->setChecked(false);
  showMagentaInt->setChecked(false);
  showAquaInt->setChecked(false);

  showBlueUnion->setChecked(true);
  showRedUnion->setChecked(false);
  showBlackUnion->setChecked(false);
  showBrownUnion->setChecked(false);
  showYellowUnion->setChecked(false);
  showMagentaUnion->setChecked(false);
  showAquaUnion->setChecked(false);

  showBlueDiff->setChecked(true);
  showRedDiff->setChecked(false);
  showBlackDiff->setChecked(false);
  showBrownDiff->setChecked(false);
  showYellowDiff->setChecked(false);
  showMagentaDiff->setChecked(false);
  showAquaDiff->setChecked(false);

  showBlueSym_Diff->setChecked(true);
  showRedSym_Diff->setChecked(false);
  showBlackSym_Diff->setChecked(false);
  showBrownSym_Diff->setChecked(false);
  showYellowSym_Diff->setChecked(false);
  showMagentaSym_Diff->setChecked(false);
  showAquaSym_Diff->setChecked(false);

  m_visible_black = false;
  showBlack ->setVisible(false);
  //drawBlack -> setVisible(false);
  showBlackComp -> setVisible(false);
  showBlackDiff -> setVisible(false);
  showBlackUnion -> setVisible(false);
  showBlackInt -> setVisible(false);
  showBlackSym_Diff -> setVisible(false);
  showBlackMink_Sum -> setVisible(false);
  

  
  m_visible_brown = false;
  showBrown ->setVisible(false);
  //drawBrown -> setVisible(false);
  showBrownComp -> setVisible(false);
  showBrownDiff -> setVisible(false);
  showBrownUnion -> setVisible(false);
  showBrownInt -> setVisible(false);
  showBrownSym_Diff -> setVisible(false);
  showBrownMink_Sum -> setVisible(false);


  m_visible_yellow = false;
  showYellow ->setVisible(false);
  //drawYellow -> setVisible(false);
  showYellowComp -> setVisible(false);
  showYellowDiff -> setVisible(false);
  showYellowUnion -> setVisible(false);
  showYellowInt -> setVisible(false);
  showYellowSym_Diff -> setVisible(false);
  showYellowMink_Sum -> setVisible(false);


  m_visible_magenta = false;
  showMagenta ->setVisible(false);
  //drawMagenta -> setVisible(false);
  showMagentaComp -> setVisible(false);
  showMagentaDiff -> setVisible(false);
  showMagentaUnion -> setVisible(false);
  showMagentaInt -> setVisible(false);
  showMagentaSym_Diff -> setVisible(false);
  showMagentaMink_Sum -> setVisible(false);  
  

  m_visible_aqua = false;
  showAqua ->setVisible(false);
  //drawAqua -> setVisible(false);
  showAquaComp -> setVisible(false);
  showAquaDiff -> setVisible(false);
  showAquaUnion -> setVisible(false);
  showAquaInt -> setVisible(false);
  showAquaSym_Diff -> setVisible(false);
  showAquaMink_Sum -> setVisible(false);
  
  actionMinusColor -> setText("Color Removal Not Allowed");

  line7 -> setVisible(false);
  line8 -> setVisible(false);
  line9 -> setVisible(false);
  line10 -> setVisible(false);
  line11 -> setVisible(false);
  

  line1 -> setGeometry(QRect(111,0,7,65));
  line2 -> setGeometry(QRect(200,0,7,65));
  line3 -> setGeometry(QRect(290,0,7,65));
  line4 -> setGeometry(QRect(345,0,7,65));
  line5 -> setGeometry(QRect(420,0,7,65));
  line6 -> setGeometry(QRect(530,0,7,65));

  
  modelChanged();
}

void MainWindow::on_actionDelete_triggered()
{
  bool lDone = false;
  bool lProceed=result_set().is_empty() ? true:ask_user_yesno("Store result", "All polygons of the selected type will be deleted\n continue anyway?\n") ;
  if (lProceed) {
    switch(m_color_active) {
     /*case 0: blue_set().assign(result_set()); break;
     case 1: red_set().assign(result_set()); break;
     case 2: black_set().assign(result_set()); break;
     case 3: brown_set().assign(result_set()); break;
     case 4: yellow_set().assign(result_set()); break;
     case 5: magenta_set().assign(result_set()); break;
     case 6: aqua_set().assign(result_set()); break;*/
     case 0: blue_set().clear();blue_circular_sources().clear();blue_bezier_sources().clear();blue_linear_sources().clear();break;
     case 1: red_set().clear();red_circular_sources().clear();red_bezier_sources().clear();red_linear_sources().clear();break;
     case 2: black_set().clear();black_circular_sources().clear();black_bezier_sources().clear();black_linear_sources().clear();break;
     case 3: brown_set().clear();brown_circular_sources().clear();brown_bezier_sources().clear();brown_linear_sources().clear();break;
     case 4: yellow_set().clear();yellow_circular_sources().clear();yellow_bezier_sources().clear();yellow_linear_sources().clear();break;
     case 5: magenta_set().clear();magenta_circular_sources().clear();magenta_bezier_sources().clear();magenta_linear_sources().clear();break;
     case 6: aqua_set().clear();aqua_circular_sources().clear();aqua_bezier_sources().clear();aqua_linear_sources().clear();break;
     case 8: universal_set().clear();univ_circular_sources().clear();univ_bezier_sources().clear();univ_linear_sources().clear();break;

     default: break;    //! \todo Handle default case.
    }
    //result_set().clear();
    lDone = true;
  }
  if (lDone) modelChanged();
}

void MainWindow::on_actionDeleteAll_triggered()
{
  bool lDone = false;
  bool lProceed = result_set().is_empty() ? true : ask_user_yesno("Store result","All polygons will be deleted\n continue anyway?\n");
  if (lProceed) {
    blue_set().clear();blue_circular_sources().clear();blue_bezier_sources().clear();blue_linear_sources().clear();
    red_set().clear();red_circular_sources().clear();red_bezier_sources().clear();red_linear_sources().clear();
    black_set().clear();black_circular_sources().clear();black_bezier_sources().clear();black_linear_sources().clear();
    brown_set().clear();brown_circular_sources().clear();brown_bezier_sources().clear();brown_linear_sources().clear();
    yellow_set().clear();yellow_circular_sources().clear();yellow_bezier_sources().clear();yellow_linear_sources().clear();
    magenta_set().clear();magenta_circular_sources().clear();magenta_bezier_sources().clear();magenta_linear_sources().clear();
    aqua_set().clear();aqua_circular_sources().clear();aqua_bezier_sources().clear();aqua_linear_sources().clear();
    universal_set().clear();univ_circular_sources().clear();univ_bezier_sources().clear();univ_linear_sources().clear();
    //result_set().clear();
  }
  lDone = true;

  if (lDone) modelChanged();
}

void MainWindow::on_actionDeleteResult_triggered()
{
    bool lDone  = false;
    bool lProceed = result_set().is_empty() ? true : ask_user_yesno("Store result","Result will be deleted\n continue anyway?\n");
    if (lProceed)
    {
      result_set().clear();
    }
    lDone = true;
    if (lDone) modelChanged();
}

void MainWindow::on_drawBlue_toggled(bool /* a_check */) { m_color_active = 0; }
void MainWindow::on_drawRed_toggled(bool /* a_check */) { m_color_active = 1; }
void MainWindow::on_drawBlack_toggled(bool /* a_check */) { m_color_active = 2; }
void MainWindow::on_drawBrown_toggled(bool /* a_check */) { m_color_active = 3; }
void MainWindow::on_drawYellow_toggled(bool /* a_check */) { m_color_active = 4; }
void MainWindow::on_drawMagenta_toggled(bool /* a_check */) { m_color_active = 5; }
void MainWindow::on_drawAqua_toggled(bool /* a_check */) { m_color_active = 6; }
void MainWindow::on_drawUniversalPolygon_toggled(bool /* a_check */) { m_color_active = 8; } 

//extra utilities
void MainWindow::on_actionRecenter_triggered() { zoomToFit(); }

void MainWindow::dragEnterEvent(QDragEnterEvent* event)
{
  if (event->mimeData()->hasFormat("text/uri-list"))
    event->acceptProposedAction();
}

void MainWindow::dropEvent(QDropEvent* event)
{
  QString filename = event->mimeData()->urls().at(0).path();
  open(filename);
  event->acceptProposedAction();
}

void MainWindow::on_actionOpenLinear_triggered()
{
  open(QFileDialog::getOpenFileName(this,
                                    tr("Open Linear Polygon"), "../data",
                                    tr("Linear Curve files (*.lps)")));
}

void MainWindow::on_actionOpenDXF_triggered()
{
  open(QFileDialog::getOpenFileName(this,
                                    tr("Open Circular Polygon"), "../data",
                                    tr("Circular curve files (*.dxf)")));
}

void MainWindow::on_actionOpenBezier_triggered()
{
  open(QFileDialog::getOpenFileName(this, 
                                    tr("Open Bezier Polygon"), "../data", 
                                    tr("Bezier Curve files (*.bps)") ));
}


//for converting linear part of circular polygon to circular part
/*Circular_polygon linearPart_2_circ(Circular_Linear_polygon const& pgn)
{
  CGAL::Cartesian_converter<Kernel,Kernel> convert;
  Circular_polygon rCP;
  for (auto ei = pgn.edges_begin(); ei != pgn.edges_end(); ++ei) {
    if (ei->source() != ei->target())
      rCP.push_back(Circular_X_monotone_curve(convert(ei->source()),
                                              convert(ei->target())));
  }
  return rCP;
}
//for converting linear part of circular polygon with holes to circular part
Circular_polygon_with_holes
linearPart_2_circ(Circular_Linear_polygon_with_holes const& pwh)
{
  Circular_polygon_with_holes rCP(linearPart_2_circ(pwh.outer_boundary()));
  for (auto hi = pwh.holes_begin(); hi != pwh.holes_end(); ++ hi)
    rCP.add_hole(linearPart_2_circ(*hi) );
  return rCP;
}*/

/*Circular_polygon linear_2_circ( Linear_polygon const& pgn )
{
  CGAL::Cartesian_converter<Linear_kernel,Gps_circular_kernel> convert ;
  
  Circular_polygon rCP;
  
  for( Linear_polygon::Edge_const_iterator ei = pgn.edges_begin(); ei != pgn.edges_end(); ++ei )
  {
    if  ( ei->source() != ei->target() )
      rCP.push_back( Circular_X_monotone_curve( convert(ei->source()), convert(ei->target())) );
  }  
  return rCP;
}
Circular_polygon_with_holes linear_2_circ( Linear_polygon_with_holes const& pwh )
{
  Circular_polygon_with_holes rCP( linear_2_circ(pwh.outer_boundary()) ) ;
  
  for( Linear_polygon_with_holes::Hole_const_iterator hi = pwh.holes_begin(); hi != pwh.holes_end(); ++ hi )
    rCP.add_hole( linear_2_circ(*hi)  );
  return rCP;
}
bool read_linear( QString aFileName, Circular_polygon_set& rSet, Circular_region_source_container& rSources )
{
  bool rOK = false ;
  
  std::ifstream in_file (qPrintable(aFileName));
  if ( in_file )
  {
    unsigned int n_regions ;
    in_file >> n_regions;
    
    for ( unsigned int r = 0 ; r < n_regions ; ++ r )
    {
      unsigned int n_boundaries;
      in_file >> n_boundaries;
      
      Circular_polygon outer ;
      std::vector<Circular_polygon> holes ;
      
      for ( unsigned int r = 0 ; r < n_boundaries ; ++ r )
      {
        Linear_polygon p ;
        in_file >> p ;
        
        if ( r == 0 )
             outer = linear_2_circ(p); 
        else holes.push_back( linear_2_circ(p) );
      }
      
      Circular_polygon_with_holes pwh(outer,holes.begin(),holes.end());
      rSources.push_back(pwh);
      rSet.join(pwh) ;    
      rOK = true ;
    }
    
  }
  
  return rOK ;
}
bool read_dxf ( QString aFileName, Circular_polygon_set& rSet, Circular_region_source_container& rSources )
{
  bool rOK = false ;
  
  std::ifstream in_file (qPrintable(aFileName));
  if ( in_file )
  {
    CGAL::Dxf_bsop_reader<Gps_circular_kernel>   reader;
    std::vector<Circular_polygon>            circ_polygons;
    std::vector<Circular_polygon_with_holes> circ_polygons_with_holes;
    
    reader(in_file
          ,std::back_inserter(circ_polygons)
          ,std::back_inserter(circ_polygons_with_holes)
          ,false
          );
          
    for ( std::vector<Circular_polygon>::iterator pit = circ_polygons.begin() ; pit != circ_polygons.end() ; ++ pit )
      circ_polygons_with_holes.push_back( Circular_polygon_with_holes(*pit) ) ;
    rSet.join( circ_polygons_with_holes.begin(), circ_polygons_with_holes.end() ) ;
    std::copy(circ_polygons_with_holes.begin(), circ_polygons_with_holes.end(), std::back_inserter(rSources) );
    rOK = true ;
  }
  
  return rOK ;
}
*/

/*bool read_bezier ( QString aFileName, Bezier_polygon_set& rSet, Bezier_region_source_container& rSources  )
{
  
  bool rOK = false ;
  
  std::ifstream in_file (qPrintable(aFileName));
  
  if ( in_file )
  {
    try
    {
      
      std::string format ;
      std::getline(in_file,format);
      
      bool lDoubleFormat = ( format.length() >= 6 && format.substr(0,6) == "DOUBLE") ;
                              
      // Red the number of bezier polygon with holes
      unsigned int n_regions ;
      in_file >> n_regions;
      
      for ( unsigned int r = 0 ; r < n_regions ; ++ r )
      {
        Bezier_polygon_vector bezier_polygons ;
        Bezier_region_source  br_source ;
        // Read the number of bezier curves.
        unsigned int n_boundaries;
        in_file >> n_boundaries;
      
        for ( unsigned int b = 0 ; b < n_boundaries ; ++ b )
        {
          Bezier_boundary_source bb_source ;        
          
          // Read the number of bezier curves.
          unsigned int n_curves;
          in_file >> n_curves;
          
          // Read the curves one by one, and construct the general polygon these
          // curve form (the outer boundary and the holes inside it).
          
          std::list<Bezier_X_monotone_curve> xcvs;
        
          for ( unsigned int k = 0; k < n_curves; ++ k ) 
          {
            // Read the current curve and subdivide it into x-monotone subcurves.
            
            std::list<CGAL::Object>                 x_objs;
            std::list<CGAL::Object>::const_iterator xoit;
            Bezier_X_monotone_curve                 xcv;
            Bezier_traits                           traits;
            Bezier_traits::Make_x_monotone_2        make_x_monotone = traits.make_x_monotone_2_object();
        
            Bezier_curve b = read_bezier_curve(in_file, lDoubleFormat);
            if ( b.number_of_control_points() >= 2 )
            {
              bb_source.push_back(b);
              //TRACE( "region " << r << " boundary " << b << " curve " << k );
                
              make_x_monotone (b, std::back_inserter (x_objs));
              
              for (xoit = x_objs.begin(); xoit != x_objs.end(); ++xoit) 
              {
                if (CGAL::assign (xcv, *xoit))
                {
                  //TRACE( " X montonote: " << xcv.source() << " -> " << xcv.target() << ( xcv.is_directed_right() ? " RIGHT":" LEFT") << ( xcv.is_vertical() ? " VERTICAL" : "")) ;
                  xcvs.push_back (xcv);
                }  
              }
            }
          }  
            
          Bezier_polygon  pgn (xcvs.begin(), xcvs.end());
          
          CGAL::Orientation  orient = pgn.orientation();
          //TRACE( "  Orientation: " << orient ) ;
            
          if (( b == 0 && orient == CGAL::CLOCKWISE) || ( b > 0 && orient == CGAL::COUNTERCLOCKWISE))
          {
            //TRACE( "Reversing orientation: " ) ;
            pgn.reverse_orientation();
          }
            
          br_source.push_back(bb_source);
          bezier_polygons.push_back (pgn);
        }
      
        if ( bezier_polygons.size() > 0 )
        {
          Bezier_polygon_with_holes pwh(bezier_polygons.front());
          
          if ( bezier_polygons.size() > 1 )
          {
            for ( Bezier_polygon_vector::const_iterator it = std::next(bezier_polygons.begin())
                ; it != bezier_polygons.end()
                ; ++ it 
                )
              pwh.add_hole(*it);    
          }
          
          if ( is_valid_polygon_with_holes(pwh, rSet.traits() ) )
          {
            rSet.join(pwh) ;      
            rSources.push_back(br_source);
          }
          else
          {
            show_warning( "Bezier polygon is not valid" );
          }
        }
        
        rOK = true ;
      }
      
    }
    catch(...)
    {
      show_error("Exception ocurred during reading of bezier polygon set.");
    } 
  }
  
  return rOK ;
}
Bezier_curve read_bezier_curve ( std::istream& is, bool aDoubleFormat )
{
  // Read the number of control points.
  unsigned int  n;
  is >> n;
  
  // Read the control points.
  std::vector<Bezier_rat_point> ctrl_pts;
  
  for ( unsigned int k = 0; k < n; k++)
  {
    Bezier_rat_point p ;
    if ( aDoubleFormat )
    {
      double x,y ;
      is >> x >> y ;
      Bezier_rational rx(static_cast<int> (1000 * x + 0.5), 1000);
      Bezier_rational ry(static_cast<int> (1000 * y + 0.5), 1000); 
      p = Bezier_rat_point(rx,ry);
    }
    else
    {
      is >> p ;
    }
    
    if ( k == 0 || ctrl_pts[k-1] != p ) 
    {
      ctrl_pts.push_back(p) ;
    }
  }
  std::vector<Bezier_rat_point> ctrl_pts2;
  typedef std::vector<Bezier_rat_point>::const_iterator cp_const_iterator ;
  cp_const_iterator beg  = ctrl_pts.begin();
  cp_const_iterator end  = ctrl_pts.end  ();
  cp_const_iterator last = end - 1 ;
  ctrl_pts2.push_back(*beg);
  if ( ctrl_pts.size() > 2 )
  {
    cp_const_iterator curr = beg ;
    cp_const_iterator next1 = curr  + 1 ;
    cp_const_iterator next2 = next1 + 1 ;
    do
    {
      CGAL::Orientation lOrient = orientation(*curr,*next1,*next2);
      if ( lOrient != CGAL::COLLINEAR )
        ctrl_pts2.push_back(*next1);
      ++ curr  ;
      ++ next1 ; 
      ++ next2 ;
    }
    while ( next2 != end ) ;
  }
  ctrl_pts2.push_back(*last);
  return Bezier_curve(ctrl_pts2.begin(),ctrl_pts2.end());
}*/


bool save_linear ( QString aFileName, Linear_polygon_set& rSet )
{
  bool rOK = false ;
  
  return rOK ;
}

bool save_circular ( QString aFileName, Circular_polygon_set& rSet )
{
  bool rOK = false ;
  
  return rOK ;
}

void save_bezier_polygon( std::ostream& out_file, Bezier_polygon const& aBP )
{
  typedef std::vector<Bezier_rat_point> Bezier_rat_point_vector ;
  
  int cc = aBP.size() ;
  int lc = cc - 1 ;
  
  out_file << "  " <<  cc << std::endl ;
  
  Bezier_rat_point lFirstP, lPrevP ;
 
  int i = 0 ;
  
  for ( Bezier_polygon::Curve_const_iterator cit = aBP.curves_begin() ; cit != aBP.curves_end() ; ++ cit, ++ i  )
  {
    Bezier_rat_point_vector lQ ;

    CGAL::Qt::Bezier_helper::clip(*cit,lQ);  
    
    out_file << "   " << lQ.size() << std::endl ;
    
    if ( i == 0 )
      lFirstP = lQ.front();
      
    if ( i == lc )  
      lQ.back() = lFirstP ;
      
    for ( Bezier_rat_point_vector::const_iterator pit = lQ.begin() ; pit != lQ.end() ; ++ pit )
    {
      Bezier_rat_point lP = pit == lQ.begin() && i > 0 ? lPrevP : *pit ;
      
      out_file << "    " << CGAL::to_double(lP.x()) << " " << CGAL::to_double(lP.y()) << std::endl ;
      
      lPrevP = lP ;
    }
  }
}

bool save_bezier_result ( QString aFileName, Bezier_polygon_set const& aSet )
{
  bool rOK = false ;
  
  std::ofstream out_file( qPrintable(aFileName) ) ;
  if ( out_file )
  {
    out_file << "DOUBLE" << std::endl ;
    
    std::vector<Bezier_polygon_with_holes> bpwh_container;
  
    aSet.polygons_with_holes( std::back_inserter(bpwh_container) ) ;
  
    out_file << bpwh_container.size() << std::endl ;

    for( std::vector<Bezier_polygon_with_holes>::const_iterator rit = bpwh_container.begin(); rit != bpwh_container.end() ; ++ rit )
    {
      Bezier_polygon_with_holes bpwh = *rit ;
      
      out_file << " " << ( 1 + bpwh.number_of_holes() ) << std::endl ;
  
      save_bezier_polygon( out_file, bpwh.outer_boundary() ) ;
      
      for ( Bezier_polygon_with_holes::Hole_const_iterator hit = bpwh.holes_begin() ; hit != bpwh.holes_end() ; ++ hit )
        save_bezier_polygon(out_file, *hit);
      
      rOK = true ;
    }
  }
  
  return rOK ;
  
}

bool save_bezier_sources ( QString aFileName, Bezier_region_source_container const& aSources )
{
  bool rOK = false ;
  
  std::ofstream out_file( qPrintable(aFileName) ) ;
  if ( out_file )
  {
    out_file << std::setprecision(19);
    
    out_file << "DOUBLE" << std::endl ;
    
    out_file << aSources.size() << std::endl ;
    
    for( Bezier_region_source_container::const_iterator rit = aSources.begin(); rit != aSources.end() ; ++ rit )
    {
      Bezier_region_source const& br = *rit ;
      
      out_file << "  " << br.size() << std::endl ;
      
      for( Bezier_region_source::const_iterator bit = br.begin(); bit != br.end() ; ++ bit )
      {
        Bezier_boundary_source const& bb = *bit ;
        
        out_file << "   " << bb.size() << std::endl ;
        
        for ( Bezier_boundary_source::const_iterator cit = bb.begin() ; cit != bb.end() ; ++ cit )
        {
          Bezier_curve const& bc = *cit ;

          out_file << "    " << bc.number_of_control_points() << std::endl ;
          
          for ( Bezier_curve::Control_point_iterator pit = bc.control_points_begin() ; pit != bc.control_points_end() ; ++ pit )
          {
            out_file << "     " << CGAL::to_double(pit->x()) << " " << CGAL::to_double(pit->y()) << std::endl ;
          }
        }
      } 
    }
    
    rOK = true ;
  }
  
  return rOK ;
  
}

void MainWindow::on_actionSaveResult_triggered()
{
  if ( m_circular_active )
  {
    if ( !save_circular(QFileDialog::getSaveFileName(this, tr("Save Result Circular Polygon Set"), "../data", tr("Circular Curve files (*.lps)") ) 
                       ,result_set().circular()
                       )
       )
    {
      show_error("Cannot save circular polygon set.");
    }
       
  }
  else if (m_bezier_active)
  {
    if ( !save_bezier_result(QFileDialog::getSaveFileName(this, tr("Save Result Bezier Polygon Set"), "../data", tr("Bezier Curve files (*.bps)") )
                            ,result_set().bezier() 
                            )
       )
    {
      show_error("Cannot save bezier polygon set.");
    }
  }

  else 
  {
    if ( !save_linear(QFileDialog::getSaveFileName(this, tr("Save Result Linear Polygon Set"), "../data", tr("Linear Curve files (*.lps)") ) 
                       ,result_set().linear()
                       )
       )
    {
      show_error("Cannot save circular polygon set.");
    }
  }


}



//check out
void MainWindow::switch_set_type(Curve_set& aSet, int aType)
{
  unlink_GI(aSet.gi());
  aSet.reset_type(aType);
  link_GI(aSet.gi());
  modelChanged();
}

void MainWindow::switch_sets_type(int aType)
{
  switch_set_type(blue_set(), aType);
  switch_set_type(red_set(), aType);
  switch_set_type(black_set(), aType);
  switch_set_type(brown_set(), aType);
  switch_set_type(yellow_set(), aType);
  switch_set_type(magenta_set(), aType);
  switch_set_type(aqua_set(), aType);
  switch_set_type(result_set(), aType);
  switch_set_type(universal_set(), aType);
}

bool MainWindow::ensure_circular_mode()
{

  if (! m_circular_active) {
    bool lProceed = blue_set().is_empty() && red_set().is_empty() &&
      black_set().is_empty() && brown_set().is_empty() &&
      yellow_set().is_empty() && magenta_set().is_empty() &&
      aqua_set().is_empty() && universal_set().is_empty();

    if (! lProceed)
      lProceed = ask_user_yesno("Circular mode switch",
                                "You are about to load a circular poygon, but there are linear or bezier curves already loaded.\n" \
                                "Both types are not interoperable. In order to proceed, the polygons must be removed first.\n" \
                                "Yes to remove and proceed?\n"
                              );

    if (lProceed) {
      switch_sets_type(2);
      m_circular_active = true;
      m_bezier_active  = false;
    }
  }
  return m_circular_active;
}

bool MainWindow::ensure_bezier_mode()
{
  if (! m_bezier_active )
  {
    bool lProceed = blue_set().is_empty() && red_set().is_empty() &&
      black_set().is_empty() && brown_set().is_empty() &&
      yellow_set().is_empty() && magenta_set().is_empty() &&
      aqua_set().is_empty() && universal_set().is_empty();
    
    if ( ! lProceed )
      lProceed = ask_user_yesno("Bezier mode switch"
                               ,"You are about to load a Bezier curve, but there are linear and/or circular polygons already loaded.\n" \
                                "Both types are not interoperable. In order to proceed, the polygons must be removed first.\n" \
                                "Yes to remove and proceed?\n"
                               ) ;
      
    if ( lProceed )
    {
      switch_sets_type(3);
      m_bezier_active = true;
      m_circular_active = false;
    }
  }
  return m_bezier_active ;
}

bool MainWindow::ensure_linear_mode()
{
  if (m_circular_active || m_bezier_active) {
    bool lProceed = blue_set().is_empty() && red_set().is_empty() &&
      black_set().is_empty() && brown_set().is_empty() &&
      yellow_set().is_empty() && magenta_set().is_empty() &&
      aqua_set().is_empty() && universal_set().is_empty();

    if (! lProceed)
      lProceed = ask_user_yesno("Linear/Circular mode switch",
                                "You are about to load a linear poygon, but there are circular curves already loaded.\n" \
                                "Both types are not interoperable. In order to proceed, the polygons must be removed first.\n" \
                                "Yes to remove and proceed?\n"
                              );

    if (lProceed) {
      switch_sets_type(1);
      m_circular_active = false;
      m_bezier_active = false;
    }
  }
  return !m_circular_active;
}

//check out
//bool read_linear(QString /* aFileName */, Linear_polygon_set& /* rSet */,
//                 Linear_region_source_container& /* rSources */)
//{
//  bool rOK = false;
//  return rOK;
//}

//bool read_circular(QString /* aFileName */, Circular_polygon_set& /* rSet */,
 //                  Circular_region_source_container& /* rSources */)
//{
//  bool rOK = false;
//}

void MainWindow::open(QString fileName)
{
  /*if(! fileName.isEmpty()) {
    bool lRead = false;
    if(fileName.endsWith(".lps"))
    {
      if ( ensure_circular_mode() )
        lRead = read_linear(fileName,active_set().circular(), active_circular_sources() ) ;
    }
    else if (fileName.endsWith(".bps"))
    {
      if ( ensure_bezier_mode() )
        lRead = read_bezier(fileName,active_set().bezier(), active_bezier_sources() ) ;
    }
    if (lRead) {
      modelChanged();
      zoomToFit();
      this->addToRecentFiles(fileName);
    }
  }*/
}



void MainWindow::on_actionInsertConicCircle_triggered()
{}


void MainWindow::on_actionInsertConicEclipse_triggered()
{}


void MainWindow::on_actionInsertCircular_toggled(bool aChecked)
{
	if(aChecked)
	{
	  this->graphicsView->setDragMode(QGraphicsView::NoDrag);
	  if(ensure_circular_mode()) 
	  {
	  	actionPAN->setChecked(false);
	  	actionInsertLinear->setChecked( false );
		actionInsertBezier->setChecked( false ); 
	  	m_scene.installEventFilter(m_circular_input);
	  }
	}
}

void MainWindow::on_actionInsertBezier_toggled(bool aChecked)
{
	if(aChecked)
	{
	  this->graphicsView->setDragMode(QGraphicsView::NoDrag);
	  if(ensure_bezier_mode()) 
	  {	
	  	actionPAN->setChecked(false);
	  	actionInsertLinear->setChecked( false );
		actionInsertCircular->setChecked( false );
	  	m_scene.installEventFilter(m_bezier_input);
	  }
	}
}

void MainWindow::on_actionInsertLinear_toggled(bool aChecked)
{
	if(aChecked)
	{
	  this->graphicsView->setDragMode(QGraphicsView::NoDrag);
	  if (ensure_linear_mode())
	  {
	  	actionPAN->setChecked(false);
	  	actionInsertCircular->setChecked( false );
		actionInsertBezier->setChecked( false ); 
	  	m_scene.installEventFilter(m_linear_input);
	  }
	}
}

void MainWindow::on_actionMinkowski_Sum_triggered()
{
  bool lDone = false;
  QCursor old = this->cursor();
  this->setCursor(Qt::WaitCursor);

  actionComplement->setChecked(false);
  actionUnion->setChecked(false);
  actionIntersection->setChecked(false);
  actionDifference->setChecked(false); 
  actionSymmetric_Difference->setChecked(false); 
  //actionMinkowski_Sum->setChecked(false);

  if (!blue_set().is_empty()) lDone = true;
  this->setCursor(old);
  if (lDone) modelChanged();
}

//only blue complement


void MainWindow::on_actionComplement_triggered()
{
  bool lDone = false;
  QCursor old = this->cursor();
  this->setCursor(Qt::WaitCursor);
  //actionComplement->setChecked(false);
  actionUnion->setChecked(false);
  actionIntersection->setChecked(false);
  actionDifference->setChecked(false); 
  actionSymmetric_Difference->setChecked(false); 
  actionMinkowski_Sum->setChecked(false);

  if(!universal_set().is_empty()) 
  {
	  switch(m_color_complement)
	  {
	    case 0: if(!blue_set().is_empty()) {result_set().assign(universal_set()); result_set().difference(blue_set());} break;
	    case 1: if(!red_set().is_empty()) {result_set().assign(universal_set()); result_set().assign(result_set());result_set().difference(red_set());} break;
	    case 2: if(!black_set().is_empty()) {result_set().assign(universal_set());  result_set().assign(result_set());result_set().difference(black_set());} break;
	    case 3: if(!brown_set().is_empty()) {result_set().assign(universal_set());  result_set().assign(result_set());result_set().difference(brown_set());} break;
	    case 4: if(!yellow_set().is_empty()) {result_set().assign(universal_set());  result_set().assign(result_set());result_set().difference(yellow_set());} break;
	    case 5: if(!magenta_set().is_empty()) {result_set().assign(universal_set()); result_set().assign(result_set());result_set().difference(magenta_set());} break;
	    case 6: if(!aqua_set().is_empty()) {result_set().assign(universal_set()); result_set().assign(result_set());result_set().difference(aqua_set());} break;
	  }
  }

  else
  {
  	bool lProceed = ask_user_yesno("Complement Operation Error", "Universal Polygon Set Not Defined\nWould like to define Universal Polygon");
  }

  lDone = true;
  this->setCursor(old);
  if (lDone) modelChanged();
}

void MainWindow::on_actionIntersection_triggered()
{
  bool lDone = false;
  QCursor old = this->cursor();
  this->setCursor(Qt::WaitCursor);

  actionComplement->setChecked(false);
  actionUnion->setChecked(false);
  //actionIntersection->setChecked(false);
  actionDifference->setChecked(false); 
  actionSymmetric_Difference->setChecked(false); 
  actionMinkowski_Sum->setChecked(false);

  if (!blue_set().is_empty() && m_blue_int) result_set().assign(blue_set());
  else if (!red_set().is_empty() && m_red_int) result_set().assign(red_set());
  else if (!black_set().is_empty() && m_black_int) result_set().assign(black_set());
  else if (!brown_set().is_empty() && m_brown_int) result_set().assign(brown_set());
  else if (!yellow_set().is_empty() && m_yellow_int) result_set().assign(yellow_set());
  else if (!magenta_set().is_empty() && m_magenta_int) result_set().assign(magenta_set());
  else if (!aqua_set().is_empty() && m_aqua_int) result_set().assign(aqua_set());

  if (!blue_set().is_empty() && m_blue_int) result_set().intersect(blue_set());
  if (!red_set().is_empty() && m_red_int) result_set().intersect(red_set());
  if (!black_set().is_empty() && m_black_int) result_set().intersect(black_set());
  if (!brown_set().is_empty() && m_brown_int) result_set().intersect(brown_set());
  if (!yellow_set().is_empty() && m_yellow_int) result_set().intersect(yellow_set());
  if (!magenta_set().is_empty() && m_magenta_int) result_set().intersect(magenta_set());
  if (!aqua_set().is_empty() && m_aqua_int) result_set().intersect(aqua_set());
    
  lDone = true;

  this->setCursor(old);
  if (lDone) modelChanged();
}

//blue - red
//no idea
void MainWindow::on_actionDifference_triggered()
{
  bool lDone = false;
  QCursor old = this->cursor();
  this->setCursor(Qt::WaitCursor);

  actionComplement->setChecked(false);
  actionUnion->setChecked(false);
  actionIntersection->setChecked(false);
  //actionDifference->setChecked(false); 
  actionSymmetric_Difference->setChecked(false); 
  actionMinkowski_Sum->setChecked(false);

  /*switch(m_color_active) {
    case 0: result_set().assign(red_set());result_set().join(black_set());result_set().join(brown_set());result_set().join(magenta_set());result_set().join(yellow_set());result_set().join(aqua_set());if(!result_set().is_empty())result_set().difference(blue_set());break;
    case 1: result_set().assign(blue_set());result_set().join(black_set());result_set().join(brown_set());result_set().join(magenta_set());result_set().join(yellow_set());result_set().join(aqua_set());if(!result_set().is_empty())result_set().difference(red_set());break;
    case 2: result_set().assign(red_set());result_set().join(blue_set());result_set().join(brown_set());result_set().join(magenta_set());result_set().join(yellow_set());result_set().join(aqua_set());if(!result_set().is_empty())result_set().difference(black_set());break;
    case 3: result_set().assign(red_set());result_set().join(black_set());result_set().join(blue_set());result_set().join(magenta_set());result_set().join(yellow_set());result_set().join(aqua_set());if(!result_set().is_empty())result_set().difference(brown_set());break;
    case 4: result_set().assign(red_set());result_set().join(black_set());result_set().join(brown_set());result_set().join(magenta_set());result_set().join(blue_set());result_set().join(aqua_set());if(!result_set().is_empty())result_set().difference(yellow_set());break;
    case 5: result_set().assign(red_set());result_set().join(black_set());result_set().join(brown_set());result_set().join(blue_set());result_set().join(yellow_set());result_set().join(aqua_set());if(!result_set().is_empty())result_set().difference(magenta_set());break;
    case 6: result_set().assign(red_set());result_set().join(black_set());result_set().join(brown_set());result_set().join(magenta_set());result_set().join(yellow_set());result_set().join(blue_set());if(!result_set().is_empty())result_set().difference(aqua_set());break;
    /*if (!blue_set().is_empty() && !red_set().is_empty()) {
      result_set().assign(blue_set());
      result_set().difference(red_set());
      lDone = true;
      case 0: blue_set().assign(result_set()); break;
     case 1: red_set().assign(result_set()); break;
     case 2: black_set().assign(result_set()); break;
     case 3: brown_set().assign(result_set()); break;
     case 4: yellow_set().assign(result_set()); break;
     case 5: magenta_set().assign(result_set()); break;
     case 6: aqua_set().assign(result_set()); break;
  }*/


  bool lProceed;
  if(m_color_diff_A == m_color_diff_B)
  {
   	lProceed =  ask_user_yesno("Difference Operation Error", "Single colorued polygon set is not supported by this operation. If you want to proceed further, please select one more coloured polygon set.\nWould you select one more color polygon set?\n");
  }
  else
  {
  		if(m_color_diff_A == 0)
  		{
  			result_set().assign(blue_set());
  		}

  		if(m_color_diff_A == 1)
  		{
  			result_set().assign(red_set());
  		}

  		if(m_color_diff_A == 2)
  		{
  			result_set().assign(black_set());
  		}

  		if(m_color_diff_A == 3)
  		{
  			result_set().assign(brown_set());
  		}


  		if(m_color_diff_A == 4)
  		{
  			result_set().assign(yellow_set());
  		}

  		if(m_color_diff_A == 5)
  		{
  			result_set().assign(magenta_set());
  		}

  		if(m_color_diff_B == 1)
  		{
  			result_set().difference(red_set());
  		}

  		if(m_color_diff_B == 2)
  		{
  			result_set().difference(black_set());
  		}


  		if(m_color_diff_B == 3)
  		{
  			result_set().difference(brown_set());
  		}

  		if(m_color_diff_B == 4)
  		{
  			result_set().difference(yellow_set());
  		}

  		if(m_color_diff_B == 5)
  		{
  			result_set().difference(magenta_set());
  		}

  		if(m_color_diff_B == 6)
  		{
  			result_set().difference(aqua_set());
  		}
  }
  lDone = true;
  this->setCursor(old);
  if (lDone) modelChanged();
}

void MainWindow::on_actionSymmetric_Difference_triggered()
{
  bool lDone = false;
  QCursor old = this->cursor();
  this->setCursor(Qt::WaitCursor);

  actionComplement->setChecked(false);
  actionUnion->setChecked(false);
  actionIntersection->setChecked(false);
  actionDifference->setChecked(false); 
  //actionSymmetric_Difference->setChecked(false); 
  actionMinkowski_Sum->setChecked(false);

  if (!blue_set().is_empty() && m_blue_sym_diff) result_set().assign(blue_set());
  else if (!red_set().is_empty() && m_red_sym_diff) result_set().assign(red_set());
  else if (!black_set().is_empty() && m_black_sym_diff) result_set().assign(black_set());
  else if (!brown_set().is_empty() && m_brown_sym_diff) result_set().assign(brown_set());
  else if (!yellow_set().is_empty() && m_yellow_sym_diff) result_set().assign(yellow_set());
  else if (!magenta_set().is_empty() && m_magenta_sym_diff) result_set().assign(magenta_set());
  else if (!aqua_set().is_empty() && m_aqua_sym_diff) result_set().assign(aqua_set());

  if (!red_set().is_empty() && m_red_sym_diff) result_set().symmetric_difference(red_set());
  if (!black_set().is_empty() && m_black_sym_diff) result_set().symmetric_difference(black_set());
  if (!brown_set().is_empty() && m_brown_sym_diff) result_set().symmetric_difference(brown_set());
  if (!yellow_set().is_empty() && m_yellow_sym_diff) result_set().symmetric_difference(yellow_set());
  if (!magenta_set().is_empty() && m_magenta_sym_diff) result_set().symmetric_difference(magenta_set());
  if (!aqua_set().is_empty() && m_aqua_sym_diff) result_set().symmetric_difference(aqua_set());
  lDone = true;
  this->setCursor(old);
  if (lDone) modelChanged();
}

void MainWindow::on_actionUnion_triggered()
{
  bool lDone = false;
  QCursor old = this->cursor();
  this->setCursor(Qt::WaitCursor);

  actionComplement->setChecked(false);
  //actionUnion->setChecked(false);
  actionIntersection->setChecked(false);
  actionDifference->setChecked(false); 
  actionSymmetric_Difference->setChecked(false); 
  actionMinkowski_Sum->setChecked(false);

  result_set().clear();

  if(m_red_union) result_set().assign(red_set());
  if(m_blue_union) result_set().join(blue_set());
  if(m_black_union) result_set().join(black_set());
  if(m_brown_union) result_set().join(brown_set());
  if(m_magenta_union) result_set().join(magenta_set());
  if(m_yellow_union) result_set().join(yellow_set());
  if(m_aqua_union) result_set().join(aqua_set());
  lDone = true;
  this->setCursor(old);

  if (lDone) modelChanged();
}

//to change which polygons to see on the screen
void MainWindow::ToogleView(size_t aGROUP, bool a_check)
{
  if (a_check) set(aGROUP).gi()->show();
  else set(aGROUP).gi()->hide();
}

void MainWindow::on_actionPAN_toggled(bool aChecked)
{
	if(aChecked)
	{
	  if (!m_circular_active && !m_bezier_active) 
	  	{
	  		m_scene.removeEventFilter(m_linear_input); 
	  		m_scene.removeEventFilter(m_bezier_input);
	  		m_scene.removeEventFilter(m_circular_input);
	  		actionInsertLinear->setChecked( false );
	  		this->graphicsView->setDragMode(QGraphicsView::ScrollHandDrag); 
	  		//m_scene.installEventFilter(m_linear_input);
	  	}
	  else if(!m_bezier_active) 
	  	{ 
	  		m_scene.removeEventFilter(m_linear_input); 
	  		m_scene.removeEventFilter(m_bezier_input);
	  		m_scene.removeEventFilter(m_circular_input);
	  		actionInsertCircular->setChecked( false );
	  		this->graphicsView->setDragMode(QGraphicsView::ScrollHandDrag); 
	  		//m_scene.installEventFilter(m_circular_input);
	  	}
	  else 
	  	{ 	
	  		m_scene.removeEventFilter(m_linear_input); 
	  		m_scene.removeEventFilter(m_bezier_input);
	  		m_scene.removeEventFilter(m_circular_input);
	  		actionInsertBezier->setChecked( false );
	  		this->graphicsView->setDragMode(QGraphicsView::ScrollHandDrag); 
	  		//m_scene.installEventFilter(m_bezier_input);
	  	}
	}
  
}

void MainWindow::zoomToFit()
{
  boost::optional<QRectF> lTotalRect;

  for (auto si = m_curve_sets.begin(); si != m_curve_sets.end(); ++ si) 
  {
    if (!si->is_empty()) 
    {
      QRectF lRect = si->bounding_rect();
      if (lTotalRect) lTotalRect = *lTotalRect | lRect;
      else lTotalRect = lRect;
    }
  }

  if (lTotalRect) {
    this->graphicsView->setSceneRect(*lTotalRect);
    this->graphicsView->fitInView(*lTotalRect, Qt::KeepAspectRatio);
  }
}


void MainWindow::processInput(CGAL::Object o)
{
    std::pair<Bezier_polygon,Bezier_boundary_source>     lBI ;
    Linear_polygon lLI;
    Circular_polygon lCI; 

    if(CGAL::assign(lBI, o))
    {
      if ( ensure_bezier_mode() )
      {
        CGAL::Orientation o = lBI.first.orientation();
        if ( o == CGAL::CLOCKWISE )
          lBI.first.reverse_orientation();
          
        active_set().bezier().join( Bezier_polygon_with_holes(lBI.first) ) ;  
        
        Bezier_region_source br ; br.push_back (lBI.second);
        
        active_bezier_sources().push_back(br);
        
      }
    }

    if (CGAL::assign(lLI, o)) 
    {
      if (ensure_linear_mode()) {
        CGAL::Orientation orient = lLI.orientation();
        if (orient == CGAL::CLOCKWISE) 
          lLI.reverse_orientation();
        Linear_polygon_with_holes lCPWH(lLI);
        active_set().linear().join(lCPWH);
        active_linear_sources().push_back(lCPWH);
    }
  }



    else if ( CGAL::assign(lCI, o) )
    {
      if ( ensure_circular_mode() )
      {
        CGAL::Orientation o = lCI.orientation();
        if ( o == CGAL::CLOCKWISE )
          lCI.reverse_orientation();
          
        Circular_polygon_with_holes lCPWH(lCI);
        active_set().circular().join(lCPWH) ;  
        
        active_circular_sources().push_back(lCPWH);
      }
    }
  modelChanged();  


}


//Main part
#include "Boolean_set_operations_2.moc"
#include <CGAL/Qt/resources.h>
int main(int argc, char* argv[])
{
  //QApplication a(argc, argv);
  QApplication app(argc, argv);

  app.setOrganizationDomain("geometryfactory.com");
  app.setOrganizationName("GeometryFactory");
  app.setApplicationName("Boolean_operations_2 demo");
  CGAL_QT_INIT_RESOURCES;
  try
  {
    MainWindow w;
    w.show();
    return app.exec();
  }
  catch (const std::exception& e)
  {
    std::string s = e.what();
    show_error("Exception throne during run of the program:\n" + s);
  }
}
