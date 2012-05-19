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
// $URL: svn+ssh://scm.gforge.inria.fr/svn/cgal/branches/features/gsoc2012-Arrangement_on_surface_2-demo-atsui/Arrangement_on_surface_2/demo/Arrangement_on_surface_2/forms.cpp $
// $Id: forms.cpp 67117 2012-01-13 18:14:48Z lrineau $
//
//
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>

#include <CGAL/basic.h>


////////////////////////////////////////////////////////////////////////

#include "forms.h"
#include "demo_tab.h"

/*! constructor - build the properties dialog form */
PropertiesForm::PropertiesForm(QTabWidget * bar, QWidget* parent,
                               int /* number_of_tabs */,
                               Qt_widget_base_tab * w_demo_p, double scale,
                               bool colors_flag):
  QDialog(parent),
  myBar(bar)
{
  setCaption( "Properties -- Options" );
  resize( 420, 290 );

  optionsFormLayout = new QVBoxLayout( this, 11, 6 );

  arrLayout1 = new QHBoxLayout( 0, 0, 6 );
  textLabel1 = new QLabel( "Width", this );
  arrLayout1->addWidget( textLabel1 );
  box1 = new QSpinBox( 300, 1000, 50, this, "box1" );
  box1->setValue(parent->width());
  arrLayout1->addWidget( box1 );
  optionsFormLayout->addLayout( arrLayout1 );

  arrLayout2 = new QHBoxLayout( 0, 0, 6 );
  textLabel2 = new QLabel( "Height", this );
  arrLayout2->addWidget( textLabel2 );
  box2 = new QSpinBox( 300, 1000, 50, this, "box2" );
  box2->setValue(parent->height());
  arrLayout2->addWidget( box2 );
  optionsFormLayout->addLayout( arrLayout2 );

  arrLayout3 = new QHBoxLayout( 0, 0, 6 );
  textLabel3 = new QLabel( "Line Width", this );
  arrLayout3->addWidget( textLabel3 );
  box3 = new QSpinBox( 1, 5, 1, this, "box3" );
  box3->setValue(w_demo_p->m_line_width);
  arrLayout3->addWidget( box3 );
  optionsFormLayout->addLayout( arrLayout3 );

  arrLayout4 = new QHBoxLayout( 0, 0, 6 );
  textLabel4 = new QLabel( "Scaling Factor", this );
  arrLayout4->addWidget( textLabel4 );
  box4 = new MySpinBox( 10, 100, 1, this, "box4" );
  box4->setValue(static_cast<int>(scale*10));
  arrLayout4->addWidget( box4 );
  optionsFormLayout->addLayout( arrLayout4 );

  arrLayout5 = new QHBoxLayout( 0, 0, 6 );
  textLabel5 = new QLabel( "Display Mode", this );
  arrLayout5->addWidget( textLabel5 );
  box5 = new QComboBox( FALSE, this );
  box5->insertItem( "Different Colors At Overlay" );
  box5->insertItem( "Uniform Color At Overlay" );
  arrLayout5->addWidget( box5 );
  optionsFormLayout->addLayout( arrLayout5 );
  if (!colors_flag)
    box5->setCurrentItem(1);

  arrLayout6 = new QHBoxLayout( 0, 0, 6 );
  textLabel6 = new QLabel( "Grid Cube Size", this );
  arrLayout6->addWidget( textLabel6 );
  box6 = new QSpinBox( 1, 100, 1, this, "box6" );
  box6->setValue(w_demo_p->cube_size);
  arrLayout6->addWidget( box6 );
  optionsFormLayout->addLayout( arrLayout6 );

  arrLayout7 = new QHBoxLayout( 0, 0, 6 );
  textLabel7 = new QLabel( "Remove Curve Mode", this );
  arrLayout7->addWidget( textLabel7 );
  box7 = new QComboBox( FALSE, this );
  box7->insertItem( "Remove entire original curve" );
  box7->insertItem( "Remove Edge" );
  arrLayout7->addWidget( box7 );
  optionsFormLayout->addLayout( arrLayout7 );
  if (!w_demo_p->remove_org_curve)
    box7->setCurrentItem(1);

  arrLayout8 = new QHBoxLayout( 0, 0, 6 );
  textLabel8 = new QLabel( "Vertex Radius", this );
  arrLayout8->addWidget( textLabel8 );
  box8 = new QSpinBox( 1, 5, 1, this, "box8" );
  box8->setValue(w_demo_p->m_vertex_width);
  arrLayout8->addWidget( box8 );
  optionsFormLayout->addLayout( arrLayout8 );

  arrLayout9 = new QHBoxLayout( 0, 0, 6 );
  textLabel9 = new QLabel( "Draw vertex not in intersection", this );
  arrLayout9->addWidget( textLabel9 );
  box9 = new QComboBox( FALSE, this );
  box9->insertItem( "Draw" );
  box9->insertItem( "Don't draw" );
  arrLayout9->addWidget( box9 );
  optionsFormLayout->addLayout( arrLayout9 );
  if (!w_demo_p->draw_vertex)
    box9->setCurrentItem(1);

  buttonsLayout = new QHBoxLayout( 0, 0, 6 );
  okPushButton = new QPushButton( "OK", this );
  okPushButton->setDefault( TRUE );
  buttonsLayout->addWidget( okPushButton );
  cancelPushButton = new QPushButton( "Cancel", this );
  buttonsLayout->addWidget( cancelPushButton );
  optionsFormLayout->addLayout( buttonsLayout );
  connect( okPushButton, SIGNAL( clicked() ), this, SLOT( accept() ) );
  connect( cancelPushButton, SIGNAL( clicked() ), this, SLOT( reject() ) );

  textLabel1->setBuddy( box1 );
  textLabel2->setBuddy( box2 );
  textLabel3->setBuddy( box3 );
  textLabel4->setBuddy( box4 );

}

////////////////////////////////////////////////////////////////////////

/*! greem_icon - used in the overlay form */
const char* green_icon[]={
  "16 16 2 1",
  "g c green",
  ". c None",
  "................",
  "................",
  "..gggggggggggg..",
  "..gggggggggggg..",
  "..gggggggggggg..",
  "..ggg......ggg..",
  "..ggg......ggg..",
  "..ggg......ggg..",
  "..ggg......ggg..",
  "..ggg......ggg..",
  "..ggg......ggg..",
  "..gggggggggggg..",
  "..gggggggggggg..",
  "..gggggggggggg..",
  "................",
  "................"};

/*! white_icon - used in the overlay form */
const char* white_icon[]={
  "16 16 2 1",
  "g c green",
  ". c None",
  "................",
  "................",
  "................",
  "................",
  "................",
  "................",
  "................",
  "................",
  "................",
  "................",
  "................",
  "................",
  "................",
  "................",
  "................",
  "................"};

/*! OverlayForm constructor - build the overlay dialog form */
OverlayForm::OverlayForm(  QTabWidget * bar, QWidget* parent ,int tab_number ,
                           const char* name, bool modal, WFlags f  ):
  QDialog( parent, name, modal, f ),
  myBar(bar)
{
  setCaption( "Planar Maps  --  Overlay" );
  resize( 590, 390 );

  optionsFormLayout = new QVBoxLayout( this, 11, 6 );

  split = new QSplitter(this);
  listBox1 = new DDListBox( split );
  listBox2 = new DDListBox( split );
  QString traits;
  Qt_widget_base_tab    *w_demo_p;
  for (int i=0; i < tab_number; i++)
  {
    if ( myBar->isTabEnabled( myBar->page(i) ) )
    {
      // We peform downcasting from QWigdet* to Qt_widget_base_tab*,
      // as we know that only
      // Qt_widget_base_tab objects are stored in the tab pages.
      w_demo_p = static_cast<Qt_widget_base_tab *> (myBar->page(i));
      switch ( w_demo_p->traits_type ) {
       case SEGMENT_TRAITS:
        traits = " ( segment traits )";
        break;
       case POLYLINE_TRAITS:
        traits = " ( polyline traits )";
        break;
       case CONIC_TRAITS:
        traits = " ( conic traits )";
        break;
      }
      listBox1->insertItem( QPixmap( green_icon ) ,
                            myBar->label(i) + traits );
    }
  }
  listBox1->set_max_items(listBox1->count());

  arrLayout = new QHBoxLayout();

  textLabel1 = new QLabel( "Possible Planar Maps", this );
  arrLayout->addWidget( textLabel1 );

  textLabel2 = new QLabel( "Chosen Planar Maps", this );
  arrLayout->addWidget( textLabel2 );

  buttonsLayout = new QHBoxLayout( 0, 0, 6 );
  okPushButton = new QPushButton( "OK", this );
  okPushButton->setDefault( TRUE );
  buttonsLayout->addWidget( okPushButton );

  cancelPushButton = new QPushButton( "Cancel", this );
  buttonsLayout->addWidget( cancelPushButton );

  optionsFormLayout->addLayout( arrLayout );
  optionsFormLayout->addWidget( split );
  optionsFormLayout->addLayout( buttonsLayout );

  connect( okPushButton, SIGNAL( clicked() ), this, SLOT( accept() ) );
  connect( cancelPushButton, SIGNAL( clicked() ), this, SLOT( reject() ) );

  setAcceptDrops(TRUE);
}
////////////////////////////////////////////////////////////////////////////////
//CheckItem::CheckItem(  QListBox * listbox, const QPixmap & pix,
//                                                               const QString & text ):
//  QListBoxPixmap( listbox, pix, text )
//{
//  check_box = new QCheckBox(listbox);
//}


//////////////////////////////////////////////////////////////////////////////

/*! DDListBox constructor */
DDListBox::DDListBox( QWidget * parent, const char * name, WFlags f ) :
  QListBox( parent, name, f ),
  max_items(0),
  flag(false)
{
  setAcceptDrops( TRUE );
  dragging = FALSE;
}

/*! dragEnterEvent - accept drag event */
void DDListBox::dragEnterEvent( QDragEnterEvent *evt )
{
  if (  QTextDrag::canDecode( evt ) )
    evt->accept();
}

/*! dropEvent - insert new item. if we are in the first listBox
 *  and all items are there, make all selectable
 */
void DDListBox::dropEvent( QDropEvent *evt )
{
  QString text;
  if (  QTextDrag::decode( evt, text ) )
    insertItem( QPixmap( green_icon ) , text );
  if (count() == max_items && max_items != 0)
  {
    flag = true;
    for (unsigned int i = 0; i < count(); i++)
      item(i)->setSelectable( true );
  }
}

/*! mousePressEvent - mouse click on the list box */
void DDListBox::mousePressEvent( QMouseEvent *evt )
{
  QListBox::mousePressEvent( evt );
  dragging = TRUE;
}

/*! mouseMoveEvent - mouse move on the list box */
void DDListBox::mouseMoveEvent( QMouseEvent * )
{
  if (count() == max_items && max_items != 0 && flag)
  {
    for (unsigned int i = 0; i < count(); i++)
      changeItem( QPixmap( green_icon ) , text(i) , i);
    flag = false;
  }

  if ( dragging && item(currentItem())->isSelectable() )
  {
    QDragObject *d = new  QTextDrag( currentText() , this );
    d->dragCopy(); // do NOT delete d.
    dragging = FALSE;
    unsigned int current = currentItem();
    if (count() == max_items && max_items != 0)
    {
      char s[100];
      strcpy(s, currentText());
      char * traits;
      traits = strtok(s," ");
      traits = strtok(NULL, " ");
      traits = strtok(NULL, " ");
      traits = strtok(NULL, " ");

      for (unsigned int i = 0; i < max_items; i++)
      {
        char s_i[100];
        strcpy(s_i, text(i));
        char * traits_i;
        traits_i = strtok(s_i," ");
        traits_i = strtok(NULL, " ");
        traits_i = strtok(NULL, " ");
        traits_i = strtok(NULL, " ");
        bool b = (strcmp(traits,traits_i) == 0);
        if (!b && i != current )
        {
          changeItem( QPixmap( white_icon ) , text(i) , i);
          item(i)->setSelectable( b );
        }
      }
    }

    removeItem ( current );
  }
}

/*! set_max_items - access to private date member */
void DDListBox::set_max_items(int num)
{
  max_items = num;
}

/*! OptionsForm constructor */
OptionsForm::OptionsForm(QWidget * parent, int /* number_of_tabs */,
                         const char * name, bool modal, WFlags f):
  QDialog( parent, name, modal, f )
{
  setCaption( "Conic Type - Options" );
  resize( 320, 290 );

  optionsFormLayout = new QVBoxLayout( this, 11, 6 );

  arrLayout1 = new QHBoxLayout( 0, 0, 6 );

  textLabel1 = new QLabel( "Conic Type", this );
  arrLayout1->addWidget( textLabel1 );

  arrComboBox1 = new QComboBox( FALSE, this );

  arrComboBox1->insertItem( "Circle" );
  arrComboBox1->insertItem( "Segment" );
  arrComboBox1->insertItem( "Ellipse" );
  arrComboBox1->insertItem( "Parabula" );
  arrComboBox1->insertItem( "Hyperbula" );

  arrLayout1->addWidget( arrComboBox1 );
  optionsFormLayout->addLayout( arrLayout1 );

  buttonsLayout = new QHBoxLayout( 0, 0, 6 );
  okPushButton = new QPushButton( "OK", this );
  okPushButton->setDefault( TRUE );
  buttonsLayout->addWidget( okPushButton );
  cancelPushButton = new QPushButton( "Cancel", this );
  buttonsLayout->addWidget( cancelPushButton );
  optionsFormLayout->addLayout( buttonsLayout );

  connect( okPushButton, SIGNAL( clicked() ), this, SLOT( accept() ) );
  connect( cancelPushButton, SIGNAL( clicked() ), this, SLOT( reject() ) );

  textLabel1->setBuddy( arrComboBox1 );

}

/*! CheckForm constructor */
CheckForm::CheckForm( OverlayForm *overlay_form , QWidget* parent ):
  QDialog( parent )
{
  setCaption( "Overlay - paint intersections" );
  resize( 320, 290 );

  optionsFormLayout = new QVBoxLayout( this, 11, 6 );
  layout = new QHBoxLayout( 0, 0, 6 );
  button_group =
          new QVButtonGroup("Check to paint Planar Maps intersections", this);

  for (unsigned int i = 0; i < overlay_form->listBox2->count(); i++)
  {
    overlay_form->listBox2->setCurrentItem(i);
    QCheckBox *b = new QCheckBox(overlay_form->listBox2->currentText() , button_group);
        b->setChecked( true );
        button_group->insert( b , i );
  }

  layout->addWidget( button_group );
  optionsFormLayout->addLayout( layout );

  buttonsLayout = new QHBoxLayout( 0, 0, 6 );
  okPushButton = new QPushButton( "OK", this );
  okPushButton->setDefault( TRUE );
  buttonsLayout->addWidget( okPushButton );
  cancelPushButton = new QPushButton( "Cancel", this );
  buttonsLayout->addWidget( cancelPushButton );
  optionsFormLayout->addLayout( buttonsLayout );

  connect( okPushButton, SIGNAL( clicked() ), this, SLOT( accept() ) );
  connect( cancelPushButton, SIGNAL( clicked() ), this, SLOT( reject() ));

}

/*! FileOpenOptionsForm constructor */
FileOpenOptionsForm::FileOpenOptionsForm( bool flag ,QWidget* parent ,
                                     const char* name, bool modal, WFlags f ):
  QDialog( parent, name, modal, f )
{
  setCaption( "File Open - Options" );
  resize( 320, 290 );

  optionsFormLayout = new QVBoxLayout( this, 11, 6 );

  buttonGroup = new QButtonGroup( 3, Qt::Vertical ,"Do you want to:",
                                  this, "buttonGroup" );

  b1 = new QRadioButton( buttonGroup, "b1");
  b1->setText( "open file in a new tab" );

  b2 = new QRadioButton( buttonGroup, "b2");
  b2->setText( "open file in current tab (delete current Pm)" );

  if (flag)
  {
    b3 = new QRadioButton( buttonGroup, "b3");
    b3->setText( "merge file into current tab" );
  }

  buttonGroup->setButton(0);

  optionsFormLayout->addWidget( buttonGroup );

  buttonsLayout = new QHBoxLayout( 0, 0, 6 );
  okPushButton = new QPushButton( "OK", this );
  okPushButton->setDefault( TRUE );
  buttonsLayout->addWidget( okPushButton );
  cancelPushButton = new QPushButton( "Cancel", this );
  buttonsLayout->addWidget( cancelPushButton );
  optionsFormLayout->addLayout( buttonsLayout );

  connect( okPushButton, SIGNAL( clicked() ), this, SLOT( accept() ) );
  connect( cancelPushButton, SIGNAL( clicked() ), this, SLOT( reject() ) );

}


/*! PointLocationStrategyForm constructor */
PointLocationStrategyForm::PointLocationStrategyForm(QWidget * parent ,
                                                     int /* number_of_tabs */,
                                                     const char * name,
                                                     bool modal, WFlags f):
QDialog( parent, name, modal, f )
{
  setCaption( "Point Location - Strategy" );
  optionsFormLayout = new QVBoxLayout( this, 11, 6 );
  arrLayout1 = new QHBoxLayout( 0, 0, 6 );

  textLabel1 = new QLabel( "Strategy", this );
  arrLayout1->addWidget( textLabel1 );

  arrComboBox1 = new QComboBox( FALSE, this );

  arrComboBox1->insertItem( "Simple" );
  arrComboBox1->insertItem( "Land marks" );
  arrComboBox1->insertItem( "Trapezoiedal" );
  arrComboBox1->insertItem( "Walk" );

  arrLayout1->addWidget( arrComboBox1 );
  optionsFormLayout->addLayout( arrLayout1 );

  buttonsLayout = new QHBoxLayout( 0, 0, 6 );
  okPushButton = new QPushButton( "OK", this );
  okPushButton->setDefault( TRUE );
  buttonsLayout->addWidget( okPushButton );
  cancelPushButton = new QPushButton( "Cancel", this );
  buttonsLayout->addWidget( cancelPushButton );
  optionsFormLayout->addLayout( buttonsLayout );

  connect( okPushButton, SIGNAL( clicked() ), this, SLOT( accept() ) );
  connect( cancelPushButton, SIGNAL( clicked() ), this, SLOT( reject() ) );

  textLabel1->setBuddy( arrComboBox1 );

}


#include "forms.moc"

