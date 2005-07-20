#ifndef FORMS_H
#define FORMS_H

/*! the forms.h and forms.C files contains all the program 
 *  dialog forms and drag-drop class.
 */
#include <qdialog.h>
#include <qcombobox.h>
#include <qlabel.h>  
#include <qlayout.h>  
#include <qpushbutton.h> 
#include <qinputdialog.h> 
#include <qslider.h> 
#include <qspinbox.h> 
#include <qtabwidget.h>
#include <qdragobject.h> 
#include <qlistbox.h> 
#include <qiconview.h> 
#include <qsplitter.h> 
#include <qbuttongroup.h>  
#include <qbutton.h>  
#include <qradiobutton.h>  
#include <qcheckbox.h> 
#include <qvbuttongroup.h> 

#include "cgal_types.h"

class QDragEnterEvent; 
class QDragDropEvent; 
class DDListBox; 
class QWidget;

/*! class MySpinBox is used for the scaling factor properties 
 *  dialog. we need to use it to support duoble as a value. 
 */ 
class MySpinBox : public QSpinBox 
{ 
  Q_OBJECT 
public: 
  /*! constructor */ 
  MySpinBox ( int minValue, int maxValue, int step = 1,
              QWidget * parent = 0, const char * name = 0 ): 
    QSpinBox ( minValue, maxValue, step , parent , name ) 
  {} 
  
  /*! distructor */ 
  ~MySpinBox () {} 
  
  /*! mapValueToText - map value to text  
   *\ param value - the real spin box value 
   *\ return a text of a double number we want to represent 
   */ 
  QString mapValueToText( int value ) 
  { 
    return QString( "%1.%2" ) // 0.0 to 10.0 
      .arg( value / 10 ).arg( value % 10 ); 
  } 
  
  /*! mapTextToValue - map text to value  
   *\ return the corresponding int value 
   */ 
  int mapTextToValue( bool *ok ) 
  { 
    return (int) ( 10 * text().toFloat() ); // 0 to 100 
  } 
}; 

/*! class PropertiesForm is the dialog form that allow the user
 *  to set the program properties.
 */
class Qt_widget_base_tab;
class PropertiesForm : public QDialog
{
  Q_OBJECT
public:
  PropertiesForm( QTabWidget * bar = 0 , QWidget* parent = 0 ,
                  int number_of_tabs = 0 , Qt_widget_base_tab *w_demo_p = 0, 
				  double scale = 0 , bool colors_flag = true);
  ~PropertiesForm() {}
  
  QLabel *textLabel1;
  QLabel *textLabel2;
  QLabel *textLabel3;
  QLabel *textLabel4;
  QLabel *textLabel5;
  QLabel *textLabel6;
  QLabel *textLabel7;
  QLabel *textLabel8;
  QLabel *textLabel9;
  QSpinBox *box1;
  QSpinBox *box2;
  QSpinBox *box3;
  MySpinBox *box4;
  QComboBox *box5;
  QSpinBox *box6;
  QComboBox *box7;
  QSpinBox *box8;
  QComboBox *box9;
  QPushButton *okPushButton;
  QPushButton *cancelPushButton;
  
protected:
  QVBoxLayout *optionsFormLayout;
  QHBoxLayout *arrLayout1;
  QHBoxLayout *arrLayout2;
  QHBoxLayout *arrLayout3;
  QHBoxLayout *arrLayout4;
  QHBoxLayout *arrLayout5;
  QHBoxLayout *arrLayout6;
  QHBoxLayout *arrLayout7;
  QHBoxLayout *arrLayout8;
  QHBoxLayout *arrLayout9;
  QHBoxLayout *buttonsLayout;
  
private:
  QTabWidget *myBar;
};

/*! OverlayForm - the dialog form of the overlay operation */
class OverlayForm : public QDialog
{
  Q_OBJECT
public:
  OverlayForm( QTabWidget * bar = 0 , QWidget* parent = 0 ,
               int number_of_tabs = 0 , const char* name = "options form",
               bool modal = FALSE, WFlags f = 0  );
  ~OverlayForm() {}
  
  QLabel *textLabel1;
  QLabel *textLabel2;
  QSplitter *split;
  DDListBox *listBox1;
  DDListBox *listBox2;
  QPushButton *okPushButton;
  QPushButton *cancelPushButton;
  
protected:
  QVBoxLayout *optionsFormLayout;
  QHBoxLayout *arrLayout;
  QHBoxLayout *buttonsLayout;
  
private:
  QTabWidget *myBar;
}; 

/*! class DDListBox used for the drag/drop action in the overlay form */ 
class DDListBox : public QListBox 
{ 
  Q_OBJECT 
public: 
  DDListBox( QWidget * parent = 0, const char * name = 0, WFlags f = 0 ); 
  // Low-level drag and drop 
  void dragEnterEvent( QDragEnterEvent *evt ); 
  void dropEvent( QDropEvent *evt ); 
  void mousePressEvent( QMouseEvent *evt ); 
  void mouseMoveEvent( QMouseEvent * ); 
  void set_max_items(int num); 
private: 
  int dragging; 
  unsigned int max_items; 
  bool flag; 
  
}; 

///*! class CheckItem used for the drag/drop action in the overlay form */ 
//class CheckItem : public QListBoxPixmap  
//{   
//public: 
//  CheckItem(  QListBox * listbox, const QPixmap & pix, const QString & text );  
//private: 
//  QCheckBox *check_box;
//}; 

/*! class OptionsForm used for choosing which conic type will be inserted */ 
class OptionsForm : public QDialog
{
  Q_OBJECT
public:
  OptionsForm( QWidget* parent = 0 ,int number_of_tabs = 0 ,const char*
               name = "options form", bool modal = FALSE, WFlags f = 0);
  ~OptionsForm() {}
  
  QLabel *textLabel1;
  QComboBox *arrComboBox1;
  QPushButton *okPushButton;
  QPushButton *cancelPushButton;
  
protected:
  QVBoxLayout *optionsFormLayout;
  QHBoxLayout *arrLayout1;
  QHBoxLayout *buttonsLayout;
  
}; 

/*! class CheckForm used for choosing which conic type will be inserted */ 
class CheckForm : public QDialog
{
  Q_OBJECT
public:
  CheckForm( OverlayForm *overlay_form , QWidget* parent = 0);
  ~CheckForm() {}
  
  QVButtonGroup *button_group;
  QPushButton *okPushButton;
  QPushButton *cancelPushButton;
  
protected:
  QVBoxLayout *optionsFormLayout;
  QHBoxLayout *layout;
  QHBoxLayout *buttonsLayout;
}; 

/*! class FileOpenOptionsForm used for choosing which action will be taken  
  when we open a new file 
*/ 
class FileOpenOptionsForm : public QDialog
{
  Q_OBJECT
public:
  FileOpenOptionsForm( bool flag = true ,QWidget* parent = 0 ,
                       const char* name = "file open options form", 
                       bool modal = FALSE, WFlags f = 0);
  ~FileOpenOptionsForm() {}
  
  QLabel *textLabel1;
  //QComboBox *arrComboBox1;
  QRadioButton *b1;
  QRadioButton *b2;
  QRadioButton *b3;
  QButtonGroup *buttonGroup;
  QPushButton *okPushButton;
  QPushButton *cancelPushButton;
  
protected:
  QVBoxLayout *optionsFormLayout;
  //QHBoxLayout *arrLayout1;
  QHBoxLayout *buttonsLayout;
  
}; 

/*! class RayShootingOptionsForm used for choosing the rayshoot diraction */ 
class RayShootingOptionsForm : public QDialog
{
  Q_OBJECT
public:
  RayShootingOptionsForm( QWidget* parent = 0 ,int number_of_tabs = 0 ,
	  const char* name = "options form", bool modal = FALSE, WFlags f = 0);
  ~RayShootingOptionsForm() {}
  
  QLabel *textLabel1;
  QComboBox *arrComboBox1;
  QPushButton *okPushButton;
  QPushButton *cancelPushButton;
  
protected:
  QVBoxLayout *optionsFormLayout;
  QHBoxLayout *arrLayout1;
  QHBoxLayout *buttonsLayout;
  
}; 


/*! class PointLocationStrategyForm used for choosing strategy for point location*/
class PointLocationStrategyForm : public QDialog
{
  Q_OBJECT
public:
  PointLocationStrategyForm(QWidget* parent = 0 ,int number_of_tabs = 0 ,
	  const char* name = "options form", bool modal = FALSE, WFlags f = 0);
  ~PointLocationStrategyForm() {}

  QLabel *textLabel1;
  QComboBox *arrComboBox1;
  QPushButton *okPushButton;
  QPushButton *cancelPushButton;
  
protected:
  QVBoxLayout *optionsFormLayout;
  QHBoxLayout *arrLayout1;
  QHBoxLayout *buttonsLayout;
  
}; 






#endif // FORMS_H

