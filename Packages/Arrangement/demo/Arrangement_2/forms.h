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

#include "cgal_types1.h"

class QDragEnterEvent;
class QDragDropEvent;
class DDListBox;
class QWidget;

////////////////////////////////////////////////////////////////////////

class OptionsForm : public QDialog
{
    Q_OBJECT
public:
    OptionsForm( QTabWidget * bar = 0 , QWidget* parent = 0 ,int number_of_tabs = 0 , const char* name = "options form",
		 bool modal = FALSE, WFlags f = 0  );
    ~OptionsForm() {}
  
    QLabel *textLabel1;
	QLabel *textLabel2;
    QComboBox *arrComboBox1;
	QComboBox *arrComboBox2;
    QPushButton *okPushButton;
    QPushButton *cancelPushButton;

protected:
    QVBoxLayout *optionsFormLayout;
    QHBoxLayout *arrLayout1;
	QHBoxLayout *arrLayout2;
    QHBoxLayout *buttonsLayout;

private:
	QTabWidget *myBar;
};
///////////////////////////////////////////////////////////////////////////////////
class MySpinBox : public QSpinBox
{
    Q_OBJECT
public:
	MySpinBox ( int minValue, int maxValue, int step = 1, QWidget * parent = 0, const char * name = 0 ):
		QSpinBox ( minValue, maxValue, step , parent , name )
	{}

	~MySpinBox () {}

    QString mapValueToText( int value )
    {
        return QString( "%1.%2" ) // 0.0 to 10.0
            .arg( value / 10 ).arg( value % 10 );
    }

    int mapTextToValue( bool *ok )
    {
        return (int) ( 10 * text().toFloat() ); // 0 to 100
    }
};

///////////////////////////////////////////////////////////////////////////////////
class PropertiesForm : public QDialog
{
    Q_OBJECT
public:
    PropertiesForm( QTabWidget * bar = 0 , QWidget* parent = 0 ,int number_of_tabs = 0 , const char* name = "options form",
		 bool modal = FALSE, WFlags f = 0  );
    ~PropertiesForm() {}
  
    QLabel *textLabel1;
	QLabel *textLabel2;
	QLabel *textLabel3;
	QLabel *textLabel4;
	QSpinBox *box1;
	QSpinBox *box2;
	QSpinBox *box3;
	MySpinBox *box4;
    QPushButton *okPushButton;
    QPushButton *cancelPushButton;

protected:
    QVBoxLayout *optionsFormLayout;
    QHBoxLayout *arrLayout1;
	QHBoxLayout *arrLayout2;
	QHBoxLayout *arrLayout3;
	QHBoxLayout *arrLayout4;
    QHBoxLayout *buttonsLayout;

private:
	QTabWidget *myBar;
};
////////////////////////////////////////////////////////////////////////

class OverLay : public QDialog
{
    Q_OBJECT
public:
    OverLay( QTabWidget * bar = 0 , QWidget* parent = 0 ,int number_of_tabs = 0 , const char* name = "OverLay",
		 bool modal = FALSE, WFlags f = 0  );
    ~OverLay() {}
  
    QLabel *textLabel1;
	QLabel *textLabel2;
	QLabel *textLabel3;
    QComboBox *arrComboBox1;
	QComboBox *arrComboBox2;
	QComboBox *arrComboBox3;
    QPushButton *okPushButton;
    QPushButton *cancelPushButton;

protected:
    QVBoxLayout *optionsFormLayout;
    QHBoxLayout *arrLayout1;
	QHBoxLayout *arrLayout2;
	QHBoxLayout *arrLayout3;
    QHBoxLayout *buttonsLayout;

private:
	QTabWidget *myBar;
};
////////////////////////////////////////////////////////////////////////

class OverlayForm : public QDialog
{
    Q_OBJECT
public:
    OverlayForm( QTabWidget * bar = 0 , QWidget* parent = 0 ,int number_of_tabs = 0 , const char* name = "options form",
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
//////////////////////////////////////////////////////////////////////

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

///////////////////////////////////////////////////////////////////////



