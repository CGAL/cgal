#include "Qt_widget_style_editor.h"

#include <qcolor.h>
#include <qcolordialog.h>
#include <qpushbutton.h>
#include <qspinbox.h>
#include <qcombobox.h>
#include <qscrollview.h>
#include <qlabel.h>
#include <qvbox.h>
#include <qhbox.h>
#include <qlayout.h>
#include <qgrid.h>
#include <qpixmap.h>
#include <qvariant.h>
#include <qsizepolicy.h>
#include "Qt_widget_styled_layer.h"

class Color_selector : public QPushButton
{
  Q_OBJECT
public:
  Color_selector(QColor c = Qt::black,
		 QWidget* parent = 0, const char* name = 0)
    : QPushButton(parent, name)
  {
    setColor(c);
    connect(this, SIGNAL(clicked()),
	    this, SLOT(color_dialog()) );
  }

  virtual ~Color_selector() {};

  void setColor(QColor c)
  {
    color = c;

    QPixmap pix(24,20);
    pix.fill(c);
    setPixmap(pix);
  }

  QColor value() const
  {
    return color;
  }

private slots:
  void color_dialog()
  {
    QColor c = QColorDialog::getColor(value());
    if( c.isValid() )
      setColor(c);
  }
private:
  QColor color;
};

class Int_selector : public QSpinBox
{
  Q_OBJECT
public:
  Int_selector(int i, QWidget *parent = 0, const char *name = 0)
    : QSpinBox(-INT_MAX, INT_MAX, 1, parent, name)
  {
    setValue(i);
  }

  virtual ~Int_selector() {};
};

class Bool_selector : public QComboBox
{
  Q_OBJECT
public:
  Bool_selector(bool b_, QWidget *parent = 0, const char *name = 0)
    : QComboBox(false, parent, name)
  {
    insertItem("False");
    insertItem("True");

    if(b_)
      setCurrentItem(1);
    else
      setCurrentItem(0);
  }

  virtual ~Bool_selector() {};

  bool value() const
  {
    return currentItem() == 1;
  }
};

class Point_style_selector : public QComboBox
{
  Q_OBJECT
public:
  typedef ::CGAL::PointStyle PointStyle;
  Point_style_selector(PointStyle s,
		       QWidget *parent = 0, const char *name = 0)
    : QComboBox(false, parent, name)
  {
    insertItem("Pixel");
    insertItem("Cross");
    insertItem("Plus");
    insertItem("Circle");
    insertItem("Disc");
    insertItem("Rect");
    insertItem("Box");

    setCurrentItem(static_cast<int>(s));
  }

  virtual ~Point_style_selector() {};

  PointStyle value() const
  {
    return PointStyle(currentItem());
  }
};

namespace CGAL {

Qt_widget_style_editor::Qt_widget_style_editor(Style* style,
					       QWidget *parent,
					       const char *name)
  : QFrame(parent, name)
{
  typedef Style::const_iterator iterator;

  QGridLayout* layout = new QGridLayout(this);
  layout->addColSpacing(1,5);

  const int labels_col = 0; // column number of labels
  const int selectors_col = 2; // column number of selectors

  int row = 0;
  for(iterator it=style->begin();
      it != style->end();
      ++it)
    {
      QLabel* label = new QLabel( it.key(), this);
      layout->addWidget(label, row, labels_col);

      QWidget* selector;
      switch( it.data().type() ) {
      case QVariant::Color:
	selector = new Color_selector(it.data().toColor(), this);
	break;
      case QVariant::Int:
	selector = new Int_selector(it.data().toInt(), this);
	break;
      default:
	break;
      }
      
      layout->addWidget(selector, row, selectors_col);
      ++row;
    }
}	

} // end namespace CGAL

// moc_source_file: Qt_widget_style_editor.h
#include "Qt_widget_style_editor.moc"

// moc_source_file: Qt_widget_style_editor.C
#include "Qt_widget_style_editor.C.moc"
