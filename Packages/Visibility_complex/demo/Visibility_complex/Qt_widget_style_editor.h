#ifndef QT_WIDGET_STYLE_EDITOR_H
#define QT_WIDGET_STYLE_EDITOR_H

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

namespace CGAL {

namespace Qt_widget_internals {
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

} // end namespace Qt_widget

class Qt_widget_style_editor : public QScrollView {
  Q_OBJECT
public:

  typedef Qt_widget_styled_layer::Style Style;

  Qt_widget_style_editor(Style* style,
			 QWidget *parent = 0 , const char *name = 0)
    : QScrollView(parent, name)
  {
    using namespace ::CGAL::Qt_widget_internals;
    typedef Style::const_iterator iterator;

    setSizePolicy(QSizePolicy(QSizePolicy::Maximum, QSizePolicy::Maximum);

    QGrid* layout = new QGrid(2, viewport());
    layout->setSpacing(5);

    int i = 0;
    for(iterator it = style->begin();
	it != style->end(); ++it)
      {
	QLabel* label = new QLabel(it.key(), layout);
	label->setText(it.key());

	QWidget* selector;

	switch(it.data().type()) {
	case QVariant::Color:
	  selector = new Color_selector(it.data().toColor(),
					layout);
	  break;
	case QVariant::Int:
	  selector = new Int_selector(it.data().toInt(),
				      layout);
	  break;
	case QVariant::Bool:
	  selector = new Bool_selector(it.data().toBool(),
				       layout);
	  break;
	case QVariant::UInt:
	  selector = 
	    new Point_style_selector(PointStyle(it.data().toUInt()),
				     layout);
	  break;
	default:
	  CGAL_assertion(false);
	  break;
	};
	++i;
      }
  }

  virtual ~Qt_widget_style_editor() {}
}; // end of class Qt_widget_style_editor

} // end namespace CGAL

#endif // QT_WIDGET_STYLE_EDITOR_H
