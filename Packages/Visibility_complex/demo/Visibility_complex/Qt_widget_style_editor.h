#ifndef QT_WIDGET_STYLE_EDITOR_H
#define QT_WIDGET_STYLE_EDITOR_H

#include <qcolor.h>
#include <qcolordialog.h>
#include <qpushbutton.h>
#include <qspinbox.>
#include <qcombobox.h>
#include <qscrollview.h>
#include <qlabel.h>
#include <qvbox.h
#include <qhbox.h>
#include <qgridlayout.h>

namespace CGAL {

namespace Qt_widget {
  class Color_selector : public QPushButton
  {
    Q_OBJECT
  public:
    Color_label(QColor c = Qt::black,
		QWidget* parent = 0, const char* name = 0)
      : QPushButton(parent, name)
    {
      setColor(c);
      connect(this, SIGNAL(clicked()),
	      this, SLOT(color_dialog()) );
    }

    void setColor(QColor c)
    {
      color = c;

      QPixmap pix(24,20);
      pix->fill(c);
      setPixmap(pix);
    }

    QColor value() const
    {
      return color;
    }

  private slots:
    void color_dialog()
    {
      setColor(QColorDialog::getColor());
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
  };

  class Bool_selector : public QComboBox
  {
    Q_OBJECT
  public:
    Bool_selector(bool b_, QWidget *parent = 0, const char *name = 0)
      : QComboBox(false, parent, name), b(b_)
    {
      insertItem("False");
      insertItem("True");

      if(b)
	setCurrentItem(1);
      else
	setCurrentItem(0);
    }

    bool value() const
    {
      return currentItem() == 1;
    }
  };

  class Point_style_selector : public QComboBox
  {
    Q_OBJECT
  public:
    typdef ::CGAL::PointStyle PointStyle;
    Point_style_selector(PointStyle s,
			 QWidget *parent = 0, const char *name = 0)
      : QComboBox(false, parent, name), b(b_)
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

    PointStyle value() const
    {
      return PointStyle(currentItem());
    }
  };

} // end namespace Qt_widget

class Qt_widget_style_editor : public QScrollView {
  Q_OBJECT
public:

 typedef QMap<QString,QVariant> Style;

  Qt_widget_style_editor(Style* style,
			 QWidget *parent, const char *name )
    : QScrollView(parent, name)
  {
    using namespace ::CGAL::Qt_widget;
    typedef Style::const_iterator iterator;

    QGridLayout layout = new QGridLayout(style->size(), 3, viewport());
    layout->addColSpacing(1, 5);

    for(iterator it = style.begin();
	it != style.end; ++it)
      {
	QLabel label = new QLabel(it->first, hbox);
	label->setText(it->first);
	layout->addItem(label);

	switch(it->second.type()) {
	case Color:
	  Color_selector color_sel = new Color_selector(it->second.toColor(),
							viewport());
	  break;
	case Int:
	  Int_selector int_sel = new Int_selector(it->second.toInt(),
						  viewport());
	  break;
	case Bool:
	  Bool_selector bool_sel = new Bool_selector(it->second.toBool(),
						     viewport());
	  break;
	case Uint:
	  Point_style_selector ps_sel = 
	    new Point_style_selector(PointStyle(it->second.toUint()),
				     viewport());
	default:
	  CGAL_assertion(false);
	  break;
	}
      }
  }
} // end of class Qt_widget_style_editor

} // end namespace CGAL

#endif // QT_WIDGET_STYLE_EDITOR_H
