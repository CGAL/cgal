#include <CGAL/KDS/IO/internal/KDS_Qt_examiner_viewer.h>
#include <qpushbutton.h>
#include <Inventor/nodes/SoSeparator.h>
#include <CGAL/KDS/IO/internal/KDS_pixmaps.h>

#include "KDS_Qt_examiner_viewer.moc"

namespace CGAL {
  namespace KDS {
    namespace internal {
      Qt_examiner_viewer::Qt_examiner_viewer(QWidget * parent):
	SoQtExaminerViewer(parent, NULL, TRUE,
			   SoQtFullViewer::BUILD_ALL, 
			   SoQtFullViewer::BROWSER,
			   FALSE){
	// Explicitly trigger the construction of viewer decorations.
	QWidget * widget = this->buildWidget(this->getParentWidget());
	this->setBaseWidget(widget);
	root_= new SoSeparator;
	this->setSceneGraph(root_);
      }

      void Qt_examiner_viewer::new_subgraph(SoNode *p){
	  root_->addChild(p);
	}

      void Qt_examiner_viewer::delete_subgraph(SoNode *p){
	  root_->removeChild(p);
	}

	

#define SETUP_QT_BUTTON(name) 	name ## _button_ = new QPushButton(parent); \
	    name##_button_->setFocusPolicy(QWidget::NoFocus);		\
	    name##_button_->setPixmap(QPixmap((const char **) name##_xpm)); \
	    name##_button_->adjustSize();				\
	    QObject::connect(name##_button_, SIGNAL(clicked()),		\
			     &core_, SLOT(name##_button()));		\
	    buttonlist->append(name##_button_);
      


      void Qt_examiner_viewer::createViewerButtons(QWidget * parent, SbPList * buttonlist){
	SoQtExaminerViewer::createViewerButtons(parent, buttonlist);
	// [now add your own button(s) to the buttonlist]
	SETUP_QT_BUTTON(play);
	SETUP_QT_BUTTON(pause);
	SETUP_QT_BUTTON(stop);
	SETUP_QT_BUTTON(play_to);
	SETUP_QT_BUTTON(play_through);
	SETUP_QT_BUTTON(reverse);
	SETUP_QT_BUTTON(faster);
	SETUP_QT_BUTTON(slower);
      }

	
      




    };
  }
}
