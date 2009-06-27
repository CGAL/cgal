#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QtOpenGL/qgl.h>
#include <CGAL/Qt/DemosMainWindow.h>

class QDragEnterEvent;
class QDropEvent;
class Scene;
class Viewer;
namespace Ui {
	class MainWindow;
}

#include "types.h"

class MainWindow : 
	public CGAL::Qt::DemosMainWindow
{
	Q_OBJECT
public:
	MainWindow(QWidget* parent = 0);
	~MainWindow();

	public slots:
		void updateViewerBBox();
		void open(QString filename);

		void setAddKeyFrameKeyboardModifiers(Qt::KeyboardModifiers);

		protected slots:

			// settings
			void quit();
			void readSettings();
			void writeSettings();

			// load
			void on_actionLoadPolyhedron_triggered();

			// algorithms
			void on_actionInside_points_triggered();

			// drag & drop
			void dragEnterEvent(QDragEnterEvent *event);
			void dropEvent(QDropEvent *event);
			void closeEvent(QCloseEvent *event);

private:
	QString strippedName(const QString &fullFileName);

	Scene* m_pScene;
	Viewer* m_pViewer;
	Ui::MainWindow* ui;
};

#endif // ifndef MAINWINDOW_H
