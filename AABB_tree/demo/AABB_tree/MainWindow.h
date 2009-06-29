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

			// view options
			void on_actionView_polyhedron_triggered();
			void on_actionView_points_triggered();
			void on_actionView_segments_triggered();

			// algorithms
			void on_actionInside_points_triggered();
			void on_actionBoundary_points_triggered();
			void on_actionBoundary_segments_triggered();
			void on_actionEdge_points_triggered();

			// distance functions
			void on_actionUnsigned_distance_function_to_facets_triggered();

			// benchmarks
			void on_actionBench_distances_triggered();
			void on_actionBench_intersections_triggered();

			// edit
			void on_actionClear_points_triggered();
			void on_actionClear_segments_triggered();

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
