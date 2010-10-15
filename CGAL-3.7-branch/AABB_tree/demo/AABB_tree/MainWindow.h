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

			// drag & drop
			void dropEvent(QDropEvent *event);
			void closeEvent(QCloseEvent *event);
			void dragEnterEvent(QDragEnterEvent *event);

			// file menu
			void on_actionLoadPolyhedron_triggered();

			// edit menu
      void on_actionSave_snapshot_triggered();
      void on_actionCopy_snapshot_triggered();
			void on_actionClear_points_triggered();
			void on_actionClear_segments_triggered();
      void on_actionClear_cutting_plane_triggered();

			// algorithm menu
      void on_actionRefine_loop_triggered();
			void on_actionEdge_points_triggered();
			void on_actionInside_points_triggered();
			void on_actionBoundary_points_triggered();
			void on_actionRefine_bisection_triggered();
			void on_actionBoundary_segments_triggered();
      void on_actionPoints_in_interval_triggered();
			void on_actionSigned_distance_function_to_facets_triggered();
			void on_actionUnsigned_distance_function_to_edges_triggered();
			void on_actionUnsigned_distance_function_to_facets_triggered();
      void on_actionIntersection_cutting_plane_triggered();
    void on_actionCutting_plane_none_triggered();

			// benchmark menu
			void on_actionBench_memory_triggered();
			void on_actionBench_distances_triggered();
			void on_actionBench_intersections_triggered();

			// benchmark against #triangles
			void on_actionBench_construction_triggered();
			void on_actionBench_distances_vs_nbt_triggered();
			void on_actionBench_intersections_vs_nbt_triggered();

			// view menu
			void on_actionView_points_triggered();
			void on_actionView_segments_triggered();
			void on_actionView_polyhedron_triggered();
      void on_actionView_cutting_plane_triggered();

private:
	Scene* m_pScene;
	Viewer* m_pViewer;
	Ui::MainWindow* ui;
};

#endif // ifndef MAINWINDOW_H
