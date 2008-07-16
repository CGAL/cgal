#include "MainWindow.h"
#include "Scene.h"

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Nef_polyhedron_3.h> 

// quick hacks to convert polyhedra from exact to inexact and vice-versa

void from_exact(Exact_polyhedron& in,
		Polyhedron& out)
{
	std::ofstream out_stream("tmp.off");
	out_stream << in;
	std::ifstream in_stream("tmp.off");
	in_stream >> out;
}

void to_exact(Polyhedron& in,
	      Exact_polyhedron& out)
{
	std::ofstream out_stream("tmp.off");
	out_stream << in;
	std::ifstream in_stream("tmp.off");
	in_stream >> out;
}

void MainWindow::on_actionUnion_triggered()
{
	boolean_operation(BOOLEAN_UNION);
}

void MainWindow::on_actionIntersection_triggered()
{
	boolean_operation(BOOLEAN_INTERSECTION);
}

void MainWindow::on_actionDifference_triggered()
{
	boolean_operation(BOOLEAN_DIFFERENCE);
}

void MainWindow::boolean_operation(const int operation)
{
	QApplication::setOverrideCursor(Qt::WaitCursor);

	// PA: to be done by LR?
	//int indexA = scene->getPolygonAIndex();
	//int indexB = scene->getPolygonBIndex();

	// to remove
	int indexA = 0;
	int indexB = 1;

	Polyhedron* polyA = scene->polyhedron(indexA);
	Polyhedron* polyB = scene->polyhedron(indexB);
	if(!polyA) return;
	if(!polyB) return;

	typedef CGAL::Nef_polyhedron_3<Exact_Kernel> Nef_polyhedron; 

	Exact_polyhedron exact_polyA; 
	to_exact(*polyA,exact_polyA);

	Exact_polyhedron exact_polyB; 
	to_exact(*polyB,exact_polyB);

	// convert to nef polyhedra
	QTime time;
	time.start();
	std::cout << "Convert to nef polyhedra...";
	Nef_polyhedron n1(exact_polyA); 
	Nef_polyhedron n2(exact_polyB); 
	std::cout << "ok (" << time.elapsed() << " ms)" << std::endl;

	// perform Boolean operation
	std::cout << "Boolean operation...";
	time.start();
	switch(operation)
	{
		case BOOLEAN_UNION:
			n1 += n2;
			break;
		case BOOLEAN_INTERSECTION:
			n1 *= n2;
			break;
		case BOOLEAN_DIFFERENCE:
			n1 -= n2;
	}
	std::cout << "ok (" << time.elapsed() << " ms)" << std::endl;

	// save the exact resulting mesh
	Exact_polyhedron exact_result;
	n1.convert_to_Polyhedron(exact_result);

	// reload as inexact one
	Polyhedron *pResult = new Polyhedron;
	from_exact(exact_result,*pResult);

	switch(operation)
	{
		case BOOLEAN_UNION:
			scene->addPolyhedron(pResult,
				tr("%1 union %2").arg(scene->polyhedronName(indexA),scene->polyhedronName(indexB)),
				Qt::yellow,
				true,
				scene->polyhedronRenderingMode(indexA));
			break;
		case BOOLEAN_INTERSECTION:
			scene->addPolyhedron(pResult,
				tr("%1 intersection %2").arg(scene->polyhedronName(indexA),scene->polyhedronName(indexB)),
				Qt::yellow,
				true,
				scene->polyhedronRenderingMode(indexA));
			break;
		case BOOLEAN_DIFFERENCE:
			scene->addPolyhedron(pResult,
				tr("%1 minus %2").arg(scene->polyhedronName(indexA),scene->polyhedronName(indexB)),
				Qt::yellow,
				true,
				scene->polyhedronRenderingMode(indexA));
	}

	QApplication::setOverrideCursor(Qt::ArrowCursor);
}

