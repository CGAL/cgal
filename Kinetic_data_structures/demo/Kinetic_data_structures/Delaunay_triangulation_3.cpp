#include <CGAL/Kinetic/Insert_event.h>
#include <CGAL/Kinetic/Exact_simulation_traits.h>
#include <algorithm>
#include <CGAL/Kinetic/Delaunay_triangulation_3.h>

#ifdef CGAL_USE_COIN
#include "include/SoQt_widget_3.h"
#include "include/SoQt_moving_points_3.h"
#include "include/SoQt_triangulation_3.h"
#include <CGAL/Kinetic/Enclosing_box_3.h>
#include <Inventor/nodes/SoIndexedLineSet.h>
#include <Inventor/nodes/SoIndexedFaceSet.h>
#endif

#ifdef CGAL_USE_BOOST_PROGRAM_OPTIONS
#include <boost/program_options.hpp>
#endif

int main(int argc, char *argv[])
{
#ifdef CGAL_USE_COIN
    int n=10;
    int d=2;
    bool print_help=false;
    std::string file;
    bool verbose=false;
#ifdef CGAL_USE_BOOST_PROGRAM_OPTIONS
    boost::program_options::options_description desc("Allowed options");
    desc.add_options()
        ("help", boost::program_options::bool_switch(&print_help), "produce help message")
        ("verbose,v", boost::program_options::bool_switch(&verbose), "produce lots of output")
        ("num-points,n", boost::program_options::value<int>(&n), "Number of points to use.")
        ("degree,d", boost::program_options::value<int>(&d), "The degree of the motions to use.")
        ("file,f", boost::program_options::value<std::string>(&file), "Read points from a file.");

    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::command_line_parser(argc, argv).
        options(desc).run(), vm);
    boost::program_options::notify(vm);

    if (print_help) {
        std::cout << desc << "\n";
        return EXIT_FAILURE;
    }
#endif

    typedef CGAL::Kinetic::Exact_simulation_traits Traits;
    typedef CGAL::Kinetic::Delaunay_triangulation_3<Traits> KDel;
    typedef Traits::Simulator::Time Time;
    typedef CGAL::Kinetic::Insert_event<Traits::Active_points_3_table> IE;
    typedef Traits::Active_points_3_table::Data MP;
    typedef Traits::Kinetic_kernel::Motion_function MF;
    typedef Traits::Kinetic_kernel::Point_3 MP;
    typedef CGAL::Kinetic::SoQt_widget_3<Traits::Simulator> Qt_gui;
    typedef CGAL::Kinetic::SoQt_moving_points_3<Traits, Qt_gui> Qt_mpt;
    typedef CGAL::Kinetic::SoQt_triangulation_3<KDel, Qt_gui, Qt_mpt> Qt_del;

    Traits tr(0,100000);
    Qt_gui::Handle qtsim= new Qt_gui(argc, argv, tr.simulator_handle());
    Qt_mpt::Handle qtmpt= new Qt_mpt(tr, qtsim);
    KDel::Handle kdel= new KDel(tr);
    Qt_del::Handle cd= new Qt_del(kdel, qtsim, qtmpt);

    typedef Traits::Simulator::Time Time;
    typedef CGAL::Kinetic::Insert_event<Traits::Active_points_3_table> MOI;

    Traits::Kinetic_kernel::Function_kernel::Construct_function cf= tr.kinetic_kernel_object().function_kernel_object().construct_function_object();
//sim->end_time();
    Traits::Simulator::Handle sim= tr.simulator_handle();
    Traits::Active_points_3_table::Handle mpt= tr.active_points_3_table_handle();
    CGAL::Kinetic::Enclosing_box_3<Traits> eb(tr,-10,10,-10,10,-10,10);

    {
      CGAL::SoQt_handle<SoShapeKit> shape= new SoShapeKit;
      shape->setName("Edges");
      CGAL::SoQt_handle<SoShapeKit> walls= new SoShapeKit;
      walls->setName("Walls");
      {
	CGAL::SoQt_handle<SoAppearanceKit> ap= new SoAppearanceKit;
	{
	  CGAL::SoQt_handle<SoMaterial> mat= new SoMaterial;
	  mat->setName("Facet_material");
	  mat->ambientColor.setValue(SbColor(0.6,0.6,0.6));
	  mat->specularColor.setValue(SbColor(0.0, 0.0, 0.0));
	  mat->emissiveColor.setValue(SbColor(0.0,0.0,0.0));
	  mat->shininess.setValue(.2);
	  mat->transparency.setValue(0.5);
	  mat->diffuseColor.setValue(SbColor(0.6, 0.6, 0.6));
	  ap->setPart("material", mat.get());
	}
	shape->setPart("appearance", ap.get());
	walls->setPart("appearance", ap.get());
      }
      CGAL::SoQt_handle<SoCoordinate3> coords= new SoCoordinate3();
      SbVec3f cids[8]={SbVec3f(-10,-10,-10),
		       SbVec3f(-10,-10,10),
		       SbVec3f(10, -10, -10),
		       SbVec3f(-10, 10, -10),
		       SbVec3f(10, 10, 10),
		       SbVec3f(-10, 10, 10),
		       SbVec3f(10, -10, 10),
		       SbVec3f(10,10,-10)};
      coords->point.setValues(0,8,cids);


      shape->setPart("coordinate3", coords.get());
      walls->setPart("coordinate3", coords.get());
      {
	CGAL::SoQt_handle<SoShapeHints> hint= new SoShapeHints;
	hint->vertexOrdering.setValue(SoShapeHints::COUNTERCLOCKWISE);
	hint->shapeType.setValue(SoShapeHints::UNKNOWN_SHAPE_TYPE);
	hint->faceType.setValue(SoShapeHints::CONVEX);
	hint->creaseAngle.setValue(0);
	walls->setPart("shapeHints", hint.get());
      }
      {
	CGAL::SoQt_handle<SoIndexedLineSet> ls= new SoIndexedLineSet();
	int lines[]={0,2,7,3,5,4,6,1,0,-1,
		     0,3,-1,
		     1,5,-1,
		       2,6,-1,
		       4,7, -1,
	};
	//CGAL_assertion(sizeof(lines)/sizeof(int)==19);
	ls->coordIndex.setNum(sizeof(lines)/sizeof(int));
	ls->coordIndex.setValues(0, sizeof(lines)/sizeof(int), lines);
	//ls->vertexProperty.setValue(coords.get());

	shape->setPart("shape",ls.get());
      }
      {
	CGAL::SoQt_handle<SoIndexedFaceSet> ls= new SoIndexedFaceSet();
	int faces[]={2,0,3,7,-1,
		     0,1,5,3,-1,
		     0,2,6,1,-1,
		     2,7,4,6,-1
	};
	//CGAL_assertion(sizeof(lines)/sizeof(int)==19);
	ls->coordIndex.setNum(sizeof(faces)/sizeof(int));
	ls->coordIndex.setValues(0, sizeof(faces)/sizeof(int), faces);
	//ls->vertexProperty.setValue(coords.get());

	walls->setPart("shape",ls.get());
      }
      CGAL::SoQt_handle<SoSeparator> p= new SoSeparator();
      p->addChild(shape.get());
      p->addChild(walls.get());
      qtsim->soqt_examiner_viewer_pointer()->new_subgraph(p.get());
    }


    if (verbose) {
        CGAL_KINETIC_SET_LOG_LEVEL(CGAL::Kinetic::LOG_LOTS);
    }
    else {
        CGAL_KINETIC_SET_LOG_LEVEL(CGAL::Kinetic::LOG_NONE);
    }

    if (file.empty()) {
        CGAL::Random rand;
        for (int i=0; i< n; ++i) {
            std::vector<double> coefsx, coefsy, coefsz;
            for (int j=0; j< d; ++j) {
                coefsx.push_back((rand.get_double()*10-5)/(j+1));
                coefsy.push_back((rand.get_double()*10-5)/(j+1));
                coefsz.push_back((rand.get_double()*10-5)/(j+1));
            }
            Traits::Kinetic_kernel::Point_3 mp(Traits::Kinetic_kernel::Motion_function(coefsx.begin(), coefsx.end()),
					       Traits::Kinetic_kernel::Motion_function(coefsy.begin(), coefsy.end()),
					       Traits::Kinetic_kernel::Motion_function(coefsz.begin(), coefsz.end()));
            tr.active_points_3_table_handle()->insert(mp);
        }
    }
    else {
        std::ifstream in(file.c_str());
        if (!in) {
            std::cerr << "Error opening input file: " << file << std::endl;
            return EXIT_FAILURE;
        }
        char buf[1000];
        int nread=0;
        while (true ) {
            in.getline(buf, 1000);
            if (!in) break;
            std::istringstream il(buf);
            Traits::Kinetic_kernel::Point_3 p;
            il >> p;
            tr.active_points_3_table_handle()->insert(p);
            ++nread;
        }
        std::cout << nread << " points read.\n";
    }

    kdel->set_has_certificates(true);

    std::cout << "This program displays a 3D kinetic Delaunay triangulation.\n";
    std::cout << "What is displayed can be controlled by pushing the followin keys over the Coin window and with the arrow tool selected.\n";
    std::cout << "Press 'h' to hide/show the convex hull.\n";
    std::cout << "Press 'f' to hide/show faces.\n";
    std::cout << "Press 's' to show spheres for points, 'p' to show points.\n";


    return qtsim->begin_event_loop();
#else 
    std::cerr << "Coin and SoQt are required for this demo.\n";
    return EXIT_FAILURE;
#endif
};
