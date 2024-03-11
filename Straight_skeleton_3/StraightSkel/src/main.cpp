/**
 * @file   main.cpp
 * @author Gernot Walzl
 * @date   2011-12-19
 *
 * @mainpage
 * StraightSkel is an implementation of the Straight Skeleton in 2- and
 * 3-dimensional space. It is used to animate the computation of offsets
 * of polygons and polyhedrons.
 * StraightSkel was written by Gernot Walzl
 * mainly in the years 2011, 2012, 2013.
 */

#include "debug.h"
#include "typedefs_thread.h"
#include "util/StringFuncs.h"

#include "data/2d/ptrs.h"
#include "data/2d/Polygon.h"
#include "data/3d/ptrs.h"
#include "data/3d/KernelFactory.h"
#include "data/3d/Polyhedron.h"

#include "db/2d/ptrs.h"
#include "db/2d/DAOFactory.h"
#include "db/2d/PolygonDAO.h"
#include "db/2d/StraightSkeletonDAO.h"
#include "db/2d/FLMAFile.h"
#include "db/3d/ptrs.h"
#include "db/3d/DAOFactory.h"
#include "db/3d/PolyhedronDAO.h"
#include "db/3d/StraightSkeletonDAO.h"
#include "db/3d/OBJFile.h"
#include "db/3d/FLMAFile.h"

#include "algo/ptrs.h"
#include "algo/Controller.h"
#include "algo/2d/ptrs.h"
#include "algo/2d/PolygonTransformation.h"
#include "algo/2d/SimpleStraightSkel.h"
#include "algo/2d/FastStraightSkel.h"
#include "algo/2d/SkelMeshGenerator.h"
#include "algo/3d/ptrs.h"
#include "algo/3d/PolyhedronBuilder.h"
#include "algo/3d/PolyhedronTransformation.h"
#include "algo/3d/SimpleStraightSkel.h"
#include "algo/3d/GraphChecker.h"

#include "ui/gl/ptrs.h"
#include "ui/gl/MainOpenGLWindow.h"

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <list>
#include <string>
#include <vector>


void printUsage(const char* argv0) {
    std::cout << "Usage: " << argv0 << " [2d|3d] [...] [GLUT_OPTIONS]" << std::endl;
    std::cout << std::endl;
    std::cout << "  2d options:" << std::endl;
    std::cout << "    PolygonID" << std::endl;
    std::cout << "    skel SkelID" << std::endl;
    std::cout << std::endl;
    std::cout << "  3d options:" << std::endl;
    std::cout << "    PolyhedronID" << std::endl;
    std::cout << "    load filename.obj" << std::endl;
    std::cout << "    import filename.obj" << std::endl;
    std::cout << "    skel SkelID" << std::endl;
    std::cout << std::endl;
    std::cout << "  general options:" << std::endl;
    std::cout << "    --no-window" << std::endl;
    std::cout << "    --save" << std::endl;
    std::cout << "    --save-offsets -1.0,-1.5" << std::endl;
    std::cout << "    --config StraightSkel.ini" << std::endl;
    std::cout << std::endl;
    std::cout << "Example: " << argv0 << " 3d load anything.obj" << std::endl;
}

void printCommand(int argc, const char* argv[]) {
    std::cout << "Command: ";
    for (int i = 0; i < argc; i++) {
        if (i > 0) {
            std::cout << " ";
        }
        std::cout << argv[i];
    }
    std::cout << std::endl;
}

bool isSet(const char* option, int argc, const char* argv[]) {
    bool result = false;
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], option) == 0) {
            result = true;
        }
    }
    return result;
}

const char* getOption(const char* option, int argc, const char* argv[]) {
    const char* result = 0;
    for (int i = 1; i < (argc-1); i++) {
        if (strcmp(argv[i], option) == 0) {
            result = argv[i+1];
        }
    }
    return result;
}

std::list<double> parseCSV(const char* csv) {
    std::list<double> values;
    std::vector<std::string> str_vals = util::StringFuncs::split(csv, ",", false);
    for (unsigned int cnt = 0; cnt < str_vals.size(); cnt++) {
        values.push_back(std::stod(str_vals[cnt]));
    }
    values.sort(std::greater<double>());
    return values;
}


int main(int argc, const char* argv[]) {
    printCommand(argc, argv);

    if (argc < 3) {
        printUsage(argv[0]);
        return EXIT_FAILURE;
    }

    // set number of dimensions
    int num_dims = 0;
    if (strcmp("2d", argv[1])==0 || strcmp("2", argv[1])==0) {
        num_dims = 2;
    } else if (strcmp("3d", argv[1])==0 || strcmp("3", argv[1])==0) {
        num_dims = 3;
    }
    if (num_dims < 2 || 3 < num_dims) {
        printUsage(argv[0]);
        return EXIT_FAILURE;
    }

    util::ConfigurationSPtr config = util::Configuration::getInstance();
    std::string str_conf_file;
    const char* chr_conf_file = getOption("--config", argc, argv);
    if (chr_conf_file) {
        str_conf_file = chr_conf_file;
    } else {
        str_conf_file = config->findDefaultFilename();
    }
    if (!config->load(str_conf_file)) {
        std::cout << "Error: Config file '" << str_conf_file << "' not found." << std::endl;
        if (chr_conf_file) {
            return EXIT_FAILURE;
        }
    }

    bool rand_move_points = false;
    bool rand_move_points_when_degenerated = false;
    double rand_move_points_range = 0.001;
    bool translate_and_scale_polyhedron = false;
    bool translate_and_scale_view = false;
    float translate[3];
    for (unsigned int i = 0; i < 3; i++) {
        translate[i] = 0.0f;
    }
    float scale = 1.0f;
    if (config->isLoaded()) {
        rand_move_points = config->getBool("main", "rand_move_points");
        rand_move_points_when_degenerated =
                config->getBool("main", "rand_move_points_when_degenerated");
        double value = config->getDouble("main", "rand_move_points_range");
        if (value != 0.0) {
            rand_move_points_range = value;
        }
        translate_and_scale_polyhedron =
                config->getBool("main", "translate_and_scale_polyhedron");
        translate_and_scale_view =
                config->getBool("main", "translate_and_scale_view");
    }

    // load input
    int id = atoi(argv[2]);
    data::_2d::PolygonSPtr polygon;
    data::_2d::skel::StraightSkeletonSPtr skel2d;
    data::_3d::PolyhedronSPtr polyhedron;
    data::_3d::skel::StraightSkeletonSPtr skel3d;
    if (num_dims == 2) {
        db::_2d::PolygonDAOSPtr polygon_dao =
                db::_2d::DAOFactory::getPolygonDAO();
        if (strcmp("load", argv[2]) == 0) {
            const char* filename = argv[3];
            if (util::StringFuncs::endsWith(filename, ".flma")) {
                polygon = db::_2d::FLMAFile::load(filename);
            }
            if (!polygon) {
                std::cout << "Error: Unable to open '" << filename << "'." << std::endl;
                return EXIT_FAILURE;
            }
        } else if (strcmp("skel", argv[2]) == 0) {
            int skelid = atoi(argv[3]);
            db::_2d::StraightSkeletonDAOSPtr skel_dao =
                    db::_2d::DAOFactory::getStraightSkeletonDAO();
            skel2d = skel_dao->find(skelid);
            if (!skel2d) {
                std::cout << "Error: StraightSkeleton with SkelID=" << skelid
                    << " not found." << std::endl;
                return EXIT_FAILURE;
            } else {
                DEBUG_VAR(skel2d->toString());
                id = skel_dao->findPolyID(skelid);
                polygon = polygon_dao->find(id);
                if (polygon) {
                    skel2d->setPolygon(polygon);
                }
            }
        } else {
            polygon = polygon_dao->find(id);
            if (!polygon) {
                std::cout << "Error: Polygon with PolyID=" << id
                    << " not found." << std::endl;
                return EXIT_FAILURE;
            }
        }
        if (rand_move_points_when_degenerated && !rand_move_points) {
            std::cout << "Checking for parallel lines." << std::endl;
            if (algo::_2d::PolygonTransformation::hasParallelLines(polygon)) {
                std::cout << "Warning: Polygon has parallel lines." << std::endl;
                rand_move_points = true;
            }
        }
        if (rand_move_points) {
            std::cout << "Points will be moved randomly." << std::endl;
            algo::_2d::PolygonTransformation::randMovePoints(polygon, rand_move_points_range);
        }
        if (translate_and_scale_view) {
            data::_2d::Point2SPtr p_box_min =
                    algo::_2d::PolygonTransformation::boundingBoxMin(polygon);
            data::_2d::Point2SPtr p_box_max =
                    algo::_2d::PolygonTransformation::boundingBoxMax(polygon);
            float scale_min = std::numeric_limits<float>::max();
            for (unsigned int i = 0; i < 2; i++) {
                translate[i] = -((*p_box_max)[i] + (*p_box_min)[i])/2.0;
                float scale_cur = 20.0/((*p_box_max)[i] - (*p_box_min)[i]);
                if (scale_cur < scale_min) {
                    scale_min = scale_cur;
                }
                scale = scale_min;
            }
        }
        if (!polygon->isConsistent()) {
            std::cout << "Warning: Polygon with PolyID=" << id
                << " is not consistent." << std::endl;
        }
        DEBUG_VAR(polygon->toString());
    } else if (num_dims == 3) {
        db::_3d::PolyhedronDAOSPtr polyhedron_dao =
                db::_3d::DAOFactory::getPolyhedronDAO();
        if (strcmp("test", argv[2]) == 0) {
            data::_3d::Point3SPtr p1 =
                    data::_3d::KernelFactory::createPoint3(-5.0, -5.0, -5.0);
            data::_3d::Point3SPtr p2 =
                    data::_3d::KernelFactory::createPoint3(5.0, 5.0, -5.0);
            data::_3d::Point3SPtr p3 =
                    data::_3d::KernelFactory::createPoint3(5.0, -5.0, 5.0);
            data::_3d::Point3SPtr p4 =
                    data::_3d::KernelFactory::createPoint3(-5.0, 5.0, 5.0);
            polyhedron = algo::_3d::PolyhedronBuilder::makeTetrahedron(
                    p1, p2, p3, p4);
        } else if (strcmp("load", argv[2]) == 0) {
            const char* filename = argv[3];
            if (util::StringFuncs::endsWith(filename, ".obj")) {
                polyhedron = db::_3d::OBJFile::load(filename);
            } else if (util::StringFuncs::endsWith(filename, ".flma")) {
                polyhedron = db::_3d::FLMAFile::load(filename);
            }
            if (!polyhedron) {
                std::cout << "Error: Unable to open '" << filename << "'." << std::endl;
                return EXIT_FAILURE;
            }
        } else if (strcmp("skel", argv[2]) == 0) {
            int skelid = atoi(argv[3]);
            db::_3d::StraightSkeletonDAOSPtr skel_dao =
                    db::_3d::DAOFactory::getStraightSkeletonDAO();
            skel3d = skel_dao->find(skelid);
            if (!skel3d) {
                std::cout << "Error: StraightSkeleton with SkelID=" << skelid
                    << " not found." << std::endl;
                return EXIT_FAILURE;
            } else {
                DEBUG_VAR(skel3d->toString());
                id = skel_dao->findPolyhedronID(skelid);
                polyhedron = polyhedron_dao->find(id);
                if (polyhedron) {
                    skel3d->setPolyhedron(polyhedron);
                }
            }
        } else if (strcmp("import", argv[2]) == 0) {
            polyhedron = db::_3d::OBJFile::load(argv[3]);
            if (!polyhedron) {
                std::cout << "Error: Polyhedron '" << argv[3] << "' not found."
                    << std::endl;
                return EXIT_FAILURE;
            }
            if (polyhedron->isConsistent()) {
                if (polyhedron_dao->insert(polyhedron) > 0) {
                    std::cout << argv[3] << " imported: PolyhedronID="
                        << polyhedron->getID() << std::endl;
                    return EXIT_SUCCESS;
                } else {
                    std::cout << "Error: Not able to import polyhedron." << std::endl;
                    return EXIT_FAILURE;
                }
            } else {
                std::cout << "Error: Polyhedron is not consistent and "
                    << "will not be imported." << std::endl;
                return EXIT_FAILURE;
            }
        } else {
            polyhedron = polyhedron_dao->find(id);
            if (!polyhedron) {
                std::cout << "Error: Polyhedron with PolyhedronID=" << id
                    << " not found." << std::endl;
                return EXIT_FAILURE;
            }
        }
        if (rand_move_points_when_degenerated && !rand_move_points) {
            std::cout << "Checking if all combinations of 3 facet supporting planes intersect in a point." << std::endl;
            std::cout << "In case this takes too long, "
                << "you may disable 'rand_move_points_when_degenerated'." << std::endl;
            if (!algo::_3d::PolyhedronTransformation::doAll3PlanesIntersect(polyhedron)) {
                std::cout << "Warning: Not all combinations of 3 planes intersect." << std::endl;
                rand_move_points = true;
            }
        }
        if (translate_and_scale_polyhedron) {
            data::_3d::Point3SPtr p_box_min =
                    data::_3d::KernelFactory::createPoint3(-10.0, -10.0, -10.0);
            data::_3d::Point3SPtr p_box_max =
                    data::_3d::KernelFactory::createPoint3(10.0, 10.0, 10.0);
            algo::_3d::PolyhedronTransformation::translateNscale(
                    polyhedron, p_box_min, p_box_max);
        }
        if (rand_move_points) {
            std::cout << "Points will be moved randomly. "
                 << "(rand_move_points_range=" << rand_move_points_range << ")" << std::endl;
            algo::_3d::PolyhedronTransformation::randMovePoints(polyhedron, rand_move_points_range);
            std::string description = polyhedron->getDescription();
            polyhedron = algo::_3d::SimpleStraightSkel::shiftFacets(polyhedron, 0.0);
            polyhedron->clearData();
            polyhedron->setDescription(description);
        }
        if (translate_and_scale_view) {
            data::_3d::Point3SPtr p_box_min =
                    algo::_3d::PolyhedronTransformation::boundingBoxMin(polyhedron);
            data::_3d::Point3SPtr p_box_max =
                    algo::_3d::PolyhedronTransformation::boundingBoxMax(polyhedron);
            float scale_min = std::numeric_limits<float>::max();
            for (unsigned int i = 0; i < 3; i++) {
                translate[i] = -((*p_box_max)[i] + (*p_box_min)[i])/2.0;
                float scale_cur = 20.0/((*p_box_max)[i] - (*p_box_min)[i]);
                if (scale_cur < scale_min) {
                    scale_min = scale_cur;
                }
                scale = scale_min;
            }
        }
        if (!polyhedron->isConsistent()) {
            std::cout << "Warning: Polyhedron with PolyhedronID=" << id
                 << " is not consistent." << std::endl;
        }
        DEBUG_VAR(polyhedron->toString());
    }

    // create OpenGL window
    bool no_window = isSet("--no-window", argc, argv);
    algo::ControllerSPtr controller;
    ui::gl::MainOpenGLWindowSPtr window;
    if (!no_window) {
        controller = algo::Controller::create();
        controller->togglePause();
        int width = config->getInt("ui_gl_MainOpenGLWindow", "width");
        if (width == 0) {
            width = 800;
        }
        int height = config->getInt("ui_gl_MainOpenGLWindow", "height");
        if (height == 0) {
            height = 600;
        }
        window = ui::gl::MainOpenGLWindow::create(argc, argv, width, height, controller);
    }

    std::list<double> save_offsets;
    const char* chr_save_offsets = getOption("--save-offsets", argc, argv);
    if (chr_save_offsets) {
        save_offsets = parseCSV(chr_save_offsets);
    }

    // run algorithm
    // has to be in this scope otherwise the algorithm will be destroyed
    ThreadSPtr thread_algo;
    ThreadSPtr thread_window;
    algo::_2d::SimpleStraightSkelSPtr algoskel2d;
    algo::_3d::SimpleStraightSkelSPtr algoskel3d;
    if (num_dims == 2) {
        algoskel2d = algo::_2d::SimpleStraightSkel::create(
                polygon, controller);
        if (window) {
            window->setPolygon(polygon);
            if (skel2d) {
                window->setSkel2d(skel2d);
            }
            thread_algo = algoskel2d->startThread();
        } else {
            algoskel2d->run();
        }
    } else if (num_dims == 3) {
        algoskel3d = algo::_3d::SimpleStraightSkel::create(
                polyhedron, controller, save_offsets);
        if (window) {
            window->setPolyhedron(polyhedron);
            if (skel3d) {
                window->setSkel3d(skel3d);
            }
            thread_algo = algoskel3d->startThread();
        } else {
            algoskel3d->run();
        }
    }
    if (window) {
        if (translate_and_scale_view) {
            window->setTranslate(translate);
            window->setScale(scale);
        }
        thread_window = window->startThread();
        thread_algo->join();
    }

    bool save = isSet("--save", argc, argv);
    if (save) {
        if (num_dims == 2) {
            skel2d = algoskel2d->getResult();
            if (skel2d) {
                db::_2d::DAOFactory::getPolygonDAO()->insert(polygon);
                db::_2d::DAOFactory::getStraightSkeletonDAO()->insert(skel2d);
                std::cout << "SkelID=" << skel2d->getID() << std::endl;
            }
        } else if (num_dims == 3) {
            skel3d = algoskel3d->getResult();
            if (skel3d) {
                db::_3d::DAOFactory::getPolyhedronDAO()->insert(polyhedron);
                db::_3d::DAOFactory::getStraightSkeletonDAO()->insert(skel3d);
                std::cout << "SkelID=" << skel3d->getID() << std::endl;
            }
        }
    }

    algo::_2d::SkelMeshGeneratorSPtr algomesh2d;
    bool create_mesh = isSet("--mesh", argc, argv);
    if (create_mesh) {
        controller->wait();
        if (num_dims == 2) {
            skel2d = algoskel2d->getResult();
            algomesh2d = algo::_2d::SkelMeshGenerator::create(skel2d, controller);
            if (window) {
                window->setMesh2d(algomesh2d->getResult());
                thread_algo = algomesh2d->startThread();
            } else {
                algomesh2d->run();
            }
        }
    }

    if (isSet("--check-graph", argc, argv) && num_dims == 3) {
        algo::_3d::GraphCheckerSPtr graphchecker = algo::_3d::GraphChecker::create();
        skel3d = algoskel3d->getResult();
        graphchecker->check(skel3d);
    }

    if (window) {
        thread_window->join();
    }

    return EXIT_SUCCESS;
}
