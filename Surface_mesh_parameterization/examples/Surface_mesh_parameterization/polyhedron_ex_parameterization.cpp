// ----------------------------------------------------------------------------
// USAGE EXAMPLES
// ----------------------------------------------------------------------------

//----------------------------------------------------------
// Discrete Conformal Map parameterization
// circle border
// OpenNL solver
// output is a eps map
// input file is mesh.off
//----------------------------------------------------------
// polyhedron_ex_parameterization -t conformal -b circle mesh.off mesh.eps

//----------------------------------------------------------
// floater parameterization
// square border
// TAUCS solver
// output is a eps map
// input file is mesh.off
//----------------------------------------------------------
// polyhedron_ex_parameterization -t floater -b square -s taucs mesh.off mesh.eps

//----------------------------------------------------------
// Least Squares Conformal Maps parameterization
// two pinned vertices (automatically picked)
// OpenNL solver
// output is a .obj
// input file is mesh.off
//----------------------------------------------------------
// polyhedron_ex_parameterization -t lscm -b 2pts mesh.off mesh.obj


#include <CGAL/Cartesian.h>
#include <CGAL/Timer.h>
#include <CGAL/parameterize.h>
#include <CGAL/Parameterization_mesh_patch_3.h>
#include <CGAL/Circular_border_parameterizer_3.h>
#include <CGAL/Square_border_parameterizer_3.h>
#include <CGAL/Two_vertices_parameterizer_3.h>
#include <CGAL/Barycentric_mapping_parameterizer_3.h>
#include <CGAL/Discrete_conformal_map_parameterizer_3.h>
#include <CGAL/Discrete_authalic_parameterizer_3.h>
#include <CGAL/Mean_value_coordinates_parameterizer_3.h>
#include <CGAL/LSCM_parameterizer_3.h>
#include <CGAL/Parameterization_mesh_feature_extractor.h>

#include <CGAL/OpenNL/linear_solver.h>
#ifdef CGAL_USE_TAUCS
    #include <CGAL/Taucs_solver_traits.h>
#endif

#include "Polyhedron_ex.h"
#include "Mesh_cutter.h"
#include "Parameterization_polyhedron_adaptor_ex.h"

#include <iostream>
#include <string.h>
#include <ctype.h>
#include <fstream>
#include <cassert>

#ifdef CGAL_USE_BOOST_PROGRAM_OPTIONS
    #include <boost/program_options.hpp>
    namespace po = boost::program_options;
#endif


// ----------------------------------------------------------------------------
// Private types
// ----------------------------------------------------------------------------

typedef Polyhedron_ex                                       Polyhedron;

// Mesh adaptors
typedef Parameterization_polyhedron_adaptor_ex              Parameterization_polyhedron_adaptor;
typedef CGAL::Parameterization_mesh_patch_3<Parameterization_polyhedron_adaptor>
                                                            Mesh_patch_polyhedron;

// Type describing a border or seam as a vertex list
typedef std::list<Parameterization_polyhedron_adaptor::Vertex_handle>
                                                            Seam;


// ----------------------------------------------------------------------------
// Private functions
// ----------------------------------------------------------------------------

// Cut the mesh to make it homeomorphic to a disk
// or extract a region homeomorphic to a disc.
// Return the border of this region (empty on error)
//
// CAUTION:
// This method is provided "as is". It is very buggy and simply part of this example.
// Developers using this package should implement a more robust cut algorithm!
static Seam cut_mesh(Parameterization_polyhedron_adaptor& mesh_adaptor)
{
    // Helper class to compute genus or extract borders
    typedef CGAL::Parameterization_mesh_feature_extractor<Parameterization_polyhedron_adaptor_ex>
                                            Mesh_feature_extractor;
    typedef Mesh_feature_extractor::Border  Border;
    typedef Mesh_cutter::Backbone           Backbone;

    Seam seam;              // returned list

    // Get refererence to Polyhedron_3 mesh
    Polyhedron& mesh = mesh_adaptor.get_adapted_mesh();

    // Extract mesh borders and compute genus
    Mesh_feature_extractor feature_extractor(mesh_adaptor);
    int nb_borders = feature_extractor.get_nb_borders();
    int genus = feature_extractor.get_genus();

    // If mesh is a topological disk
    if (genus == 0 && nb_borders > 0)
    {
        // Pick the longest border
        seam = feature_extractor.get_longest_border();
    }
    else // if mesh is *not* a topological disk, create a virtual cut
    {
        Backbone seamingBackbone;           // result of cutting
        Backbone::iterator he;

        // Compute a cutting path that makes the mesh a "virtual" topological disk
        mesh.compute_facet_centers();
        Mesh_cutter cutter(mesh);
        if (genus == 0)
        {
            // no border, we need to cut the mesh
            assert (nb_borders == 0);
            cutter.cut(seamingBackbone);    // simple cut
        }
        else // genus > 0 -> cut the mesh
        {
            cutter.cut_genus(seamingBackbone);
        }

        // The Mesh_cutter class is quite buggy
        // => we check that seamingBackbone is valid
        //
        // 1) Check that seamingBackbone is not empty
        if (seamingBackbone.begin() == seamingBackbone.end())
            return seam;                    // return empty list
        //
        // 2) Check that seamingBackbone is a loop and
        //    count occurences of seam halfedges
        mesh.tag_halfedges(0);              // Reset counters
        for (he = seamingBackbone.begin(); he != seamingBackbone.end(); he++)
        {
            // Get next halfedge iterator (looping)
            Backbone::iterator next_he = he;
            next_he++;
            if (next_he == seamingBackbone.end())
                next_he = seamingBackbone.begin();

            // Check that seamingBackbone is a loop: check that
            // end of current HE == start of next one
            if ((*he)->vertex() != (*next_he)->opposite()->vertex())
                return seam;                // return empty list

            // Increment counter (in "tag" field) of seam halfedges
            (*he)->tag( (*he)->tag()+1 );
        }
        //
        // 3) check that the seamingBackbone is a two-way list
        for (he = seamingBackbone.begin(); he != seamingBackbone.end(); he++)
        {
            // Counter of halfedge and opposite halfedge must be 1
            if ((*he)->tag() != 1 || (*he)->opposite()->tag() != 1)
                return seam;                // return empty list
        }

        // Convert list of halfedges to a list of vertices
        for (he = seamingBackbone.begin(); he != seamingBackbone.end(); he++)
            seam.push_back((*he)->vertex());
    }

    return seam;
}

// Call appropriate parameterization method based on command line parameters
template<
    class ParameterizationMesh_3,   // 3D surface
    class GeneralSparseLinearAlgebraTraits_d,
                                    // Traits class to solve a general sparse linear system
    class SymmetricSparseLinearAlgebraTraits_d
                                    // Traits class to solve a symmetric sparse linear system
>
typename CGAL::Parameterizer_traits_3<ParameterizationMesh_3>::Error_code
parameterize(ParameterizationMesh_3& mesh,  // Mesh parameterization adaptor
             const std::string& type,              // type of parameterization (see usage)
             const std::string& border)            // type of border parameterization (see usage)
{
    typename CGAL::Parameterizer_traits_3<ParameterizationMesh_3>::Error_code err;

    if ( (type == std::string("floater"))  && (border == std::string("circle")) )
    {
        err = CGAL::parameterize(
            mesh,
            CGAL::Mean_value_coordinates_parameterizer_3<
                ParameterizationMesh_3,
                CGAL::Circular_border_arc_length_parameterizer_3<ParameterizationMesh_3>,
                GeneralSparseLinearAlgebraTraits_d
            >());
    }
    else if ( (type == std::string("floater")) && (border == std::string("square")) )
    {
        err = CGAL::parameterize(
            mesh,
            CGAL::Mean_value_coordinates_parameterizer_3<
                ParameterizationMesh_3,
                CGAL::Square_border_arc_length_parameterizer_3<ParameterizationMesh_3>,
                GeneralSparseLinearAlgebraTraits_d
            >());
    }
    else if ( (type == std::string("barycentric")) && (border == std::string("circle")) )
    {
        err = CGAL::parameterize(
            mesh,
            CGAL::Barycentric_mapping_parameterizer_3<
                ParameterizationMesh_3,
                CGAL::Circular_border_uniform_parameterizer_3<ParameterizationMesh_3>,
                GeneralSparseLinearAlgebraTraits_d
            >());
    }
    else if ( (type == std::string("barycentric")) && (border == std::string("square")) )
    {
        err = CGAL::parameterize(
            mesh,
            CGAL::Barycentric_mapping_parameterizer_3<
                ParameterizationMesh_3,
                CGAL::Square_border_uniform_parameterizer_3<ParameterizationMesh_3>,
                GeneralSparseLinearAlgebraTraits_d
            >());
    }
    else if ( (type == std::string("conformal")) && (border == std::string("circle")) )
    {
        err = CGAL::parameterize(
            mesh,
            CGAL::Discrete_conformal_map_parameterizer_3<
                ParameterizationMesh_3,
                CGAL::Circular_border_arc_length_parameterizer_3<ParameterizationMesh_3>,
                GeneralSparseLinearAlgebraTraits_d
            >());
    }
    else if ( (type == std::string("conformal")) && (border == std::string("square")) )
    {
        err = CGAL::parameterize(
            mesh,
            CGAL::Discrete_conformal_map_parameterizer_3<
                ParameterizationMesh_3,
                CGAL::Square_border_arc_length_parameterizer_3<ParameterizationMesh_3>,
                GeneralSparseLinearAlgebraTraits_d
            >());
    }
    else if ( (type == std::string("authalic")) && (border == std::string("circle")) )
    {
        err = CGAL::parameterize(
            mesh,
            CGAL::Discrete_authalic_parameterizer_3<
                ParameterizationMesh_3,
                CGAL::Circular_border_arc_length_parameterizer_3<ParameterizationMesh_3>,
                GeneralSparseLinearAlgebraTraits_d
            >());
    }
    else if ( (type == std::string("authalic")) && (border == std::string("square")) )
    {
        err = CGAL::parameterize(
            mesh,
            CGAL::Discrete_authalic_parameterizer_3<
                ParameterizationMesh_3,
                CGAL::Square_border_arc_length_parameterizer_3<ParameterizationMesh_3>,
                GeneralSparseLinearAlgebraTraits_d
            >());
    }
    else if ( (type == std::string("lscm")) && (border == std::string("2pts")) )
    {
        err = CGAL::parameterize(
            mesh,
            CGAL::LSCM_parameterizer_3<
                ParameterizationMesh_3,
                CGAL::Two_vertices_parameterizer_3<ParameterizationMesh_3>,
                SymmetricSparseLinearAlgebraTraits_d
            >());
    }
    else
    {
        std::cerr << "Error: invalid parameters combination " << type << " + " << border << std::endl;
        err = CGAL::Parameterizer_traits_3<ParameterizationMesh_3>::ERROR_WRONG_PARAMETER;
    }

    return err;
}


// ----------------------------------------------------------------------------
// main()
// ----------------------------------------------------------------------------

#ifdef CGAL_USE_BOOST_PROGRAM_OPTIONS
int main(int argc, char * argv[])
#else
int main()
#endif
{
    CGAL::Timer total_timer;
    total_timer.start();

    std::cerr << "PARAMETERIZATION" << std::endl;

    //***************************************
    // Read options on the command line
    //***************************************

    std::string type;               // default: Floater param
    std::string border;             // default: circular border param.
    std::string solver;             // default: OpenNL solver
    std::string input;              // required
    std::string output;             // default: out.eps
    try
    {
#ifdef CGAL_USE_BOOST_PROGRAM_OPTIONS
        po::options_description desc("Allowed options");
        desc.add_options()
            ("help,h", "prints this help message")
            ("type,t", po::value<std::string>(&type)->default_value("floater"),
            "parameterization method: floater, conformal, barycentric, authalic or lscm")
            ("border,b", po::value<std::string>(&border)->default_value("circle"),
            "border shape: circle, square or 2pts (lscm only)")
            ("solver,s", po::value<std::string>(&solver)->default_value("opennl"),
            "solver: opennl or taucs")
            ("input,i", po::value<std::string>(&input)->default_value(""),
            "input mesh (OFF)")
            ("output,o", po::value<std::string>(&output)->default_value("out.eps"),
            "output file (EPS or OBJ)")
            ;

        po::positional_options_description p;
        p.add("input", 1);
        p.add("output", 1);

        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
        po::notify(vm);

        if (vm.count("help")) {
            std::cout << desc << "\n";
            return 1;
        }
#else
        std::cerr << "Command-line options require Boost.ProgramOptions" << std::endl;
        std::cerr << "Use hard-coded options" << std::endl;
        border = "square";
        type = "floater";
        solver = "opennl";
        input = "data/rotor.off";
        output = "rotor_floater_square_opennl_parameterized.obj";
#endif
    }
    catch(std::exception& e) {
      std::cerr << "error: " << e.what() << "\n";
      return 1;
    }
    catch(...) {
      std::cerr << "Exception of unknown type!\n";
    }

    //***************************************
    // Read the mesh
    //***************************************

    CGAL::Timer task_timer;
    task_timer.start();

    // Read the mesh
    std::ifstream stream(input.c_str());
    Polyhedron mesh;
    stream >> mesh;
    if(!stream || !mesh.is_valid() || mesh.empty())
    {
        std::cerr << "Error: cannot read OFF file " << input << std::endl;
        return EXIT_FAILURE;
    }

    std::cerr << "Read file " << input << ": "
              << task_timer.time() << " seconds "
              << "(" << mesh.size_of_facets() << " facets, "
              << mesh.size_of_vertices() << " vertices)" << std::endl;
    task_timer.reset();

    //***************************************
    // Create mesh adaptor
    //***************************************

    // The Surface_mesh_parameterization package needs an adaptor to handle Polyhedron_ex meshes
    Parameterization_polyhedron_adaptor mesh_adaptor(mesh);

    // The parameterization methods support only meshes that
    // are topological disks => we need to compute a cutting path
    // that makes the mesh a "virtual" topological disk
    //
    // 1) Cut the mesh
    Seam seam = cut_mesh(mesh_adaptor);
    if (seam.empty())
    {
        std::cerr << "Input mesh not supported: the example cutting algorithm is too simple to cut this shape" << std::endl;
        return EXIT_FAILURE;
    }
    //
    // 2) Create adaptor that virtually "cuts" a patch in a Polyhedron_ex mesh
    Mesh_patch_polyhedron   mesh_patch(mesh_adaptor, seam.begin(), seam.end());
    if (!mesh_patch.is_valid())
    {
        std::cerr << "Input mesh not supported: non manifold shape or invalid cutting" << std::endl;
        return EXIT_FAILURE;
    }

    std::cerr << "Mesh cutting: " << task_timer.time() << " seconds." << std::endl;
    task_timer.reset();

    //***************************************
    // switch parameterization
    //***************************************

    std::cerr << "Parameterization..." << std::endl;

    // Defines the error codes
    typedef CGAL::Parameterizer_traits_3<Mesh_patch_polyhedron> Parameterizer;
    Parameterizer::Error_code err;

    if (solver == std::string("opennl"))
    {
        err = parameterize<Mesh_patch_polyhedron,
                           OpenNL::DefaultLinearSolverTraits<double>,
                           OpenNL::SymmetricLinearSolverTraits<double>
                          >(mesh_patch, type, border);
    }
    else if (solver == std::string("taucs"))
    {
#ifdef CGAL_USE_TAUCS
        err = parameterize<Mesh_patch_polyhedron,
                           CGAL::Taucs_solver_traits<double>,
                           CGAL::Taucs_symmetric_solver_traits<double>
                          >(mesh_patch, type, border);
#else
        std::cerr << "Error: TAUCS is not installed" << std::endl;
        err = Parameterizer::ERROR_WRONG_PARAMETER;
#endif
    }
    else
    {
        std::cerr << "Error: invalid solver parameter " << solver << std::endl;
        err = Parameterizer::ERROR_WRONG_PARAMETER;
    }

    // Report errors
    switch(err) {
    case Parameterizer::OK: // Success
        break;
    case Parameterizer::ERROR_EMPTY_MESH: // Input mesh not supported
    case Parameterizer::ERROR_NON_TRIANGULAR_MESH:
    case Parameterizer::ERROR_NO_TOPOLOGICAL_DISC:
    case Parameterizer::ERROR_BORDER_TOO_SHORT:
        std::cerr << "Input mesh not supported: " << Parameterizer::get_error_message(err) << std::endl;
        return EXIT_FAILURE;
        break;
    default: // Error
        std::cerr << "Error: " << Parameterizer::get_error_message(err) << std::endl;
        return EXIT_FAILURE;
        break;
    };

    std::cerr << "Parameterization: " << task_timer.time() << " seconds." << std::endl;
    task_timer.reset();

    //***************************************
    // Output
    //***************************************

    // get output file's extension
    std::string extension = output.substr(output.find_last_of('.'));

    // Save mesh
    if (extension == ".eps" || extension == ".EPS")
    {
        // write Postscript file
        if ( ! mesh.write_file_eps(output.c_str()) )
        {
            std::cerr << "Error: cannot write file " << output << std::endl;
            return EXIT_FAILURE;
        }
    }
    else if (extension == ".obj" || extension == ".OBJ")
    {
        // write Wavefront obj file
        if ( ! mesh.write_file_obj(output.c_str()) )
        {
            std::cerr << "Error: cannot write file " << output << std::endl;
            return EXIT_FAILURE;
        }
    }
    else
    {
        std::cerr << "Error: output format not supported" << output << std::endl;
        err = Parameterizer::ERROR_WRONG_PARAMETER;
        return EXIT_FAILURE;
    }

    std::cerr << "Write file " << output << ": "
              << task_timer.time() << " seconds " << std::endl;

    return EXIT_SUCCESS;
}
