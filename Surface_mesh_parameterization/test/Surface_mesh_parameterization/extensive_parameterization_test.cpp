// extensive_parameterization_test.cpp


// ----------------------------------------------------------------------------
// USAGE EXAMPLES
// ----------------------------------------------------------------------------

//----------------------------------------------------------
// Test all parameterization methods and all solvers
// No output
// Input files are .off
//----------------------------------------------------------
// extensive_parameterization_test mesh1.off mesh2.off...


// CGAL
#include <CGAL/basic.h> // include basic.h before testing #defines
#include <CGAL/Cartesian.h>
#include <CGAL/Timer.h>

// This package
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
#ifdef CGAL_EIGEN3_ENABLED
    #include <CGAL/Eigen_solver_traits.h>
#endif

// This test
#include "Polyhedron_ex.h"
#include "Mesh_cutter.h"
#include "Parameterization_polyhedron_adaptor_ex.h"

// STL stuff
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <cassert>


// ----------------------------------------------------------------------------
// Private types
// ----------------------------------------------------------------------------

// Mesh true type and Surface_mesh_parameterization adaptors
typedef Polyhedron_ex                                       Polyhedron;
typedef Parameterization_polyhedron_adaptor_ex              Parameterization_polyhedron_adaptor;
typedef CGAL::Parameterization_mesh_patch_3<Parameterization_polyhedron_adaptor>
                                                            Mesh_patch_polyhedron;

// Parameterizer for this kind of mesh
typedef CGAL::Parameterizer_traits_3<Mesh_patch_polyhedron> Parameterizer;

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
            assert(nb_borders == 0);
            cutter.cut(seamingBackbone);   // simple cut
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
        // 3) check that the seamingBackbone is a 2-way list
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


// ----------------------------------------------------------------------------
// main()
// ----------------------------------------------------------------------------

int main(int argc, char * argv[])
{
    std::cerr << "PARAMETERIZATION" << std::endl;
    std::cerr << "Test all parameterization methods and all solvers" << std::endl;
    std::cerr << "No output" << std::endl;

    //***************************************
    // decode parameters
    //***************************************

    if (argc-1 == 0)
    {
        std::cerr << "Usage: " << argv[0] << " input_file1.off input_file2.obj ..." << std::endl;
        return(EXIT_FAILURE);
    }

    // Accumulated errors
    int accumulated_fatal_err = EXIT_SUCCESS;

    // Parameterize each input file and accumulate errors
    for (int arg_index = 1; arg_index <= argc-1; arg_index++)
    {
        std::cerr << std::endl << std::endl;

        // File name is:
        const char* input_filename  = argv[arg_index];

        //***************************************
        // Read the mesh
        //***************************************

        CGAL::Timer task_timer;
        task_timer.start();

        // Read the mesh
        std::ifstream stream(input_filename);
        Polyhedron mesh;
        stream >> mesh;
        if(!stream || !mesh.is_valid() || mesh.empty())
        {
            std::cerr << "Error: cannot read OFF file " << input_filename << std::endl;
            accumulated_fatal_err = EXIT_FAILURE;
            continue;
        }
        std::cerr << "Read file " << input_filename << ": "
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
            // this is not a bug => do not set accumulated_fatal_err
            continue;
        }
        //
        // 2) Create adaptor that virtually "cuts" a patch in a Polyhedron_ex mesh
        Mesh_patch_polyhedron   mesh_patch(mesh_adaptor, seam.begin(), seam.end());
        if (!mesh_patch.is_valid())
        {
            std::cerr << "Input mesh not supported: non manifold shape or invalid cutting" << std::endl;
            // this is not a bug => do not set accumulated_fatal_err
            continue;
        }

        std::cerr << "Mesh cutting: " << task_timer.time() << " seconds." << std::endl << std::endl;
        task_timer.reset();

        //***************************************
        // Tutte Barycentric Mapping parameterization
        // with square uniform border parameterization
        // OpenNL solver
        //***************************************

        Parameterizer::Error_code err;

        std::cerr << "Tutte Barycentric Mapping parameterization" << std::endl;
        std::cerr << "with square uniform border parameterization" << std::endl;
        std::cerr << "OpenNL solver" << std::endl;

        err = CGAL::parameterize(
            mesh_patch,
            CGAL::Barycentric_mapping_parameterizer_3<
                Mesh_patch_polyhedron,
                CGAL::Square_border_uniform_parameterizer_3<Mesh_patch_polyhedron>,
                OpenNL::DefaultLinearSolverTraits<double>
            >());
        switch(err) {
        case Parameterizer::OK: // Success
            break;
        case Parameterizer::ERROR_EMPTY_MESH: // Input mesh not supported
        case Parameterizer::ERROR_NON_TRIANGULAR_MESH:
        case Parameterizer::ERROR_NO_TOPOLOGICAL_DISC:
        case Parameterizer::ERROR_BORDER_TOO_SHORT:
            std::cerr << "Input mesh not supported: " << Parameterizer::get_error_message(err) << std::endl;
            // this is not a bug => do not set accumulated_fatal_err
            break;
        default: // Error
            std::cerr << "Error: " << Parameterizer::get_error_message(err) << std::endl;
            accumulated_fatal_err = EXIT_FAILURE;
            break;
        };

        std::cerr << "Parameterization: " << task_timer.time() << " seconds." << std::endl << std::endl;
        task_timer.reset();

        //***************************************
        // Floater Mean Value Coordinates parameterization
        // with circular arc length border parameterization
        // OpenNL solver
        //***************************************

        std::cerr << "Floater Mean Value Coordinates parameterization" << std::endl;
        std::cerr << "with circular arc length border parameterization" << std::endl;
        std::cerr << "OpenNL solver" << std::endl;

        err = CGAL::parameterize(
            mesh_patch,
            CGAL::Mean_value_coordinates_parameterizer_3<
                Mesh_patch_polyhedron,
                CGAL::Circular_border_arc_length_parameterizer_3<Mesh_patch_polyhedron>,
                OpenNL::DefaultLinearSolverTraits<double>
            >());
        switch(err) {
        case Parameterizer::OK: // Success
            break;
        case Parameterizer::ERROR_EMPTY_MESH: // Input mesh not supported
        case Parameterizer::ERROR_NON_TRIANGULAR_MESH:
        case Parameterizer::ERROR_NO_TOPOLOGICAL_DISC:
        case Parameterizer::ERROR_BORDER_TOO_SHORT:
            std::cerr << "Input mesh not supported: " << Parameterizer::get_error_message(err) << std::endl;
            // this is not a bug => do not set accumulated_fatal_err
            break;
        default: // Error
            std::cerr << "Error: " << Parameterizer::get_error_message(err) << std::endl;
            accumulated_fatal_err = EXIT_FAILURE;
            break;
        };

        std::cerr << "Parameterization: " << task_timer.time() << " seconds." << std::endl << std::endl;
        task_timer.reset();

        //***************************************
        // Least Squares Conformal Maps parameterization
        // OpenNL solver
        //***************************************

        std::cerr << "Least Squares Conformal Maps parameterization" << std::endl;
        std::cerr << "OpenNL solver" << std::endl;

        err = CGAL::parameterize(
            mesh_patch,
            CGAL::LSCM_parameterizer_3<
                Mesh_patch_polyhedron,
                CGAL::Two_vertices_parameterizer_3<Mesh_patch_polyhedron>,
                OpenNL::SymmetricLinearSolverTraits<double>
            >());
        switch(err) {
        case Parameterizer::OK: // Success
            break;
        case Parameterizer::ERROR_EMPTY_MESH: // Input mesh not supported
        case Parameterizer::ERROR_NON_TRIANGULAR_MESH:
        case Parameterizer::ERROR_NO_TOPOLOGICAL_DISC:
        case Parameterizer::ERROR_BORDER_TOO_SHORT:
            std::cerr << "Input mesh not supported: " << Parameterizer::get_error_message(err) << std::endl;
            // this is not a bug => do not set accumulated_fatal_err
            break;
        default: // Error
            std::cerr << "Error: " << Parameterizer::get_error_message(err) << std::endl;
            accumulated_fatal_err = EXIT_FAILURE;
            break;
        };

        std::cerr << "Parameterization: " << task_timer.time() << " seconds." << std::endl << std::endl;
        task_timer.reset();

#ifdef CGAL_EIGEN3_ENABLED

        //***************************************
        // Discrete Conformal Map parameterization
        // with circular arc length border parameterization
        // Eigen solver (if installed)
        //***************************************

        std::cerr << "Discrete Conformal Map parameterization" << std::endl;
        std::cerr << "with circular arc length border parameterization" << std::endl;
        std::cerr << "Eigen solver" << std::endl;

        err = CGAL::parameterize(
            mesh_patch,
            CGAL::Discrete_conformal_map_parameterizer_3<
                Mesh_patch_polyhedron,
                CGAL::Circular_border_arc_length_parameterizer_3<Mesh_patch_polyhedron>,
                CGAL::Eigen_solver_traits<>
            >());
        switch(err) {
        case Parameterizer::OK: // Success
            break;
        case Parameterizer::ERROR_EMPTY_MESH: // Input mesh not supported
        case Parameterizer::ERROR_NON_TRIANGULAR_MESH:
        case Parameterizer::ERROR_NO_TOPOLOGICAL_DISC:
        case Parameterizer::ERROR_BORDER_TOO_SHORT:
            std::cerr << "Input mesh not supported: " << Parameterizer::get_error_message(err) << std::endl;
            // this is not a bug => do not set accumulated_fatal_err
            break;
        default: // Error
            std::cerr << "Error: " << Parameterizer::get_error_message(err) << std::endl;
            accumulated_fatal_err = EXIT_FAILURE;
            break;
        };

        std::cerr << "Parameterization: " << task_timer.time() << " seconds." << std::endl << std::endl;
        task_timer.reset();

        //***************************************
        // Discrete Authalic Parameterization
        // with square arc length border parameterization
        // Eigen solver (if installed)
        //***************************************

        std::cerr << "Discrete Authalic Parameterization" << std::endl;
        std::cerr << "with square arc length border parameterization" << std::endl;
        std::cerr << "Eigen solver" << std::endl;

        err = CGAL::parameterize(
            mesh_patch,
            CGAL::Discrete_authalic_parameterizer_3<
	    Mesh_patch_polyhedron,
	    CGAL::Square_border_arc_length_parameterizer_3<Mesh_patch_polyhedron> >()
				 );
        switch(err) {
        case Parameterizer::OK: // Success
            break;
        case Parameterizer::ERROR_EMPTY_MESH: // Input mesh not supported
        case Parameterizer::ERROR_NON_TRIANGULAR_MESH:
        case Parameterizer::ERROR_NO_TOPOLOGICAL_DISC:
        case Parameterizer::ERROR_BORDER_TOO_SHORT:
            std::cerr << "Input mesh not supported: " << Parameterizer::get_error_message(err) << std::endl;
            // this is not a bug => do not set accumulated_fatal_err
            break;
        default: // Error
            std::cerr << "Error: " << Parameterizer::get_error_message(err) << std::endl;
            accumulated_fatal_err = EXIT_FAILURE;
            break;
        };

        std::cerr << "Parameterization: " << task_timer.time() << " seconds." << std::endl << std::endl;
        task_timer.reset();

        //***************************************
        // Least Squares Conformal Maps parameterization
        // Eigen solver (if installed)
        //***************************************

        std::cerr << "Least Squares Conformal Maps parameterization" << std::endl;
        std::cerr << "Eigen solver" << std::endl;
        
        typedef CGAL::Eigen_sparse_matrix<double>::EigenType EigenMatrix;
        err = CGAL::parameterize(
            mesh_patch,
            CGAL::LSCM_parameterizer_3<
                Mesh_patch_polyhedron,
                CGAL::Two_vertices_parameterizer_3<Mesh_patch_polyhedron>,
                CGAL::Eigen_solver_traits< Eigen::ConjugateGradient<EigenMatrix> >
            >());
        switch(err) {
        case Parameterizer::OK: // Success
            break;
        case Parameterizer::ERROR_EMPTY_MESH: // Input mesh not supported
        case Parameterizer::ERROR_NON_TRIANGULAR_MESH:
        case Parameterizer::ERROR_NO_TOPOLOGICAL_DISC:
        case Parameterizer::ERROR_BORDER_TOO_SHORT:
            std::cerr << "Input mesh not supported: " << Parameterizer::get_error_message(err) << std::endl;
            // this is not a bug => do not set accumulated_fatal_err
            break;
        default: // Error
            std::cerr << "Error: " << Parameterizer::get_error_message(err) << std::endl;
            accumulated_fatal_err = EXIT_FAILURE;
            break;
        };

        std::cerr << "Parameterization: " << task_timer.time() << " seconds." << std::endl << std::endl;
        task_timer.reset();

#else

        std::cerr << "Skip EIGEN tests as EIGEN is not installed" << std::endl << std::endl;
        // this is not a bug => do not set accumulated_fatal_err

#endif // CGAL_USE_EIGEN

    } // for each input file

    // Return accumulated fatal error
    std::cerr << "Tool returned " << accumulated_fatal_err << std::endl;
    return accumulated_fatal_err;
}


