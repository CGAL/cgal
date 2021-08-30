#define BOOST_FILESYSTEM_VERSION 2

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/refine_mesh_3.h>

/* OUTPUT */
#include <CGAL/IO/File_medit.h>

/* tools */
#include <sstream>
#include <stdlib.h>
#include <algorithm>
#include <memory>

#include "thread_queue.h"

/* OPTIONS and PARAMETERS */
#include <boost/program_options.hpp>
namespace po = boost::program_options;

using namespace CGAL::parameters; //to avoid verbose function and named parameters call

/* FILE SYSTEM */
#include <iostream>
#include <boost/filesystem.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/convenience.hpp>
namespace fs = boost::filesystem;


struct K: public CGAL::Exact_predicates_inexact_constructions_kernel {};


template <class Domain>
class Domain_builder {};


template<typename C3t3> void save_histogram(const std::string& filename,const C3t3& c3t3);
template<typename MeshCriteria> MeshCriteria get_parameters(const std::string& param_line,
                                                            po::variables_map& vm);

template <class C3T3, class Domain> struct Optimizer;


template <class C3T3, class MeshCriteria, class Domain>
struct Mesher
{
  Mesher(const std::shared_ptr<Domain_builder<Domain> >& pdomain_builder,
         const int mesh_nb,
         const std::string& filename,
         const std::string& output,
         const std::vector<std::string>& lines,
         MessageQueue& msg_queue)
    : pdomain_builder_(pdomain_builder)
    , mesh_nb_(mesh_nb)
    , filename_(filename)
    , output_prefix_(output)
    , lines_(lines)
    , msg_queue_(&msg_queue) {}


  void operator()(unsigned int nb, const std::string& str)
  {
    if ( (unsigned int)PRINT_ERROR_MSG == nb ) return print(str);

    std::vector<std::string>::iterator lit = lines_.begin();

    std::stringstream cout_loc;
    cout_loc << "[" << nb << "] Launch mesher[" << filename_ << "_" << mesh_nb_ << "] with: "
             << *lit << std::endl;
    std::cout << cout_loc.str();

    //Prog out
    std::ofstream file_out(progout_filename().c_str());
    file_out << *lit << std::endl;

    po::variables_map vm_p_mesh;
    MeshCriteria mcp = get_parameters<MeshCriteria>(*lit++, vm_p_mesh);

    // Mesh generation
    if ( ! vm_p_mesh.count("mesh") )
    {
      std::cerr << "[" << nb << "] ERROR: Bad mesh args: --mesh is missing\n";
      file_out << "ERROR: Bad mesh args: --mesh is missing\n";
      return;
    }

    CGAL::Timer timer;
    timer.start();

    // we keep c3t3 between lines
    std::shared_ptr<C3T3> pc3t3_save (new C3T3());

    // Generate Mesh
    file_out << "Generate mesh...";
    std::flush(file_out);
    // Initialize c3t3
    if ( !vm_p_mesh.count("off_vertices") )
    {
      *pc3t3_save = CGAL::make_mesh_3<C3T3>(pdomain_builder_->domain(), MeshCriteria(), no_exude(), no_perturb());
    }
    else
    {
      pc3t3_save->insert_surface_points(pdomain_builder_->points_begin(),
                                        pdomain_builder_->points_end(),
                                        pdomain_builder_->points_index());
    }

    // Mesh
    CGAL::refine_mesh_3<C3T3>(*pc3t3_save, pdomain_builder_->domain(), mcp, no_exude(), no_perturb());
    file_out << "done (" << timer.time() << "s)   [" << pc3t3_save->triangulation().number_of_vertices()
            << " vertices, " << pc3t3_save->number_of_facets() << " facets, " << pc3t3_save->number_of_cells() << " cells]\n";

    //save mesh
    file_out << "Save mesh...";
    std::flush(file_out);
    std::stringstream output_filename;
    output_filename << output_prefix_ << "_" << mesh_nb_ << "-out.mesh";
    std::ofstream medit_file(output_filename.str().c_str());
    pc3t3_save->output_to_medit(medit_file);

    //save histogram
    file_out << "done.\nSave histogram...";
    std::stringstream histo_filename;
    histo_filename << output_prefix_ << "_" << mesh_nb_ << "-histo.txt";
    save_histogram<C3T3>(histo_filename.str(), *pc3t3_save);
    file_out << "done.\n";

    int i=0;
    for ( ; lit != lines_.end() ; ++lit )
    {
      po::variables_map vm_p;
      MeshCriteria dummy_mcp = get_parameters<MeshCriteria>(*lit, vm_p);

      std::stringstream opt_out;
      opt_out << output_prefix_ << "_" << mesh_nb_ << "_";
      Optimizer<C3T3,Domain> optimizer(pc3t3_save, pdomain_builder_, i++,
                                       opt_out.str(), *lit);

      if(vm_p.count("lloyd"))
        optimizer.set_lloyd(vm_p["lloyd"].as<int>());

      if(vm_p.count("odt"))
        optimizer.set_odt(vm_p["odt"].as<int>());

      if(vm_p.count("perturb"))
        optimizer.set_perturb(vm_p["perturb"].as<double>());

      if(vm_p.count("exude"))
        optimizer.set_exude();

      std::stringstream cout_loc;
      cout_loc << "+ Queue thread. Optimize[" << filename_ << "_"
               << mesh_nb_ << "] with: " << *lit << "\n";
      std::cout << cout_loc.str();

      msg_queue_->send(optimizer);
    }
  }

  void print(const std::string& str) const
  {
    std::ofstream file_out(progout_filename().c_str(), std::ios_base::app);
    file_out << str;
  }

private:
  std::string progout_filename() const
  {
    std::size_t pos = output_prefix_.find_last_of('/');
    std::string out_pref = output_prefix_.substr(0,pos);
    std::string filename = output_prefix_.substr(pos+1);
    std::stringstream progout_filename;
    progout_filename << out_pref << "/ProgramOutput." << filename << "_" << mesh_nb_ << ".txt";

    return progout_filename.str();
  }

private:
  std::shared_ptr<Domain_builder<Domain> > pdomain_builder_;
  int mesh_nb_;
  std::string filename_;
  std::string output_prefix_;
  std::vector<std::string> lines_;
  MessageQueue* msg_queue_;
};



template <class C3T3, class Domain>
struct Optimizer
{
  Optimizer(const std::shared_ptr<C3T3>& pc3t3,
            const std::shared_ptr<Domain_builder<Domain> >& pdomain_builder,
            const int mesh_nb,
            const std::string& output,
            const std::string& command_line)
    : pc3t3_(pc3t3)
    , pdomain_builder_(pdomain_builder)
    , mesh_nb_("")
    , output_prefix_(output)
    , command_line_(command_line)
    , lloyd_(false)
    , odt_(false)
    , perturb_(false)
    , exude_(false)
  {
    std::stringstream mesh_nb_sstr;
    if ( mesh_nb < 10 )
      mesh_nb_sstr << "0";

    mesh_nb_sstr << mesh_nb;
    mesh_nb_ = mesh_nb_sstr.str();
  }


  void operator()(unsigned int nb, const std::string& str)
  {
    if ( (unsigned int)PRINT_ERROR_MSG == nb ) return print(str);

    std::size_t pos = output_prefix_.find_last_of('/');
    std::string filename = output_prefix_.substr(pos+1);

    std::stringstream cout_loc;
    cout_loc << "[" << nb << "] Launch optimizer[" << filename << mesh_nb_
             << "] with: " << command_line_ << "\n";
    std::cout << cout_loc.str();

    // Output
    std::ofstream file_out(progout_filename().c_str());

    //std::stringstream str_out;
    file_out << command_line_ << std::endl;

    CGAL::Timer timer;
    timer.start();

    file_out << "Copy C3t3...";
    std::flush(file_out);
    C3T3 c3t3 = *pc3t3_;
    file_out << "done (" << timer.time() << "s)\n";

    CGAL::Mesh_optimization_return_code code;

    if(lloyd_)
    {
      file_out << "Lloyd optimization...";
      std::flush(file_out);
      code = CGAL::lloyd_optimize_mesh_3(c3t3, pdomain_builder_->domain(), max_iteration_number=nb_lloyd_);
      file_out << "done (" << timer.time() << "s)    [" << translate_code(code) << "]\n";
    }
    timer.stop(); timer.reset(); timer.start();

    if(odt_)
    {
      file_out << "Odt optimization...";
      std::flush(file_out);
      code = CGAL::odt_optimize_mesh_3(c3t3, pdomain_builder_->domain(), max_iteration_number=nb_odt_);
      file_out << "done (" << timer.time() << "s)    [" << translate_code(code) << "]\n";
    }
    timer.stop(); timer.reset(); timer.start();

    if(perturb_)
    {
      file_out << "Perturbation...";
      std::flush(file_out);
      code = CGAL::perturb_mesh_3(c3t3, pdomain_builder_->domain(), time_limit=nb_perturb_);
      file_out << "done (" << timer.time() << "s)    [" << translate_code(code) << "]\n";
    }
    timer.stop(); timer.reset(); timer.start();

    if(exude_)
    {
      file_out << "Exudation...";
      std::flush(file_out);
      code = CGAL::exude_mesh_3(c3t3, sliver_bound=20);
      file_out << "done (" << timer.time() << "s)    [" << translate_code(code) << "]\n";
    }
    timer.stop(); timer.reset(); timer.start();

    //save mesh
    file_out << "Save mesh...";
    std::flush(file_out);
    std::stringstream output_filename;
    output_filename << output_prefix_ << mesh_nb_ << "-out.mesh";
    std::ofstream medit_file(output_filename.str().c_str());
    c3t3.output_to_medit(medit_file);

    //save histogram
    file_out << "done.\nSave histogram...";
    std::stringstream histo_filename;
    histo_filename << output_prefix_  << mesh_nb_ << "-histo.txt";
    save_histogram<C3T3>(histo_filename.str(), c3t3);
    file_out << "done.\n";
  }

  void set_lloyd(int i) { lloyd_=true; nb_lloyd_=i; }
  void set_odt(int i) { odt_=true; nb_odt_=i; }
  void set_perturb(double d) { perturb_=true; nb_perturb_=d; }
  void set_exude() { exude_=true; }

  std::string translate_code(const CGAL::Mesh_optimization_return_code& code) const
  {
    switch ( code )
    {
    case CGAL::BOUND_REACHED: return std::string("sliver bound");
    case CGAL::TIME_LIMIT_REACHED: return std::string("time limit");
    case CGAL::CANT_IMPROVE_ANYMORE: return std::string("can't improve");
    case CGAL::CONVERGENCE_REACHED: return std::string("convergence");
    case CGAL::MAX_ITERATION_NUMBER_REACHED: return std::string("max iteration number");
    }

    return std::string("UNKNOWN CODE");
  }

  void print(const std::string& str) const
  {
    std::ofstream file_out(progout_filename().c_str(), std::ios_base::app);
    file_out << str;
  }

private:
  std::string progout_filename() const
  {
    std::size_t pos = output_prefix_.find_last_of('/');
    std::string out_pref = output_prefix_.substr(0,pos);
    std::string filename = output_prefix_.substr(pos+1);
    std::stringstream progout_filename;
    progout_filename << out_pref << "/ProgramOutput." << filename << mesh_nb_ << ".txt";

    return progout_filename.str();
  }


private:
  std::shared_ptr<C3T3> pc3t3_;
  std::shared_ptr<Domain_builder<Domain> > pdomain_builder_;
  std::string mesh_nb_;
  std::string output_prefix_;
  std::string command_line_;
  bool lloyd_;
  bool odt_;
  bool perturb_;
  bool exude_;
  int nb_lloyd_;
  int nb_odt_;
  double nb_perturb_;
};



template <typename T>
T set_arg(const std::string& param_name,
          const std::string& param_string,
          const po::variables_map& vm)
{
        if(vm.count(param_name))
        {
                T param_value = vm[param_name].as<T>();
                return param_value;
        }
        else
        {
                return T();
        }
}

std::vector<std::string> split_line(const std::string& str)
{
        std::vector<std::string> args;

        std::string::size_type lastPos = str.find_first_not_of(" ", 0);
        std::string::size_type pos = str.find_first_of(" ", lastPos);
        while(pos != std::string::npos || lastPos != std::string::npos)
        {
                args.push_back(str.substr(lastPos, pos-lastPos));
                lastPos = str.find_first_not_of(" ", pos);
                pos = str.find_first_of(" ", lastPos);
        }
        return args;
}

template<typename C3t3>
void save_histogram(const std::string& filename,
                    const C3t3& c3t3)
{
        std::vector<int> histo(181,0);

  double min_value = 180.;
  double max_value = 0.;

        for (typename C3t3::Cell_iterator cit = c3t3.cells_begin() ;
       cit != c3t3.cells_end() ;
       ++cit)
        {
                if( !c3t3.is_in_complex(cit))
                        continue;

                typedef typename K::Point_3 Point_3;
                const Point_3& p0 = cit->vertex(0)->point();
                const Point_3& p1 = cit->vertex(1)->point();
                const Point_3& p2 = cit->vertex(2)->point();
                const Point_3& p3 = cit->vertex(3)->point();

                double a = CGAL::to_double(CGAL::abs(CGAL::Mesh_3::dihedral_angle(p0,p1,p2,p3)));
                histo[static_cast<int>(std::floor(a))] += 1;
    min_value = (std::min)(min_value, a);
    max_value = (std::max)(max_value, a);

                a = CGAL::to_double(CGAL::abs(CGAL::Mesh_3::dihedral_angle(p0, p2, p1, p3)));
                histo[static_cast<int>(std::floor(a))] += 1;
    min_value = (std::min)(min_value, a);
    max_value = (std::max)(max_value, a);

                a = CGAL::to_double(CGAL::abs(CGAL::Mesh_3::dihedral_angle(p0, p3, p1, p2)));
                histo[static_cast<int>(std::floor(a))] += 1;
    min_value = (std::min)(min_value, a);
    max_value = (std::max)(max_value, a);

                a = CGAL::to_double(CGAL::abs(CGAL::Mesh_3::dihedral_angle(p1, p2, p0, p3)));
                histo[static_cast<int>(std::floor(a))] += 1;
    min_value = (std::min)(min_value, a);
    max_value = (std::max)(max_value, a);

                a = CGAL::to_double(CGAL::abs(CGAL::Mesh_3::dihedral_angle(p1, p3, p0, p2)));
                histo[static_cast<int>(std::floor(a))] += 1;
    min_value = (std::min)(min_value, a);
    max_value = (std::max)(max_value, a);

                a = CGAL::to_double(CGAL::abs(CGAL::Mesh_3::dihedral_angle(p2, p3, p0, p1)));
                histo[static_cast<int>(std::floor(a))] += 1;
    min_value = (std::min)(min_value, a);
    max_value = (std::max)(max_value, a);

        }
        std::ofstream file(filename.c_str());

  file << std::setprecision(4);
  file << min_value << "\n";
  file << max_value << "\n";
        std::copy(histo.begin(), histo.end(), std::ostream_iterator<int>(file, "\n"));
}


template<typename MeshingCriteria>
MeshingCriteria get_parameters(const std::string& param_line,
                               po::variables_map& vm)
{
        po::options_description mesh("Mesh generation parameters");
        mesh.add_options()
  ("mesh", "Generate mesh")
        ("facet_angle", po::value<double>(), "Set facet angle bound")
        ("facet_size", po::value<double>(), "Set facet size bound")
        ("facet_error", po::value<double>(), "Set approximation error bound")
        ("tet_shape", po::value<double>(), "Set tet radius-edge bound")
        ("tet_size", po::value<double>(), "Set tet size bound");

        po::options_description implicit_functions("Implicit functions");
        implicit_functions.add_options()
        ("torus", "Mesh torus function")
        ("sphere", "Mesh sphere function")
        ("chair", "Mesh chair function")
        ("tanglecube", "Mesh tanglecube function")
        ("cube", "Mesh cube function")
        ("ellipsoid", "Mesh ellipsoid function")
        ("heart", "Mesh heart function")
        ("octic", "Mesh octic function");

        po::options_description optim("Optimization parameters");
        optim.add_options()
        ("exude", "Exude mesh after refinement")
        ("perturb", po::value<double>(), "Perturb mesh after refinement (sliver removal)")
        ("lloyd", po::value<int>(), "Lloyd smoothing after refinement. arg is max nb iterations")
        ("odt", po::value<int>(), "ODT smoothing after refinement. arg is max nb iterations");

        po::options_description additional_options("Options");
        additional_options.add_options()
        ("off_vertices", "Use polyhedron vertices as initialization step")
        ("no_label_rebind", "Don't rebind cell labels in medit output")
        ("show_patches", "Show surface patches in medit output");

        po::options_description cmdline_options("Usage", 1);
        cmdline_options.add(mesh).add(implicit_functions).add(optim).add(additional_options);

        std::vector<std::string> args = split_line(param_line);
        po::store(po::command_line_parser(args).options(cmdline_options).run(), vm);
        po::notify(vm);

        double facet_angle = set_arg<double>("facet_angle", "Facet angle", vm);
        double facet_size  = set_arg<double>("facet_size", "Facet size", vm);
        double facet_error = set_arg<double>("facet_error", "Facet approximation error", vm);
        double tet_shape = set_arg<double>("tet_shape","Tet shape (radius-edge)", vm);
        double tet_size = set_arg<double>("tet_size","Tet size", vm);

        typename MeshingCriteria::Facet_criteria fc(facet_angle, facet_size, facet_error);
        typename MeshingCriteria::Cell_criteria cc(tet_shape, tet_size);
        return MeshingCriteria(fc, cc);
}



template <class C3T3, class Mesh_criteria, class Domain>
void mesh(const std::string& data, const std::string& output_dir, const int nb_threads)
{
  // Create threads
  std::cout << "Launch " << nb_threads << " threads...\n";
  MessageQueue msg_queue;
  boost::thread_group threads;
  for ( int i=0; i<nb_threads; ++i)
  {
    boost::thread* exec_thread = new boost::thread(pooledThread,&msg_queue,i);
    threads.add_thread(exec_thread);
  }

  // Verify filesystem
  if(!fs::is_directory(data))
  {
    std::cout << "ERR: Problem while reading " << data << "\n";
    return;
  }

  if( output_dir == "" )
  {
    std::cout << "ERR: Wrong output directory: '" << output_dir << "'\n";
    return;
  }

  if( !fs::is_directory(output_dir) && !fs::create_directory(output_dir))
  {
    std::cout << "ERR: Problem while creating " << output_dir << "\n";
    return;
  }

  fs::path path(data);

  // Mesh
  for(fs::directory_iterator it(path), end_it ; it != end_it ; ++it)
  {
    if(fs::is_directory(*it)
       || ( fs::extension(*it) != ".off" &&
            fs::extension(*it) != ".inr" &&
            (fs::extension(fs::basename(*it)) != ".inr" || fs::extension(*it) != ".gz")) )
      continue;

    std::string line_param;
    std::string filename(fs::basename(*it));
    if ( fs::extension(*it) == ".gz" )
      filename = fs::basename(fs::basename(*it));
    std::string filename_param(data + filename + ".txt");

    std::ifstream file_param(filename_param.data()); //parameters
    if(!file_param)
    {
      std::stringstream cout_loc;
      cout_loc << "ERR: Could not find parameter file: " << filename_param << ". Skip "
               << filename << "." << std::endl;
      std::cout << cout_loc.str();
      continue;
    }

    // Timer
    CGAL::Timer timer;
    timer.start();

    //Load the domain
    std::stringstream cout_loc;
    cout_loc << "+ [" << filename << "] Create domain...";
    std::shared_ptr<Domain_builder<Domain> > pdomain_builder(new Domain_builder<Domain>(it->path().string()));
    cout_loc << "done (" << timer.time() << "s)\n";
    std::cout << cout_loc.str();

    int mesh_nb = 0;
    std::string out_prefix = output_dir + "/" + filename;
    std::getline(file_param,line_param);

    while ( ! file_param.eof() )
    {
      std::vector<std::string> lines_param;
      lines_param.push_back(line_param);

      while (   std::getline(file_param,line_param)
             && line_param.find(std::string("--mesh")) == std::string::npos )
      {
        lines_param.push_back(line_param);
      }

      // Build new mesher
      Mesher<C3T3,Mesh_criteria,Domain> mesher(pdomain_builder, mesh_nb++, filename,
                                              out_prefix, lines_param, msg_queue);

      std::stringstream cout_loc;
      cout_loc << "+ Queue thread. Mesher[" << filename << "] with: "
               << lines_param.front() << "\n";
      std::cout << cout_loc.str();

      msg_queue.send(mesher);
    }
  }

  msg_queue.deactivate();
  threads.join_all();
}





