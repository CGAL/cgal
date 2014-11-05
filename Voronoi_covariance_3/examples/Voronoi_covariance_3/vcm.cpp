#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Voronoi_covariance_3/voronoi_covariance_3.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Fuzzy_sphere.h>

#include <boost/timer.hpp>
#include <boost/progress.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>

#include <iostream>
#include <fstream>
#include <vector>

    template <class T>
bool load(const std::string &s,
          std::vector<T> &data)
{
    std::ifstream is (s.c_str());

    data.clear();
    std::copy (std::istream_iterator<T> (is),
               std::istream_iterator<T> (),
               std::back_inserter(data));

    return true;
}

    template <class T>
bool load_if_newer (const std::string &t,
                    const std::string &s,
                    std::vector<T> &data)
{
    if (!boost::filesystem::exists(s))
        return false;

    if (boost::filesystem::last_write_time(s) <
        boost::filesystem::last_write_time(t))
        return false;

    return load (s, data);
}

    template <class T>
bool save(const std::string &s,
          const std::vector<T> &data)
{
    std::ofstream os (s.c_str());

    std::copy (data.begin(), data.end(),
               std::ostream_iterator<T> (os));

    return true;
}

    template <class K, class Covariance>
void convolve (const std::vector< typename CGAL::Point_3<K> > &points,
               const std::vector<Covariance> &cov,
               std::vector<Covariance> &ncov,
               double r)
{
    typedef typename CGAL::Point_3<K> Point;
    typedef typename CGAL::Search_traits_3<K> Traits;
    typedef typename CGAL::Kd_tree<Traits> Tree;
    typedef typename CGAL::Fuzzy_sphere<Traits> Fuzzy_sphere;

    Tree tree (points.begin(), points.end());

    boost::progress_display p(points.size(), std::cerr);

    std::map<Point, size_t> indices;
    for (size_t i = 0; i < points.size(); ++i)
        indices[points[i]] = i;

    ncov.clear();
    for (size_t i = 0; i < points.size(); ++i)
    {
        std::vector<Point> nn;
        tree.search(std::back_inserter(nn),
                    Fuzzy_sphere (points[i], r));

        Covariance m;
        for (size_t k = 0; k < nn.size(); ++k)
            m += cov[indices [nn[k]]];
        ncov.push_back(m);

        ++p;
    }
}

    template <class K, class Covariance>
void offset (const std::vector< typename CGAL::Point_3<K> > &points,
             std::vector<Covariance> &cov,
             double R, size_t N)
{
    typedef CGAL::Delaunay_triangulation_3<K> DT;

    typename CGAL::Voronoi_covariance_3::Sphere_discretization<K> sphere(R, N);
    DT dt(points.begin(), points.end());

    boost::progress_display p(points.size(), std::cerr);
    cov.clear();

    boost::timer t;
    for (size_t i = 0; i < points.size(); ++i)
    {
        typename DT::Vertex_handle vh = dt.nearest_vertex(points[i]);
        Covariance c = CGAL::Voronoi_covariance_3::voronoi_covariance_3(dt, vh, sphere);
        cov.push_back(c);
        ++p;
    }
    std::cerr << "time: " <<  t.elapsed() << "s\n";
}

namespace po = boost::program_options;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::FT FT;
typedef K::Point_3 Point;
typedef CGAL::Voronoi_covariance_3::Voronoi_covariance_3<FT> Covariance;

int main(int argc, char **argv)
{
    double R, r;
    size_t N = 20;

    po::options_description desc("Options");
    desc.add_options()
        ("help", "produce help message")
        ("N", po::value<size_t>(&N)->default_value(20),
         "number of planes used to discretize sphere:\n"
         "    20 for dodecahedron\n"
         "    12 for icosahedron")
        ("R", po::value<double>(&R)->default_value(0.1),
         "offset radius")
        ("r", po::value<double>(&r)->default_value(0),
         "convolution radius")
        ("input", po::value< std::vector<std::string> >(),
         "input cloud/xyz files");
    po::positional_options_description p;
    p.add("input", -1);

    po::variables_map args;
    po::store(po::command_line_parser(argc, argv).options(desc)
              .style (po::command_line_style::default_style |
                      po::command_line_style::allow_long_disguise)
              .positional(p).run(), args);
    po::notify(args);

    if (args.count("help"))
    {
        std::cout << desc;
        return 1;
    }

    if (!args.count ("input"))
    {
        std::cout << argv[0] << ": no input file\n";
        return 1;
    }

    std::vector< std::string > files =
        args["input"].as< std::vector<std::string> >();

    for (size_t i = 0; i < files.size(); ++i)
    {
        std::vector<Covariance> cov, ccov;
        std::vector<Point> points;

        std::string source = files[i],
            bn = boost::filesystem::basename(source),
            fncov = (boost::format("%s-R=%g.p") % bn % R).str(),
            fnccov = (boost::format("%s-R=%g-r=%g.p") % bn % R % r).str();

        load(files[i], points);

        //std::cerr << "built Delaunay triangulation in " << t.elapsed() << "s\n";
        //std::cerr<< "number of vertices: " << dt.number_of_vertices() << "\n";

        std::cerr << "offsetting " << source << " into " << fncov << "\n";

        if (!load_if_newer (source, fncov, cov) ||
            cov.size() != points.size())
        {
            offset (points, cov, R, N);
            save (fncov, cov);
        }

        if (r == 0)
            continue;

        std::cerr << "convolving " << source << " into " << fnccov << "\n";

        if (!load_if_newer (fncov, fnccov, ccov) ||
            ccov.size() != points.size())
        {
            convolve (points, cov, ccov, r);
            save (fnccov, ccov);
        }
    }

    return 0;
}

