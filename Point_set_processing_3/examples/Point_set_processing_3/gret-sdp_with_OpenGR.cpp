#include <CGAL/Simple_cartesian.h>
#include <CGAL/property_map.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/Cartesian_matrix.h>
#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/OpenGR/gret_sdp.h>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/filesystem.hpp>

#include <fstream>
#include <iostream>
#include <utility>


using namespace std;
namespace pt = boost::property_tree;
namespace fs = boost::filesystem;

typedef double Scalar;
typedef CGAL::Simple_cartesian<Scalar> K;
typedef K::Point_3 Point_3;
typedef K::Vector_3 Vector_3;
typedef boost::tuple<Point_3, int, Vector_3> Pwn;
typedef CGAL::Nth_of_tuple_property_map<0, Pwn> Point_map;
typedef CGAL::Nth_of_tuple_property_map<1, Pwn> Index_map;
typedef CGAL::Nth_of_tuple_property_map<2, Pwn> Normal_map;
typedef std::vector<Pwn> Patch;

enum {Dim = 3};
typedef Eigen::Matrix<K::FT, Dim+1, Dim+1> MatrixType;

namespace params = CGAL::parameters;

struct RegistrationProblem {
    int n;
    int m;
    int d;
    vector<Patch> patches;
};

template <typename TrRange>
void extractPatchesAndTrFromConfigFile(const string& configFilePath,  RegistrationProblem& problem, TrRange& transformations);

int main(int argc, const char** argv)
{
    RegistrationProblem problem;
    vector<MatrixType> gt_transformations;

    extractPatchesAndTrFromConfigFile("/home/felix/Workspace/GSoC20/cgal-felix/build/examples/Point_set_processing_3/gret-sdp-data/config.json", problem, gt_transformations);

    const int d = problem.d;
    const int n = problem.n;
    const int m = problem.m;
    const vector<Patch>& patches = problem.patches;

    CGAL::OpenGR::GRET_SDP matcher;
    matcher.registerPatches(patches, problem.n, params::point_map(Point_map())
                                                .normal_map(Normal_map())
                                                .vertex_index_map(Index_map()));
}

template <typename PatchRange>
void readPatch(const string& file, PatchRange& patch){
    int num_points;
    fs::ifstream file_stream(file);
    file_stream >> num_points;
    patch.reserve(num_points);

    K::FT x, y, z;
    int index;
    while(file_stream >> index){
        file_stream >> x;
        file_stream >> y;
        file_stream >> z;        
        patch.emplace_back(Point_3(x,y,z), index, Vector_3());
  }
}

void readTransformation(const string& filePath, MatrixType& trafo){
    int rows, cols;
    fs::ifstream file(filePath);
    if(file.is_open()){
        file >> rows >> cols;
        if(cols != Dim+1 || rows != Dim+1)
            throw std::runtime_error("matrices have to be of size " + to_string(Dim+1) + "x" + to_string(Dim+1));
        for(int i = 0; i < cols; i++)
            for (int j = 0; j < rows; j++)
                file >> trafo(i, j);
    }
}

template <typename TrRange>
void extractPatchesAndTrFromConfigFile(const string& configFilePath,  RegistrationProblem& problem, TrRange& transformations){
    const string workingDir = fs::path(configFilePath).parent_path().native();

    pt::ptree root;
    pt::read_json(configFilePath, root);

    int n = root.get<int>("n");
    int m = root.get<int>("m");
    int d = root.get<int>("d");

    problem.n = n;
    problem.m = m;
    problem.d = d;

    vector< string  > patchFiles;

    for (pt::ptree::value_type &item : root.get_child("patches"))
    {
        patchFiles.push_back(item.second.data());
    }

    if(patchFiles.size() != m)
        throw runtime_error("Number of patches m and number of given patch files is not the same.");

    if(d != Dim)
        throw runtime_error("Dimension of point type has to be " + to_string(Dim));

    // read patch files
    problem.patches.resize(m);
    ifstream patch_file;
    for(int i = 0; i < m; i++){
        readPatch(workingDir + "/" + patchFiles[i], problem.patches[i]);
    }
    
    vector< string  > transformationFiles;
    for (pt::ptree::value_type &item : root.get_child("gt_trafos"))
        transformationFiles.push_back(item.second.data());

    if(transformationFiles.size() != m)
        throw runtime_error("Number of transformations and number of given transformation files is not the same.");

    transformations.reserve(m);
    MatrixType trafo;
    for(int i = 0; i < m; i++){
        readTransformation(workingDir + "/" + transformationFiles[i], trafo);
        transformations.emplace_back(trafo);
    }

}

