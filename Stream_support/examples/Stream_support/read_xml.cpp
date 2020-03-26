#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <vector>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
typedef K::Point_3                                Point_3;

int main(int argc, char* argv[]){
  std::ifstream in( (argc>1)? argv[1] : "data/cloud.pol");
  boost::property_tree::ptree tree;
  boost::property_tree::read_xml(in, tree);

  std::vector<Point_3> points;

  for(boost::property_tree::ptree::value_type& node : tree.get_child("PolySet.Polygon")){
    boost::property_tree::ptree subtree = node.second;
    if( node.first == "Point" ){
      for( boost::property_tree::ptree::value_type const& v : subtree.get_child( "" ) ) {
        std::string label = v.first;

        if ( label == "<xmlattr>" ) {
          Point_3 p(subtree.get<double>( label+".X"),
                    subtree.get<double>( label+".Y"),
                    subtree.get<double>( label+".Z"));
          points.push_back(p);
        }
      }
    }
  }

  std::cout << points.size() << " points read"<< std::endl;
  return 0;
}
