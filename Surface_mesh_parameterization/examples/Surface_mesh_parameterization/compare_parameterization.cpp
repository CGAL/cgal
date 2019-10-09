#include <cstdlib>
#include <iostream>
#include <fstream>

int main(int argc, char** argv)
{
  std::string inputMesh = (argc>1) ? argv[1] : "data/mannequin-devil.off";
  std::string bash;
  int out;

  bash = "./square_border_authalic_parameterizer " + inputMesh;
  out = system(bash.c_str());
  if (out == 0)
      std::cout << "square_border_authalic_parameterizer output written in discrete_result.off\n";

  bash = "./square_border_iterative_parameterizer " + inputMesh;
  out = system(bash.c_str());
  if (out == 0)
    std::cout << "square_border_iterative_parameterizer output written in discrete_result.off\n";

}
