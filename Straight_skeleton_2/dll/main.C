#include <iostream>
#include <vector>

#include "StraightSkeleton.h"

int 
main () 
{	
  int np;
  int *np_i;

  int numFaces=-1, numVertices=-1;
  double* xf = NULL;
  double* yf = NULL;
  int* numFace_i = NULL;
  

  std::cin >> np;
  std::cerr << np << " polygons" << std::endl;

  std::vector<double> X;
  std::vector<double> Y;

  np_i = new int[np];
  for(int i = 0; i < np; i++){
    std::cin >> np_i[i];
  std::cerr << np_i[i] << " vertices" << std::endl;
    for (int j = 0; j < np_i[i]; j++){
      double x, y;
      std::cin >> x >> y;
      X.push_back(x);
      Y.push_back(y);
    }
  }
    



  bool result = StraightSkeleton(np, np_i, &X[0], &Y[0],
				 numFaces,
				 numVertices,
				 numFace_i,
				 xf, yf,
				 true);
  if(result = false){
    std::cerr << "Straight skeleton function failed" << std::endl;
    return 0;
  }

  std::cout << "numFaces = " << numFaces << std::endl;
 
  int j =0, end=0;
  for(int i=0; i < numFaces; i++){
    std::cout << "Face " << i << ":" << std::endl;
    end += numFace_i[i];
    for(;j < end; j++){
      std::cout << "  " << xf[j] << ", " << yf[j] << std::endl;
    }
  }

  StraightSkeletonFree(numFace_i, xf, yf);
  std::cout << "done" << std::endl;
  return 0;
}
