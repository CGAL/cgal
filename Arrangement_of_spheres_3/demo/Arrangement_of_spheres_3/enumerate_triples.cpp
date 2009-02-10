#include <iostream>
#include <limits>

int main(int , char *[]) {

  int max[4]={0,0,0,0};
  for (int i=0; i< std::numeric_limits<int>::max(); ++i) {
    int incr=i%4;
    ++max[incr];
    for (int j=0; j < max[(incr+1)%4]; ++j){
      for (int k=0; k < max[(incr+2)%4]; ++k){
	for (int l=0; l < max[(incr+3)%4]; ++l){
	  int vals[4];
	  vals[incr]=max[incr];
	  vals[(incr+1)%4]=j;
	  vals[(incr+2)%4]=k;
	  vals[(incr+3)%4]=l;
	  if (vals[0]*vals[0]+vals[1]*vals[1]+vals[2]*vals[2]==vals[3]*vals[3]){
	    std::cout << vals[0] << " " << vals[1] << " " << vals[2] << " " 
		      << vals[3] << std::endl;
	  }
	}
      }
    }
  }

  return EXIT_SUCCESS;
}
