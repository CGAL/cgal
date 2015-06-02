#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <string>
using namespace std;

int main(int argc, char* argv[])
{
  if ( !( argc == 3 ) ) {
    std::cout <<"usage: "<< argv[0] <<" [filename] " << "[num_of_segments]\n";
  } 
  ofstream ofs( (argc == 1 ) ? "data/test.cin" : argv[1] );

  //outputFile << "start writing" << endl;
  int num_of_segments = (argc == 3) ? atoi(argv[2]) : 10;
  //cout << num_of_segments << endl;
  int j = 0;
  for (int i=0 ; i < num_of_segments ; i++) {
    ofs << "s ";
    // initialize random seed
    srand ( i * time(0) );
    //generate random number: 
    int k1 = rand() % 10 + 10*(i+1) + 10*j + 1;
    srand ( k1 * time(0) );
    int k2 = rand() % 10 + 10*(i+1) + 10*j + 1;
    srand ( k2 * time(0) );
    int len = rand() % 10;
    j = j+2;
    if(i % 2 == 0) {
      //horizontal segment
      //k1 is left x, k2 is y and k1+len is right x
      ofs << k1 << " " << k2 << " " << k1+len << " " << k2 << endl;
    } else {
      //vertical segment
      //k1 is x, k2 is bottom y and k2+len is top y
      ofs << k1 << " " << k2 << " " << k1 << " " << k2+len << endl;
    }
  }

  ofs.close();
  cout << "Done!\n";

  return 0;
}
