#include <fstream>
#include "dxf_converter.h"

int main(int argc, char* argv[]){

Dxf_converter converter;
char* Benchfilename;
char exten[4];


if (argv[1] != NULL)
 {
                int len =strlen(argv[1]);
                for (int j=0; j < 3 ; j++)
                {
                  exten[j]=argv[1][len - 3 + j];
                }
                if (strncmp(exten,"dxf",3))
                {
                  std::cout<< "File is not correct (*.dxf is needed)." << std::endl;
                 return 0;
                }
                else{
                std::cout<< "File "<< argv[1] << " is correct."<<std::endl;
                }
                if (argc >2 and argv[2] != NULL)
                 {
                        int len =strlen(argv[2]);
                        for (int j=0; j < 4 ; j++)
                        {
                                  exten[j]=argv[2][len - 4 + j];
                        }
                        if (strncmp(exten,"arr",4) !=0)
                        {
                                 std::cout<< "File "<< argv[2] << " is not correct (*.arr is needed)." << std::endl;
                                 return 0;
                        }
                        else{
                        std::cout<< "File "<< argv[2] << " is correct." <<std::endl;
                        strcpy(Benchfilename,argv[2]);
                        }
                 }
                else
                 {strcpy(Benchfilename,"");
                int len= strlen(argv[1]);
                strncat(Benchfilename,argv[1],len-4);
                strcat(Benchfilename,".arr");
                 }

 }
else
 {
 cout<<"file "<< argv[1]<< " is not found"<<std::endl;
 return 0;
 }
  std::cout<<Benchfilename<<std::endl;

 std::ifstream fintest;
 fintest.open(argv[1]);
 if (!fintest.is_open())
  {
    cout<<"file "<<argv[1] << " is not found"<<std::endl;
fintest.close();
    return 0;
  }

 std::ifstream fin(argv[1]);
 std::ofstream fout(Benchfilename);
// std::ofstream fout("test.arr");
 fout<< "FileFormat(\"AcsBenchmark\",0,1)"<<std::endl;
 fout<<"BenchmarkName(\""<<Benchfilename<<"\")"<<std::endl;
 fout<<"Classification(\"Arrangement\",\" Line\",\" BoundedArcs\",\"Rnd\",\"1\",\"Jan-2006\")"<<std::endl;
 fout<<"Classification(\"Arrangement\",\" Circles\",\" FullCurves\",\"Rnd\",\"1\",\"Jan-2006\")"<<std::endl;
 fout<<"Classification(\"Arrangement\",\" Conics\",\" FullCurves\",\"Rnd\",\"1\",\"Jan-2006\")"<<std::endl;
 fout<<"Classification(\"Arrangement\",\" Conics\",\" BoundedArcs\",\"Rnd\",\"1\",\"Jan-2006\")"<<std::endl;
 fout <<std::endl;
 converter(fin,fout);
 std::cout<<"convertation-finished"<<std::endl;
 fin.close();
 fout.close();
 std::cout<<"files closing-finished"<<std::endl;
return 1;
}
