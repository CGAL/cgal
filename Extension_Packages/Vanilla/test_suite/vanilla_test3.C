#include <CGAL/basic.h>
#include <iostream>
#include <fstream>
#include <CEP/Vanilla/Flavored_object.h>



int main(int argc, char** argv)
{
   if (argc != 2)
   {
      std::cerr << "Usage: " << argv[0] << " <answer_file>" << std::endl;
      exit(1);
   }

   std::ifstream answers(argv[1]);

   if (!answers)
   {
      std::cerr << argv[1] << ": no such file or directory" << std::endl;
      exit(1);
   }

   Flavor orig_f, f;
   Flavor right_f;
   std::cin >> orig_f;
   f = orig_f;
   do
   {
      answers >> right_f;
      if (!answers || f != right_f)
      {
         std::cerr << f << " != " << right_f << std::endl;
         exit(1);
      }
      f = flavor_enhance(f);
   }
   while (f != orig_f);
   exit (0);
}
