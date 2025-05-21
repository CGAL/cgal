#ifndef ECG_BENCH_H
#define ECG_BENCH_H

#include <CGAL/Interval_nt.h>
#include <CGAL/Memory_sizer.h>

#include <string>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <iomanip>
#include <fstream>
#include "Input_data.h"
#include <variant>

class Bench
{
  int numof_f_fails;
  std::ofstream htmlout,texout;
  struct rusage before,after;
  struct timeval utime,stime;
  bool firsttime;
  char *description;
  int vert[4],hedg[4];
  int i;
  CGAL::Memory_sizer mem_sizer;
  char* Dxffile;
  bool ONLY_DXF;


  public:

  Bench( char *fhtml="benchmarks.html", char* ftex="benchmarks.tex",char* Dxffilename="",const bool only_dxf=false):
     htmlout(fhtml,std::ios::app),
     texout(ftex,std::ios::app)
  { i=0;
    numof_f_fails= 0;
    Dxffile=Dxffilename;
    ONLY_DXF = only_dxf;

    description = "Arrangement ";

    std::cout << ".:Starting Bench Suite:." <<  std::endl;
    std::cout << "bencmarking : " << description << std::endl;

    htmlout<<"Dxf file is : "<< Dxffile <<  std::endl;
    htmlout<<"<table border = \"3\">"<<std::endl;

    if(!ONLY_DXF){
    htmlout<< "<tr bgcolor=\"gray\"><td>Kernel</td><td>grid I(400 Circle_2)</td><td>grid II(800 Circle_2)</td><td>random double (100 Circle)</td><td>Dxf input. Dxf file is : "<< Dxffile <<"</td><td> rotation (125 Circle_2) </td></tr>"
    <<std::endl;

    texout << "\\documentclass[12pt]{article}"<<std::endl
    <<"\\usepackage{amsmath}"<<std::endl
    <<"\\begin{document}"<<std::endl
    <<"\\begin{tabular}{|c|c|c|c|c|c|}  \\hline $Kernel$ & grid I&grid II &random double & Dxf input. Dxf file is : "<< Dxffile <<" & rotation \\\\ \\hline"<<std::endl;
    }
    else
    {htmlout<< "<tr bgcolor=\"gray\"><td>Kernel</td><td>Dxf input. Dxf file is : "<< Dxffile <<"</td></tr>"
    <<std::endl;

    texout << "\\documentclass[12pt]{article}"<<std::endl
    <<"\\usepackage{amsmath}"<<std::endl
    <<"\\begin{document}"<<std::endl
    <<"\\begin{tabular}{|c|c|c|c|c|c|}  \\hline $Kernel$ & & Dxf input. Dxf file is : "<< Dxffile <<" \\\\ \\hline"<<std::endl;
    }
    firsttime = true;
  }

  ~Bench()
  {
    htmlout <<"</table>"<<std::endl;
    texout << "\\end{tabular}" << std::endl
           <<"\\end{document}"<<std::endl;

  }

    void open_cell()
  {
        htmlout <<"<td>";
          texout<<" ";
  }
   void close_cell()
  {
        htmlout <<"</td>";
          texout <<" ";
  }
  void open_row()
  {
        htmlout <<"<tr>";
          texout <<"";
  }
    void close_row()
  {
          htmlout <<"</tr>"<< std::endl;
          texout <<"\\\\ \\hline"<< std::endl;
  }

 void kernel(char* kernel)
 {   this->open_row();
     this->open_cell();
     htmlout << kernel;
     texout << kernel ;
     this->close_cell();
 }

  void empty_cell()
 {
   this->open_cell();
    htmlout <<  "<table border=\"1\"><tr><td>User time</td><td>Num. of fails</td><td>vertices</td><td>halfedges</td>" << std::endl
    << "<tr><td>"
    <<"-"<<"</td><td>"
    <<"-"<< "</td><td>"
    <<"-"<< "</td><td>"
    <<"-"<< "</td></tr>"
    << "</table></td>" ;
    texout<<"& -";
     this->close_cell();
 }

 void infotable(){

  texout<< "\\hline\\hline   & grid I&grid II &random double & Dxf input. Dxf file is : "<< Dxffile <<" & rotation \\\\ \\hline" << std::endl
  << "vertices" ;
  for (int k=1; k < 6 ;k++)
  {
     texout << " & "<< vert[k];
  }
   texout<<"\\\\ \\hline"<<std::endl;
   texout<<"halfedges";
     for (int l=1; l < 6 ;l++)
  {
     texout <<" & "<< hedg[l];
  }
  texout<<"\\\\ \\hline"<<std::endl;
 }

 void newDxfFilename(char* Dxffilename="")
 {Dxffile=Dxffilename;
     if(!ONLY_DXF){
            htmlout<< "<tr  bgcolor=\"gray\"><td>Kernel</td><td>grid I</td><td>grid II</td><td>random double</td><td>Dxf input. Dxf file is : "<< Dxffile <<"</td><td> rotation </td></tr>"
            <<std::endl;
    }
    else{
    htmlout<< "<tr bgcolor=\"gray\"><td>Kernel</td><td>Dxf input. Dxf file is : "<< Dxffile <<"</td> </tr>"
    <<std::endl;

    texout <<"\\begin{tabular}{|c|c|c|c|c|c|}  \\hline $Kernel$  & Dxf input. Dxf file is : "<< Dxffile <<" \\\\ \\hline"<<std::endl;

    }

 }

  void start(void)  {
  getrusage(RUSAGE_SELF,&before);

  }
  void stop(void)  {getrusage(RUSAGE_SELF,&after);}

  void summarize(int vertices,int hedges)
  {
   int temp;
   temp=CGAL::Interval_nt<>::number_of_failures();
   numof_f_fails =  temp - numof_f_fails;

   std::cout  << "  numbers_of_filter_fails : " << numof_f_fails << std::endl;
   timersub(&(after.ru_utime),&(before.ru_utime),&utime);
   timersub(&(after.ru_stime),&(before.ru_stime),&stime);

   htmlout <<  "<table border=\"1\"><tr><td>User time</td><td>Num. of fails</td><td>vertices</td><td>halfedges</td>" << std::endl
   << "<tr><td>"
   <<utime.tv_sec<<"."<< std::setw(6) << std::setfill('0')
                                       << utime.tv_usec <<"</td><td>"
   << numof_f_fails << "</td><td>"
   << vertices << "</td><td>"
   << hedges<< "</td></tr>"
   << "</table></td>";

    texout << " & " <<utime.tv_sec<<"."<< std::setw(6) << std::setfill('0')
                                       << utime.tv_usec;

    //It's necessary to make output of number of filter fails here!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    std::cout << "  vertices  : " << vertices     << std::endl
              << "  halfedges : " << hedges       << std::endl
              << "  Time (sec): " << utime.tv_sec << std::endl
              << "      (usec): " << std::setw(6) << std::setfill('0')
                                       << utime.tv_usec << std::endl;
  numof_f_fails = CGAL::Interval_nt<>::number_of_failures();
  }

  void fail(void)
  {
    htmlout <<  "<table border=\"1\"><tr><td>User time</td><td>Num. of fails</td><td>vertices</td><td>halfedges</td>" <<std::endl
    << "<tr><td>"
    <<"fail"<<"</td><td>"
    <<"fail"<< "</td><td>"
    <<"fail"<< "</td><td>"
    <<"fail"<< "</td></tr>"
    << "</table></td>" ;
    texout << "&" <<"fail";
    std::cout << " :: Abort:: "<< std::endl;
  }

 template <class Traits,class ArcContainer>
 void arrangement(const ArcContainer & ac){
 bool fail=false;
  typedef typename CGAL::Arrangement_2<Traits>                 Pmwx;
  typedef typename CGAL::Arr_naive_point_location<Pmwx>        Point_location;
  std::cout << "memory size before construction" << mem_sizer.virtual_size() << std::endl;
  std::cout << "memory resident size before insert()" << mem_sizer.resident_size () << std::endl;
  Pmwx _pm;
  Point_location _pl(_pm);

  std::cout << "Construction completed"<<std::endl;
  try{
    this->start();
    for (typename ArcContainer::const_iterator it=ac.begin();
         it != ac.end(); ++it) {
      insert(_pm,*it,_pl,std::false_type());
    };
    this->stop();
  }
  catch (...) {
    this->fail();
    fail=true;
  }
   vert[i]=(int)_pm.number_of_vertices();
   hedg[i]=(int)_pm.number_of_halfedges();


    //_pl.~Point_location();
   // _pm.~Pmwx();<Traits, class Dcel,Point_location>
   if (!fail){this->summarize(_pm.number_of_vertices(),_pm.number_of_halfedges());}

    _pm.clear();
    std::cout << "memory size after insert()" << mem_sizer.virtual_size () << std::endl;
    std::cout << "memory resident size after insert()" << mem_sizer.resident_size () << std::endl;
 }



template <class CK,class Traits,class ArcContainer>
 void grid(){

  i=1;

  ArcContainer ac;

  _bench_input_grid<CK,ArcContainer>(ac);

  this->open_cell();

     this->arrangement<Traits,ArcContainer>(ac);

  this->close_cell();
 }

 template <class CK,class Traits,class ArcContainer>
 void grid2(){
 i=2;

  ArcContainer ac;

  _bench_input_grid2<CK,ArcContainer>(ac);
  this->open_cell();
  this->arrangement<Traits,ArcContainer>(ac);
  this->close_cell();
 }
 template <class CK,class Traits,class ArcContainer>
 void random()
 {
 i=3;

     ArcContainer ac;
     if (firsttime)
     {
     _bench_input_random<CK,ArcContainer>(ac);
     std::cout << "Input from random generator!"<<std::endl;
     }
     else
      {
      _bench_input_file<CK,ArcContainer>(ac);
      std::cout << "Input from file!"<<std::endl;
      }

 this->open_cell();
     this->arrangement<Traits,ArcContainer>(ac);
 this->close_cell();
 firsttime = false;
 }
 template <class CK,class Traits,class ArcContainer>
 void dfx(char* Dfxfile){
 i=4;

      ArcContainer arc;
     _bench_input_dfx<CK,ArcContainer>(arc,Dfxfile);
 this->open_cell();
     this->arrangement<Traits,ArcContainer>(arc);

 this->close_cell();

 }

  template <class CK,class Traits,class ArcContainer>
 void rotation(){
 i=5;

      ArcContainer arc;
      try{

     _bench_input_rotation<CK,ArcContainer>(arc);
     }
     catch (...)
     {
     std::cout << "failed before input";
     }
     this->open_cell();
     this->arrangement<Traits,ArcContainer>(arc);

      this->close_cell();

     this->close_row();
 }

 template <class CK,class Traits,class ArcContainer>
 void Compute(char* dxffile)
 {
//    this->grid<CK,Traits,ArcContainer>();
//
//    this->grid2<CK,Traits,ArcContainer>();
//
//    this->random<CK,Traits,ArcContainer>();
   this->empty_cell();
    this->empty_cell();
   this->empty_cell();

 this->dfx<CK,Traits,ArcContainer>(dxffile);

//    this->rotation<CK,Traits,ArcContainer>();
   this->empty_cell();
   this->close_row();

 }

  template <class CK,class Traits,class ArcContainer>
 void Compute_no_dxf()
 {
  this->grid<CK,Traits,ArcContainer>();


  this->grid2<CK,Traits,ArcContainer>();

  this->random<CK,Traits,ArcContainer>();
   this->empty_cell();
//    this->empty_cell();
//    this->empty_cell();
//    this->empty_cell();
//    this->empty_cell();
//    this->close_row();
 this->rotation<CK,Traits,ArcContainer>();
 }

 void empty_row()
 { this->empty_cell();
   this->empty_cell();
   this->empty_cell();
   this->empty_cell();
   this->empty_cell();
   this->close_row();

 }
   template <class CK,class Traits,class ArcContainer>
 void Compute_dxf(char* dxffile)
 {
   this->dfx<CK,Traits,ArcContainer>(dxffile);
   this->close_row();
 }
};

#endif //#define ECG_BENCH_H

