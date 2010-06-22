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
#include <boost/variant.hpp>
#include <exception>



class Bench
{
private:

  int numof_f_fails;
  std::ofstream htmlout,texout;
  struct rusage before,after;
  struct timeval utime,stime;
  bool firsttime;
  int vert[4],hedg[4];
  int i;
typedef std::size_t size_type;
  CGAL::Memory_sizer mem_sizer;
  std::string Dxffile;
  bool ONLY_DXF;
  double MemBefore,MemAfter;
  
public:
  
  Bench(const bool only_dxf=false, const char *fhtml="benchmarks.html", const char* ftex="benchmarks.tex",const char* Dxffilename=""):
     htmlout(fhtml,std::ios::app),
     texout(ftex,std::ios::app){ 
	i=0;
    	numof_f_fails= 0;
    	Dxffile=Dxffilename;
    	ONLY_DXF = only_dxf;
    
    	std::cout << ".:Starting Bench Suite:." <<  std::endl;
    
    	htmlout<<"Dxf file is : "<< Dxffile <<  std::endl;

    	htmlout<<"<table border = \"3\">"<<std::endl;
    
    	if(!ONLY_DXF){
    		htmlout<< "<tr bgcolor=\"gray\"><td>Kernel</td><td>grid I(400 Circle_2)</td><td>grid II(800 Circle_2)</td><td>random double (100 Circle)</td><td> rotation (125 Circle_2) </td></tr>" 
    		<<std::endl;
    
    		texout << "\\documentclass[12pt]{article}"<<std::endl
    		<<"\\usepackage{amsmath}"<<std::endl
		<<"\\usepackage{rotating}"<<std::endl
    		<<"\\begin{document}"<<std::endl
		<<"\\newsavebox{\\foo}"<<std::endl
		<<"\\savebox{\\foo}"<<std::endl
		<<"{\\parbox{5cm}{"<<std::endl
    		<<"\\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|}  \\hline Kernel &\\multicolumn{4}{c|}{ grid I}&\\multicolumn{4}{c|}{grid II }&\\multicolumn{4}{c|}{random double} &\\multicolumn{4}{c|}{ rotation }\\\\ \\hline"<<std::endl;
    		}
    	
	else	{
		htmlout<< "<tr bgcolor=\"gray\"><td>Kernel</td><td>Dxf input is: "<< Dxffile <<"</td></tr>" 
    		<<std::endl;
    
    		texout << "\\documentclass[12pt]{article}"<<std::endl
    		<<"\\usepackage{amsmath}"<<std::endl
    		<<"\\begin{document}"<<std::endl
    		<<"\\begin{tabular}{|c|c|c|}  \\hline $Kernel$ & & \\multicolumn{4}{c|}{Dxf input is:"<< Dxffile <<"} \\\\ \\hline"<<std::endl;
    		}
    	firsttime = true;
  
  }

 void infotable(){

    	if(!ONLY_DXF){
  		texout<< "\\hline\\hline  & \\multicolumn{4}{c|}{ grid I}&\\multicolumn{4}{c|}{grid II }&\\multicolumn{4}{c|}{random double} &\\multicolumn{4}{c|}{ rotation }\\\\ \\hline" << std::endl
  		<< "vertices" ;
  		for (int k=1; k < 6 ;k++) { 
     			texout << " &\\multicolumn{4}{c|}{"<< vert[k]<<"}";
  		}
		texout<<"\\\\ \\hline"<<std::endl;
   		texout<<"halfedges";

     		for (int l=1; l < 6 ;l++){ 
     			texout << " &\\multicolumn{4}{c|}{"<< hedg[l]<<"}";
  		}
  		texout<<"\\\\ \\hline"<<std::endl;
    		}
    	
	else	{
		 texout<< "\\hline\\hline    &\\multicolumn{4}{c|}{Dxf input is:"<< Dxffile <<"} \\\\ \\hline" << std::endl
  		<< "vertices" ; 
     		texout <<  " &\\multicolumn{4}{c|}{"<< vert[4]<<"}";
		texout<<"\\\\ \\hline"<<std::endl;
   		texout<<"halfedges";
     		texout <<" &\\multicolumn{4}{c|}{"<< hedg[4]<<"}";
  		texout<<"\\\\ \\hline"<<std::endl;
    		}
 }

  ~Bench(){
	std::cout<<"destructor"<<std::endl;
    	htmlout <<"</table>"<<std::endl;
	infotable();

    	texout << "\\end{tabular}" << std::endl
	<<"}"<< std::endl
	<<"}"<< std::endl
	<<"\\begin{turn}{-90}\\usebox{\\foo}\\end{turn}"<< std::endl
        <<"\\end{document}"<<std::endl;
  }
private:
  void open_cell(){	
        htmlout <<"<td>";
  	texout<<" ";	
  }

  void close_cell(){	
        htmlout <<"</td>";
  	texout <<" ";
  }
  void open_row(){	
        htmlout <<"<tr>";
  	texout <<"";
  }

  void close_row(){
  	htmlout <<"</tr>"<< std::endl;
  	texout <<"\\\\ \\hline"<< std::endl;
  }
  

 void empty_cell(){
   	this->open_cell();
    	htmlout <<  "<table border=\"1\"><tr><td>User time</td><td>Num. of fails</td><td>vertices</td><td>halfedges</td><td>MemBefore</td><td>MemAfter</td>" << std::endl
    	<< "<tr><td>"
    	<<"-"<<"</td><td>"
    	<<"-"<< "</td><td>"
    	<<"-"<< "</td><td>" 
    	<<"-"<< "</td><td>"
    	<<"-"<< "</td><td>" 
    	<<"-"<< "</td></tr>"
    	<< "</table></td>" ; 
    	texout<<"& -& -&-& -";
     	this->close_cell();
 } 

  void start(void){
  	getrusage(RUSAGE_SELF,&before);
  }

  void stop(void){
	getrusage(RUSAGE_SELF,&after);
  }

  void summarize(int vertices,int hedges){
   	int temp; 
//   	temp=CGAL::Interval_nt<>::number_of_failures();	
   	numof_f_fails =  temp - numof_f_fails; 
   	std::cout  << "  numbers_of_filter_fails : " << numof_f_fails << std::endl;
   	
	timersub(&(after.ru_utime),&(before.ru_utime),&utime);
   	timersub(&(after.ru_stime),&(before.ru_stime),&stime);
   	
	htmlout <<  "<table border=\"1\"><tr><td>User time</td><td>Num. of fails</td><td>vertices</td><td>halfedges</td><td>MemBefore</td><td>MemAfter</td>" << std::endl
   	<< "<tr><td>"
   	<<utime.tv_sec<<"."<< std::setw(6) << std::setfill('0')<< utime.tv_usec <<"</td><td>"
   	<< numof_f_fails << "</td><td>"
   	<< vertices << "</td><td>" 
   	<< hedges<<"</td><td>"
  	<< MemBefore << "</td><td>" 
   	<< MemAfter<< "</td></tr>"
   	<< "</table></td>";
    	
	texout << " & " <<utime.tv_sec<<"."<< std::setw(6) << std::setfill('0')
	     			  << utime.tv_usec;
    	texout << "&"<< numof_f_fails;
	texout << "&"<< MemBefore;
	texout << "&"<<MemAfter;
    	std::cout << "  vertices  : " << vertices     << std::endl
	      << "  halfedges : " << hedges       << std::endl
	      << "  Time (sec): " << utime.tv_sec << std::endl
	      << "      (usec): " << std::setw(6) << std::setfill('0')
	     			  << utime.tv_usec << std::endl;
//  	numof_f_fails = CGAL::Interval_nt<>::number_of_failures();
  }
  
  void fail(void){
    	htmlout <<  "<table border=\"1\"><tr><td>User time</td><td>Num. of fails</td><td>vertices</td><td>halfedges</td><td>MemBefore</td><td>MemAfter</td>"<<std::endl
    	<< "<tr><td>"
    	<<"fail"<<"</td><td>"
    	<<"fail"<< "</td><td>"
    	<<"fail"<< "</td><td>" 
    	<<"fail"<< "</td><td>"
    	<<"fail"<< "</td><td>" 
    	<<"fail"<< "</td></tr>"
    	<< "</table></td>" ;
    	texout << "& fail & fail & fail & fail";
    	std::cout << " :: Abort:: "<< std::endl;
  }
  

 template <class Traits,class ArcContainer> 
 void arrangement(const ArcContainer & ac){
 	bool fail=false;
  	typedef typename CGAL::Arrangement_2<Traits>                 Pmwx;
  	typedef typename CGAL::Arr_naive_point_location<Pmwx>        Point_location;
  	
	std::cout << "memory size before construction" << mem_sizer.virtual_size() << std::endl;
  	std::cout << "memory resident size before insert()" << mem_sizer.resident_size () << std::endl;
  	MemBefore = mem_sizer.virtual_size ()/1000000;
	Pmwx _pm;
  	Point_location _pl(_pm);
	try{
    		this->start(); 
      		insert(_pm,ac.begin(),ac.end(),boost::false_type());
    		this->stop();
  	} catch (std::exception &e) {
    		this->fail();
                std::cout << "Exception:" << std::endl << e.what() << std::endl;
    		fail=true;
  	} catch(...) {
                this->fail();
    		fail=true; 
        }
   	vert[i]=(int)_pm.number_of_vertices();
   	hedg[i]=(int)_pm.number_of_halfedges();

   	if (!fail){this->summarize(_pm.number_of_vertices(),_pm.number_of_halfedges());}
   
    	
        MemAfter = mem_sizer.virtual_size ()/1000000;
	std::cout << "memory size after insert()" << mem_sizer.virtual_size () << std::endl; 
    	std::cout << "memory resident size after insert()" << mem_sizer.resident_size () << std::endl;
        _pm.clear();
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
 void random(){ 
 	i=3; 
     	ArcContainer ac;
     	if (firsttime){
     		_bench_input_random<CK,ArcContainer>(ac);
     		std::cout << "Input from random generator!"<<std::endl;
     	}
     	else{
      		_bench_input_file<CK,ArcContainer>(ac);
      		std::cout << "Input from file!"<<std::endl;
      	}
	this->open_cell();
     	this->arrangement<Traits,ArcContainer>(ac);  
 	this->close_cell();
 	firsttime = false;
 }

 template <class CK,class Traits,class ArcContainer>
 void dfx(const char* Dfxfile){
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
      	bool fail=false;
	ArcContainer arc; 
      	try{
     		_bench_input_rotation<CK,ArcContainer>(arc);
     	}
     	catch (...)
     	{	
		this->fail();
    		fail=true;
     		std::cout << "failed before input";
     	}
   	
	if (!fail){
     		this->open_cell();
     		this->arrangement<Traits,ArcContainer>(arc); 
      		this->close_cell();
     		this->close_row();
	}
 }
public:
 template <class CK,class Traits,class ArcContainer>
 void Compute(const char* dxffile){
     	if(!ONLY_DXF){
	this->grid<CK,Traits,ArcContainer>();
   	this->grid2<CK,Traits,ArcContainer>();
   	this->random<CK,Traits,ArcContainer>();
 	this->dfx<CK,Traits,ArcContainer>(dxffile);
	this->rotation<CK,Traits,ArcContainer>();
    	}
    	else{   
   	this->dfx<CK,Traits,ArcContainer>(dxffile);
   	this->close_row();
    	}
 }
 
  template <class CK,class Traits,class ArcContainer>
 void Compute_no_dxf(){
     	if(!ONLY_DXF){
  	this->grid<CK,Traits,ArcContainer>();
  	this->grid2<CK,Traits,ArcContainer>();
  	this->random<CK,Traits,ArcContainer>();
	//this->empty_cell();
	this->rotation<CK,Traits,ArcContainer>();
    	}
 }
 
 void empty_row(){ 
	this->empty_cell();
   	this->empty_cell();
   	this->empty_cell();
   	this->empty_cell();
   	this->empty_cell();
   	this->close_row();
 }

template <class CK,class Traits,class ArcContainer>
 void Compute_dxf(const char* dxffile){
   	this->dfx<CK,Traits,ArcContainer>(dxffile);
   	this->close_row(); 
 }

 void kernel(const char* kernel){ 
	this->open_row();
     	this->open_cell();
     	htmlout << kernel;
     	texout << kernel ;
     	this->close_cell();
 }

 
 void newDxfFilename(const char* Dxffilename=""){
	const char* newDxf = Dxffilename;
	
     	if(!ONLY_DXF){
    		htmlout<< "<tr  bgcolor=\"gray\"><td>Kernel</td><td>grid I</td><td>grid II</td><td>random double</td><td>Dxf input is: "<< newDxf <<"</td><td> rotation </td></tr>" 
    		<<std::endl; 
		infotable();
    		texout <<"\\end{tabular}"<<std::endl<<"\\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|}  \\hline Kernel &\\multicolumn{4}{c|}{ grid I}&\\multicolumn{4}{c|}{grid II }&\\multicolumn{4}{c|}{random double} & \\multicolumn{4}{c|}{Dxf input is:"<< newDxf<<" }&\\multicolumn{4}{c|}{ rotation }\\\\ \\hline"<<std::endl;
    	}
    	else{   
		infotable();
    		htmlout<< "<tr bgcolor=\"gray\"><td>Kernel</td><td>Dxf input is: "<< newDxf <<"</td> </tr>"
    		<<std::endl;
    		texout <<"\\end{tabular}"<<std::endl<<"\\begin{tabular}{|c|c|c|c|c|}  \\hline $Kernel$ &\\multicolumn{4}{c|}{Dxf input is:"<< newDxf <<"} \\\\ \\hline"<<std::endl; 
    	}
	Dxffile=Dxffilename;
 }
};

#endif //#define ECG_BENCH_H

