#include <iostream>

// if LEDA is not installed, a message will be issued in runtime.
#ifndef CGAL_USE_LEDA
int main(int argc, char* argv[])
{

  std::cout << "Sorry, helputil needs LEDA..";
  std::cout << std::endl;

  return 0;
}

#else

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/leda_rational.h>

#include <iostream>
#include <fstream>
#include <cassert>
#include <cstdlib>
#include <ctime>
#include "read_inp.h"
#include "inter.h"


typedef PM_input< CGAL::Cartesian<double> > pminpD;
typedef pminpD::Point pntD;

typedef PM_input< CGAL::Cartesian<leda_rational> > pminpR;
typedef pminpR::Point pntR;

void print_help()
{
  std::cout << "Utilities program for CGAL Planar Map package" << std::endl;
  std::cout << "---------------------------------------------" << std::endl;
  std::cout << "The program handles input files for the Planar Map";
  std::cout << std::endl;
  std::cout << "example program dctest." << std::endl;
  std::cout << "The Input file should be at the following format:";
  std::cout << std::endl;
  std::cout << "       <n = number of points>" << std::endl;
  std::cout << "       <x1> <y1> " << std::endl;
  std::cout << "       ..." << std::endl;
  std::cout << "       <xn> <yn>" << std::endl;
  std::cout << "       <k = number of curves>" << std::endl;
  std::cout << "       <s1> <t1>" << std::endl;
  std::cout << "       ..." << std::endl;
  std::cout << "       <sk> <tk>" << std::endl;
  std::cout << "       Where xi, yi are the coordinates of point i (for ";
  std::cout << std::endl;
  std::cout << "       rational version in rational format), and si, ti ";
  std::cout << std::endl;
  std::cout << "       are the indices (in the the above list of points) ";
  std::cout << std::endl;
  std::cout << "       of the source and the target of curve i." << std::endl;
  std::cout << "" << std::endl;
  std::cout << "Usages:" << std::endl;
  std::cout << "* The first parameter is always 'e' or 'f' indicates that ";
  std::cout << std::endl;
  std::cout << "  the input file is:" << std::endl;
  std::cout << "  'e' --- exact calculations (leda_rationals)." << std::endl;
  std::cout << "  'f' --- floating point (double)." << std::endl;
  std::cout << "" << std::endl;
  std::cout << "1. Intersect Line Segments:" << std::endl;
  std::cout << "   ****  helputil <e/f> i <input_file> <output_file>";
  std::cout << std::endl;
  std::cout << "       		" << std::endl;
  std::cout << "2. Convert rationals file to doubles file (and vice versa)";
  std::cout << std::endl;
  std::cout << "   ****  helputil <e/f> c <input_file> <output_file>";
  std::cout << std::endl;
  std::cout << "" << std::endl;
  std::cout << "3. Create random (not intersected) file with <num> segments";
  std::cout << std::endl;
  std::cout << "   ****  helputil <e/f> r <num> <output_file>" << std::endl;

  std::cin.get();
  exit(1);
}

void convert(const pminpR &i1, pminpD &i2)
  // convert from rationals to doubles
{
  i2.alloc(i1.get_num_pnts(), i1.get_num_cvs());
  pntD pD;
  pntR pR;
  int i,k, j;

        
  for(i = 0; i < i1.get_num_pnts(); i++)
    {
      i1.get(i, pR);
      pD = pntD(CGAL::to_double(pR.x()), CGAL::to_double(pR.y()));
      i2.set(i, pD);
    }

  for(i = 0; i < i1.get_num_cvs(); i++)
    {
      i1.get(i, k, j);
      i2.set(i, k, j);
    }
}

void convert(const pminpD &i1, pminpR &i2)
  // convert from rationals to doubles
{
  i2.alloc(i1.get_num_pnts(), i1.get_num_cvs());
  pntD pD;
  pntR pR;
  int i,k, j;

  for(i = 0; i < i1.get_num_pnts(); i++)
    {
      i1.get(i, pD);
      pR = pntR(leda_rational(pD.x()), leda_rational(pD.y()));
      i2.set(i, pR);
    }

  for(i = 0; i < i1.get_num_cvs(); i++)
    {
      i1.get(i, k, j);
      i2.set(i, k, j);
    }
}

void intersect(const pminpD &in, pminpD &out)
{
  pntD p;
  int k, j;
  Lines<double> l(1);
  int nc, np;
  int intr_nc, intr_np;

  nc = in.get_num_cvs();
  np = in.get_num_pnts();

  Lines<double>::pmcurve *cvs = new Lines<double>::pmcurve[nc];
  Lines<double>::pmpoint *pnts = new Lines<double>::pmpoint[np];
	
  Lines<double>::pmcurve *intr_cvs;
  Lines<double>::pmpoint *intr_pnts;

  int i;
  for(i = 0; i < in.get_num_pnts(); i++)
    {
      in.get(i, p);
      pnts[i].x = p.x();
      pnts[i].y = p.y();
      pnts[i].n = i;
    }

  for(i = 0; i < in.get_num_cvs(); i++)
    {
      in.get(i, k, j);
      cvs[i].s = pnts[k];
      cvs[i].t = pnts[j];
    }

  l.Intersect(pnts, np, cvs, nc, intr_pnts, intr_np, intr_cvs, intr_nc);

  out.alloc(intr_np, intr_nc);

  for(i = 0; i < intr_np; i++)
    {
      p = pntD(intr_pnts[i].x, intr_pnts[i].y);
      out.set(i, p);
    }

  for(i = 0; i < intr_nc; i++)
    {
      out.set(i, intr_cvs[i].s.n, intr_cvs[i].t.n);
    }

  if (cvs != NULL) delete[] cvs;
  if (pnts != NULL) delete[] pnts;
  if (intr_cvs != NULL) delete[] intr_cvs;
  if (intr_pnts != NULL) delete[] intr_pnts;

}

void intersect(const pminpR &in, pminpR &out)
{
  pntR p;
  int k, j;
  Lines<leda_rational> l(0);
  int nc, np;
  int intr_nc, intr_np;

  nc = in.get_num_cvs();
  np = in.get_num_pnts();

  Lines<leda_rational>::pmcurve *cvs = new Lines<leda_rational>::pmcurve[nc];
  Lines<leda_rational>::pmpoint *pnts = new Lines<leda_rational>::pmpoint[np];
	
  Lines<leda_rational>::pmcurve *intr_cvs;
  Lines<leda_rational>::pmpoint *intr_pnts;

  int i;
  for(i = 0; i < in.get_num_pnts(); i++)
    {
      in.get(i, p);
      pnts[i].x = p.x();
      pnts[i].y = p.y();
      pnts[i].n = i;
    }

  for(i = 0; i < in.get_num_cvs(); i++)
    {
      in.get(i, k, j);
      cvs[i].s = pnts[k];
      cvs[i].t = pnts[j];
    }

  l.Intersect(pnts, np, cvs, nc, intr_pnts, intr_np, intr_cvs, intr_nc);

  out.alloc(intr_np, intr_nc);

  for(i = 0; i < intr_np; i++)
    {
      p = pntR(intr_pnts[i].x, intr_pnts[i].y);
      out.set(i, p);
    }

  for(i = 0; i < intr_nc; i++)
    {
      out.set(i, intr_cvs[i].s.n, intr_cvs[i].t.n);
    }

  if (cvs != NULL) delete[] cvs;
  if (pnts != NULL) delete[] pnts;
  if (intr_cvs != NULL) delete[] intr_cvs;
  if (intr_pnts != NULL) delete[] intr_pnts;

}


void random(int num, pminpD &out, int max = 200)
{
  srand( (unsigned)time( NULL ) );
  out.alloc(num*2, num);
  pntD p;
  int i;

  for (i =0; i < num*2; i++)
    {
      p = pntD(rand() % max, rand() % max);
      out.set(i, p);
    }

  for (i =0; i < num; i++)
    {
      out.set(i, i*2, i*2+1);
    }
}

void random(int num, pminpR &out, int max = 200)
{
  srand( (unsigned)time( NULL ) );
  out.alloc(num*2, num);
  pntR p;

  int i;
  for (i =0; i < num*2; i++)
    {
      p = pntR(leda_rational(rand() % max), 
	       leda_rational(rand() % max));
      out.set(i, p);
    }

  for (i =0; i < num; i++)
    {
      out.set(i, i*2, i*2+1);
    }
}

int main(int argc, char *argv[])
{
  if (argc < 3)
    print_help();

  bool use_rational;

  switch (argv[1][0])
    {
    case 'e':
      use_rational = true;
      break;
    case 'f':
      use_rational = false;
      break;
    default:
      print_help();
    };


  switch (argv[2][0])
    {
    case 'r': // generate random file argv[4] with argv[3] curves
      if (argc < 5)
	print_help();

      if (use_rational)
	{
	  std::cout << "Create random file using rationals with ";
	  std::cout << CGAL_CLIB_STD::atoi(argv[3]) 
		    << "segments" << std::endl;
	  std::cout << "Attention: the segments might intersect.";
	  std::cout << std::endl;
	  std::cout << "           Split them using the i parameter of "
	               "helputil before running the demo.";
	  std::cout << std::endl;
	  pminpR rnd;
	  random(CGAL_CLIB_STD::atoi(argv[3]), rnd);
	  std::ofstream f(argv[4]);
	  f << rnd;
	  f.close();
	}
      else
	{
	  std::cout << "Create random file using doubles with ";
	  std::cout << CGAL_CLIB_STD::atoi(argv[3]) 
		    << " segments" << std::endl;
	  std::cout << "Attention: the segments might intersect.";
	  std::cout << std::endl;
	  std::cout << "Split them using the i parameter of helputil "
                       "before running the demo";
	  std::cout << std::endl;
	  pminpD rnd;
	  random(CGAL_CLIB_STD::atoi(argv[3]), rnd);
	  std::ofstream f(argv[4]);
	  f << rnd;
	  f.close();
	}
      break;

    case 'i': // intersect argv[3] to argv[4]
      if (argc < 5)
	print_help();
      std::cout << "Intersect the segments in " << argv[3] 
		<< " and write the intersected segments in ";
      std::cout << argv[4] << std::endl;

      if (use_rational)
	{
	  pminpR inpR1;
	  pminpR inpR2;
			
	  std::ifstream in(argv[3]);
	  in >> inpR1;
	  in.close();
			
	  intersect(inpR1, inpR2);

	  std::ofstream out(argv[4]);
	  out << inpR2;
	  out.close();
	}
      else
	{
	  pminpD inpD1;
	  pminpD inpD2;
			
	  std::ifstream in(argv[3]);
	  in >> inpD1;
	  in.close();
			
	  intersect(inpD1, inpD2);

	  std::ofstream out(argv[4]);
	  out << inpD2;
	  out.close();
	}
      break;

    case 'c': // convert argv[3] to argv[4] switching from rationals to doubles
      if (argc < 5)
	print_help();

      std::cout << "Convert the numbers in " << argv[3] 
		<< " and write the new file " << argv[4] << std::endl;

      if (use_rational)
	{
	  pminpR inpR;
	  pminpD inpD;
			
	  std::ifstream in(argv[3]);
	  in >> inpR;
	  in.close();
			
	  convert(inpR, inpD);

	  std::ofstream out(argv[4]);
	  out << inpD;
	  out.close();
	}
      else
	{
	  pminpR inpR;
	  pminpD inpD;
			
	  std::ifstream in(argv[3]);
	  in >> inpD;
	  in.close();
			
	  convert(inpD, inpR);

	  std::ofstream out(argv[4]);
	  out << inpR;
	  out.close();
	}

      break;

    default:
      print_help();
      break;
    };
	

  return 0;
}

#endif // CGAL_USE_LEDA
