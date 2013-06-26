#ifndef CHECK_SOLVER_H
#define CHECK_SOLVER_H
#include <CGAL/Polynomial/internal/polynomial_generators.h>
#include <CGAL/Timer.h>
#include <vector>
#include <iterator>

struct True{};
struct False{};



template <class K, class CE=False, bool CA=true>
class Check_solver
{
public:
  typedef typename K::Root_stack Root_stack;
  typedef typename K::Root_container Root_container;
  typedef std::vector<double> DV;
  typedef typename K::Root Rt;
  typedef typename K::Function Fn;
  typedef typename Fn::NT NT;
  typedef K        Kernel;

  K k_;
  bool verbose;
  typename K::Construct_function cf_;

  double total_time;
  unsigned int num_roots;
  double ce_time;
  unsigned int nce;

  Check_solver(const K& k, bool verbose): k_(k), verbose(verbose), cf_(k_.construct_function_object()){}

  /*void time_estimate(const Fn &, const Rt& , False){
  }


  void time_estimate(const Fn &q, const Rt& start, True){
    Root_stack s= k_.root_stack_object(q, start, std::numeric_limits<Rt>::infinity());
    CGAL::Timer timer;
    timer.start();
    do {
      s.compute_estimate(start);
      ++nce;
    } while (timer.time() <3 && !verbose);
    timer.stop();
    ce_time= timer.time();
    }*/

  template <class It>
  void check_polynomial(const Fn &q, It roots_b, It roots_e,
			const Rt& start = Rt(-1000.0),
			const Rt& end = Rt(1000.0) ) {

    CGAL::Timer timer;

    //typename K::Is_even_multiplicity iem = k_.is_even_multiplicity_object(q);
    //int current_root=0;
    if (verbose) {
      std::cout << "Polynomial: " << q << std::endl;
      std::cout << "Window: [";

      std::cout<< CGAL::to_double(start);
      std::cout << ", ";
      std::cout << CGAL::to_double(end);
      std::cout << ")" << std::endl;
      std::cout << "Solution: ";
      for (It c= roots_b; c!= roots_e; ++c) {
	/*if (CGAL::to_double(start) >= a[i]-.00001){
	  ++current_root;
	  } else if (CGAL::to_double(end) <= a[i]){

	  } else {*/
	std::cout << "<" << *c  << "> ";
	//}
      }
      std::cout << std::endl;

      std::cout << "Solver: ";
    }
    std::vector<Rt> roots;
    DV comp;
    int reps=0;
    timer.start();
    int total_roots=0;
    do {
      Root_stack s= k_.root_stack_object(q, start, end);
      Rt last_root= -Rt(1000.0);
      while (!s.empty()) {
	if (reps==0  && verbose) {
	  //comp.push_back(CGAL::to_double(r));
	  if (last_root > s.top()) {
	    std::cerr << "NUMERIC ISSUE last root was " << last_root << " and current root is "
		      << s.top() << std::endl;
	  }
	  assert(last_root<= s.top());
	  last_root=s.top();
	  //Rt cur= s.top();
	  //if (cur != std::numeric_limits<Rt>::infinity()) {
	  /*if (cur == std::numeric_limits<Rt>::infinity() || cur > 10000){
	    int i=0;
	    }*/
	  roots.push_back(s.top());

	  std::cout << "<" << CGAL::to_double(s.top());
	  //if (iem(s.top())) { std::cout <<"E";}
	  std::cout << "> " << std::flush;
	  /*}
	    else {

	    }*/
	}
	s.pop();
	++total_roots;
      }
      ++total_roots;
      if (reps==0 && verbose && CA) {
	std::vector<bool> taken_maple(roots_e-roots_b, false), taken_solver(roots.size(), false);
	for (unsigned int i=0; i< roots.size(); ++i) {
	  double rd= CGAL::to_double(roots[i]);
	  //bool ie= iem(roots[i]);
	  int mult=1;
	  //if (ie) mult=2;
	  for (int k=0; k<mult; ++k) {
	    int mm=-1;
	    double md=.5;
	    for (It c= roots_b; c != roots_e; ++c) {
	      if (std::abs(rd-*c) < md  && taken_maple[c-roots_b]==false) {
		md= std::abs(rd-*c);
		mm=c-roots_b;
	      }
	    }
	    if (mm != -1) {
	      taken_maple[mm]=true;
	      taken_solver[i]=true;
	    }
	  }
	}

	bool has_error=false;
	for (unsigned int i=0; i< taken_maple.size(); ++i) {
	  if (!taken_maple[i]) has_error=true;
	}

	for (unsigned int i=0; i< taken_solver.size(); ++i) {
	  if (!taken_solver[i]) has_error=true;
	}

	if (has_error) {
	  std::cout << std::endl;
	  for (unsigned int i=0; i< taken_maple.size(); ++i) {
	    if (!taken_maple[i])
	      std::cerr << "NUMERIC ISSUE Missing " << roots_b[i] << std::endl;
	  }

	  for (unsigned int i=0; i< taken_solver.size(); ++i) {
	    if (!taken_solver[i]) std::cerr << "NUMERIC ISSUE Extra " << roots[i] << std::endl;
	  }
	}

	std::cout << std::endl;
      }
      ++reps;
    } while (timer.time() <3 && !verbose);

    timer.stop();

    if (verbose) {
      std::cout << "Elapsed time: " << timer.time()/reps << std::endl;

      std::cout << std::endl << std::endl;
    }
    else {
      num_roots += total_roots;
      total_time += timer.time();
    }

    //time_estimate(q, start, CE());
  }

  void clear_timings() {
    total_time=0;
    num_roots=0;
    ce_time=0;
    nce=0;
  }
  void write_timings() {
    if (!verbose) {
      std::cout.setf(std::ios::fixed);
      std::cout.precision(2);
      std::cout << 1000000.0*total_time/static_cast<double>(num_roots) ;
      if (ce_time != 0) {
	std::cout << "(" <<  1000000.0*ce_time/static_cast<double>(nce) << ")";
      }
      std::cout << " & ";
    }
  }

  void cleaned() { {
      // polynomials whose roots have been filtered based on derivitives
      Fn q= cf_(-16.614835192865492,19.540659228538015,-4);
      DV a; a.push_back(3.788873605);
      check_polynomial(q,a.begin(), a.end(), 1.096291202);
    }
    {
      Fn q= cf_(16.614835192865492,-19.540659228538015,4);
      DV a; a.push_back(1.09630); a.push_back(3.788873605);
      check_polynomial(q,a.begin(), a.end(), 1.096291202);
    }
    {
      Fn q= cf_(3.788873605, -3.788873605);
      DV a; a.push_back(1.0);
      check_polynomial(q,a.begin(), a.end(), 1);
    }
    {
      Fn q= cf_(-3.788873605, 3.788873605);
      DV a;
      check_polynomial(q,a.begin(), a.end(), 1);
    }
    {
      Fn q= cf_(-2295485086.0,2072822157.0, 116461914.2,  -116175036.5, -10063149.87,
		-196007.0344, 3460.886000, 136.9100396, 1.0);
      DV a;
      a.push_back(-79.34012316);
      a.push_back(-45.32403162);
      a.push_back(-32.43524488);
      a.push_back(-22.45230969);
      a.push_back(-5.239841933);
      a.push_back(1.123932833);
      a.push_back(3.435646759);
      a.push_back(43.32193209);
      check_polynomial(q, a.begin(), a.end(), Rt(-79.34012316));
      check_polynomial(q, a.begin()+2, a.end(), Rt(-45.32403162));
    }
  }

  void wilkinson() {
    clear_timings();

    unsigned int max_root = 20;
    CGAL_POLYNOMIAL_NS::Wilkinson_generator<K> gen(k_);
    Fn q=gen(max_root);
    DV a;
    for (unsigned int i = 1; i <= max_root; i++) {
      a.push_back(i);
    }
    check_polynomial(q,a.begin(), a.end());
    write_timings();
  }

  void mignotte() {
    clear_timings();

    unsigned int degree = 50;
    CGAL_POLYNOMIAL_NS::Mignotte_generator<K> gen(k_);
    Fn q=gen(degree);
    DV a;
    a.push_back(-1.0925395143411487);
    a.push_back(0.20000000000000001);
    a.push_back(0.20000000000000001);
    a.push_back(1.0756542734086822);
    check_polynomial(q,a.begin(), a.end());
    write_timings();
  }

  void square_free() {
    clear_timings();
    {
      Fn q= cf_(1,3,-4,-7);
      DV a; a.push_back(-0.8746784187); a.push_back(-0.2800178096);
      a.push_back(0.5832676569);
      check_polynomial(q,a.begin(), a.end());
    }
    {
      Fn q= cf_(2,0,1);
      DV a;
      check_polynomial(q,a.begin(), a.end());
    }
    /*
      {
      Fn q= Fn::make(0);
      DV a;
      check_square_free_polynomial(q,a.begin(), a.end(),print_multiplicity);
      }
      {
      Fn q= Fn::make(-3);
      DV a;
      check_square_free_polynomial(q,a.begin(), a.end(),print_multiplicity);
      }
    */
    {
      Fn q= cf_(2,-1);
      DV a;
      a.push_back(2);
      check_polynomial(q,a.begin(), a.end());
    }
    {
      Fn q= cf_(1,3,-4);
      DV a;a.push_back(-.25); a.push_back(1);
      check_polynomial(q,a.begin(), a.end());
    }
    {
      Fn q= cf_(2,5,-4,-7,8,-14);
      DV a;a.push_back(0.7382979184);
      check_polynomial(q,a.begin(), a.end());
    }
    {
      Fn q= cf_(2, 0, -1, 1, -3, -2);
      DV a;a.push_back(-1.833007477);
      a.push_back(-.8474637161);
      a.push_back(.7987710980);
      check_polynomial(q,a.begin(), a.end());
    }
    write_timings();
  }

  void small_intervals() {
    clear_timings();

    Fn q= cf_(-2295485086.0,2072822157.0, 116461914.2,  -116175036.5, -10063149.87,
	      -196007.0344, 3460.886000, 136.9100396, 1.0);
    DV a;
    a.push_back(-79.34012316);
    a.push_back(-45.32403162);
    a.push_back(-32.43524488);
    a.push_back(-22.45230969);
    a.push_back(-5.239841933);
    a.push_back(1.123932833);
    a.push_back(3.435646759);
    a.push_back(43.32193209);
    check_polynomial(q, a.begin()+7, a.end()-1, Rt(5.0), Rt(40.0));
    check_polynomial(q, a.begin()+5, a.end(), Rt(0.0), Rt(45.0));
    check_polynomial(q, a.begin()+6, a.end()-1, Rt(1.123932834), Rt(20.0));
    write_timings();
  }

  void non_simple() {
    clear_timings();
    {
      Fn q= cf_(-4365900, -5040, 1274073, 72979, -87647, -1045,
		2275, -79, -17, 1);
      // -9., -5., -3., -3., 2., 7., 7., 10., 11.
      /*DV a;a.push_back(RP(-9,1));
	a.push_back(RP(-5,1));
	a.push_back(RP(-3,2 ));
	a.push_back(RP(2,1));
	a.push_back(RP(7,2));
	a.push_back(RP(10,1));
	a.push_back(RP(11,1));*/
      DV a;a.push_back(-9);
      a.push_back(-5);
      a.push_back(-3);
      a.push_back(-3);
      a.push_back(2);
      a.push_back(7);
      a.push_back(7);
      a.push_back(10);
      a.push_back(11);
      check_polynomial(q,a.begin(), a.end());
      check_polynomial(q,a.begin()+2, a.end()-2, Rt(-4), Rt(9));
      check_polynomial(q,a.begin()+2, a.end()-4, Rt(-4), Rt(7));
      check_polynomial(q,a.begin()+2, a.end()-2, Rt(-3), Rt(9));
      check_polynomial(q,a.begin()+2, a.end()-4, Rt(-3), Rt(7));
    }
    write_timings();
  }

  void roots() {
    assert(std::numeric_limits<Rt>::has_infinity);
    Rt inf= std::numeric_limits<Rt>::infinity();
    Rt minf= -std::numeric_limits<Rt>::infinity();
    Rt big(100000);
    Rt small_rt(-1000000);
    Rt z(0);
    assert(big < inf);
    assert(big > minf);
    assert(small_rt < big);
    assert(z > small_rt);
    assert(-big < big);
    assert(small_rt < inf);
    assert(small_rt > minf);
    assert(z > minf);
    assert(z < inf);
  }

  void all() {
    square_free();
    wilkinson();
    mignotte();
    small_intervals();
    non_simple();
  }

  void exact() { {
      Fn q= cf_(NT("1225969745589853940699882447880155625/664613997892457936451903530140172288"),
		NT("-149192919905533219658325090294590900625/2658455991569831745807614120560689152"),
		NT("-176473560706138356181046066127379288725/2658455991569831745807614120560689152"),
		NT("-59323957331081368745121550919931250845/10633823966279326983230456482242756608"),
		NT("213976196148401335217623674372070736157/10633823966279326983230456482242756608"),
		NT("10107182444809889626671506586504631263/5316911983139663491615228241121378304"),
		NT("-12967429119484147645038507190192518603/10633823966279326983230456482242756608"),
		NT("-894074013945673684338239220768251181/10633823966279326983230456482242756608"),
		NT("963347446801307512325971073199501/5316911983139663491615228241121378304"));

      DV a;
      a.push_back(-14.46710600);
      a.push_back(-3.440766418);
      a.push_back(0.031679482287);
      a.push_back(2.562471103);
      a.push_back(3.490003117);
      a.push_back(478.0766178);
      check_polynomial(q, a.begin(), a.end());
      Root_container s=
	k_.root_container_object(q, Rt(-10000.0),
				 Rt(10000.0));

      std::vector<Rt> roots;
      std::copy(s.begin(), s.end(), std::back_inserter(roots));

      //typename Root_container::iterator it=s.begin();
      for (unsigned int offset=0; offset < roots.size()-1; ++offset) {
	check_polynomial(q, a.begin()+offset, a.end(), roots[offset]);
      }
    }

    {

      Fn q= cf_(NT("2102889623460648672193194772242161145741209860309/2854495385411919762116571938898990272765493248"),
		NT("-27785292706190294695336118341445556791078034998297/5708990770823839524233143877797980545530986496"),
		NT("26399617884277473686036041581122386808181680932651/2854495385411919762116571938898990272765493248"),
		NT("-65515073263163417483134441553521565/5070602400912917605986812821504"),
		NT("6369956942652855102749740125692887/316912650057057350374175801344"),
		NT("-19558908001688696225/2251799813685248"),
		NT("5694238778856570877/562949953421312"),
		NT("-6468/1"));

      DV a;
      a.push_back(0.2307182653);
      a.push_back(0.4663233325);
      a.push_back(1.602781661);
      check_polynomial(q, a.begin(), a.end(), 0);
      Root_container s=
	k_.root_container_object(q, -Rt(10000.0),
				 Rt(100000.0));

      std::vector<Rt> roots;
      std::copy(s.begin(), s.end(), std::back_inserter(roots));

      //typename Root_container::iterator it=s.begin();
      for (unsigned int offset=0; offset < roots.size()-1; ++offset) {
	check_polynomial(q, a.begin()+offset, a.end(), roots[offset]);
      }

    }
  }

};

//typedef CGAL::Polynomial::Polynomial_kernel<Polynomial_gmpq> Kernel_gmpq;
//typedef CGAL::Polynomial::Polynomial_kernel<Polynomial_double> Kernel_double;
/*typedef CGAL::Polynomial::Filtered_kernel<Polynomial_gmpq,
  Polynomial_double,
  CGAL::Polynomial::Polynomial<CGAL::Interval_nt_advanced>,
  Polynomial_gmpq> Filtered_kernel_gmpq;*/
//typedef CGAL::Polynomial::Polynomial_kernel<Polynomial_bigint> Kernel_core;

//typedef Kernel_gmpq Exact_kernel;
//typedef Filtered_kernel_gmpq Filtered_kernel;
//typedef Kernel_double Numeric_kernel;
#endif
