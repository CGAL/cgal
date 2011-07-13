/*
	Multi-precision real number class. C++ wrapper fo MPFR library.
	Project homepage: http://www.holoborodko.com/pavel/
	Contact e-mail:   pavel@holoborodko.com

	Copyright (c) 2008-2010 Pavel Holoborodko

	This library is free software; you can redistribute it and/or
	modify it under the terms of the GNU Lesser General Public
	License as published by the Free Software Foundation; either
	version 2.1 of the License, or (at your option) any later version.

	This library is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
	Lesser General Public License for more details.

	You should have received a copy of the GNU Lesser General Public
	License along with this library; if not, write to the Free Software
	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

	Contributors:
	Brian Gladman, Helmut Jarausch, Fokko Beekhof, Ulrich Mutze, 
	Heinz van Saanen, Pere Constans, Dmitriy Gubanov
*/

#include <cstring>
#include "mpreal.h"

using std::ws;
using std::cerr;
using std::endl;
using std::string;
using std::ostream;
using std::istream;

namespace mpfr{

mp_rnd_t   mpreal::default_rnd  = mpfr_get_default_rounding_mode();	
mp_prec_t  mpreal::default_prec = mpfr_get_default_prec();	
int		   mpreal::default_base = 10;
int        mpreal::double_bits = -1;

// Default constructor: creates mp number and initializes it to 0.
mpreal::mpreal() 
{ 
	mpfr_init2(mp,default_prec); 
	mpfr_set_ui(mp,0,default_rnd);
}

mpreal::mpreal(const mpreal& u) 
{
	mpfr_init2(mp,mpfr_get_prec(u.mp));
	mpfr_set(mp,u.mp,default_rnd);
}

mpreal::mpreal(const mpfr_t u)
{
	mpfr_init2(mp,mpfr_get_prec(u));
	mpfr_set(mp,u,default_rnd);
}

mpreal::mpreal(const mpf_t u)
{
	mpfr_init2(mp,mpf_get_prec(u));
	mpfr_set_f(mp,u,default_rnd);
}

mpreal::mpreal(const mpz_t u, mp_prec_t prec, mp_rnd_t mode)
{
	mpfr_init2(mp,prec);
	mpfr_set_z(mp,u,mode);
}

mpreal::mpreal(const mpq_t u, mp_prec_t prec, mp_rnd_t mode)
{
	mpfr_init2(mp,prec);
	mpfr_set_q(mp,u,mode);
}

mpreal::mpreal(const double u, mp_prec_t prec, mp_rnd_t mode)
{
    if(double_bits == -1 || fits_in_bits(u, double_bits))
    {
    	mpfr_init2(mp,prec);
	    mpfr_set_d(mp,u,mode);
    }
    else
        throw conversion_overflow();
}

mpreal::mpreal(const long double u, mp_prec_t prec, mp_rnd_t mode)
{ 
    mpfr_init2(mp,prec);
	mpfr_set_ld(mp,u,mode);
}

mpreal::mpreal(const unsigned long int u, mp_prec_t prec, mp_rnd_t mode)
{ 
	mpfr_init2(mp,prec);
	mpfr_set_ui(mp,u,mode);
}

mpreal::mpreal(const unsigned int u, mp_prec_t prec, mp_rnd_t mode)
{ 
	mpfr_init2(mp,prec);
	mpfr_set_ui(mp,u,mode);
}

mpreal::mpreal(const long int u, mp_prec_t prec, mp_rnd_t mode)
{ 
	mpfr_init2(mp,prec);
	mpfr_set_si(mp,u,mode);
}

mpreal::mpreal(const int u, mp_prec_t prec, mp_rnd_t mode)
{ 
	mpfr_init2(mp,prec);
	mpfr_set_si(mp,u,mode);
}

mpreal::mpreal(const char* s, mp_prec_t prec, int base, mp_rnd_t mode)
{
	mpfr_init2(mp,prec);
	mpfr_set_str(mp, s, base, mode); 
}

mpreal::~mpreal() 
{ 
	mpfr_clear(mp); 
}                           

// Operators - Assignment
mpreal& mpreal::operator=(const char* s)
{
	mpfr_t t;
	if(0==mpfr_init_set_str(t,s,default_base,default_rnd))
	{
		mpfr_set(mp,t,mpreal::default_rnd);				
		mpfr_clear(t);
	}else{
		mpfr_clear(t);
		// cerr<<"fail to convert string"<<endl;
	}

	return *this;
}

const mpreal fma (const mpreal& v1, const mpreal& v2, const mpreal& v3, mp_rnd_t rnd_mode)
{
	mpreal a;
	mp_prec_t p1, p2, p3;

	p1 = v1.get_prec(); 
	p2 = v2.get_prec(); 
	p3 = v3.get_prec(); 

	a.set_prec(p3>p2?(p3>p1?p3:p1):(p2>p1?p2:p1));

	mpfr_fma(a.mp,v1.mp,v2.mp,v3.mp,rnd_mode);
	return a;
}

const mpreal fms (const mpreal& v1, const mpreal& v2, const mpreal& v3, mp_rnd_t rnd_mode)
{
	mpreal a;
	mp_prec_t p1, p2, p3;

	p1 = v1.get_prec(); 
	p2 = v2.get_prec(); 
	p3 = v3.get_prec(); 

	a.set_prec(p3>p2?(p3>p1?p3:p1):(p2>p1?p2:p1));

	mpfr_fms(a.mp,v1.mp,v2.mp,v3.mp,rnd_mode);
	return a;
}

const mpreal agm (const mpreal& v1, const mpreal& v2, mp_rnd_t rnd_mode)
{
	mpreal a;
	mp_prec_t p1, p2;

	p1 = v1.get_prec(); 
	p2 = v2.get_prec(); 

	a.set_prec(p1>p2?p1:p2);

	mpfr_agm(a.mp, v1.mp, v2.mp, rnd_mode);

	return a;
}

const mpreal hypot (const mpreal& x, const mpreal& y, mp_rnd_t rnd_mode)
{
	mpreal a;
	mp_prec_t yp, xp;

	yp = y.get_prec(); 
	xp = x.get_prec(); 

	a.set_prec(yp>xp?yp:xp);

	mpfr_hypot(a.mp, x.mp, y.mp, rnd_mode);

	return a;
}

const mpreal sum (const mpreal tab[], unsigned long int n, mp_rnd_t rnd_mode)
{
	mpreal x;
	mpfr_ptr* t;
	unsigned long int i;

	t = new mpfr_ptr[n];
	for (i=0;i<n;i++) t[i] = (mpfr_ptr)tab[i].mp;
	mpfr_sum(x.mp,t,n,rnd_mode);
	delete[] t;
	return x;
}

const mpreal remainder (const mpreal& x, const mpreal& y, mp_rnd_t rnd_mode)
{	
	mpreal a;
	mp_prec_t yp, xp;

	yp = y.get_prec(); 
	xp = x.get_prec(); 

	a.set_prec(yp>xp?yp:xp);

	mpfr_remainder(a.mp, x.mp, y.mp, rnd_mode);

	return a;
}

const mpreal remquo (long* q, const mpreal& x, const mpreal& y, mp_rnd_t rnd_mode)
{
	mpreal a;
	mp_prec_t yp, xp;

	yp = y.get_prec(); 
	xp = x.get_prec(); 

	a.set_prec(yp>xp?yp:xp);

	mpfr_remquo(a.mp,q, x.mp, y.mp, rnd_mode);

	return a;
}

template <class T>
std::string to_string(T t, std::ios_base & (*f)(std::ios_base&))
{
	std::ostringstream oss;
	oss << f << t;
	return oss.str();
}

mpreal::operator std::string() const
{
	return to_string();
}

string mpreal::to_string(size_t n, int b, mp_rnd_t mode) const
{
	char *s, *ns = NULL;	
	size_t slen, nslen;
	mp_exp_t exp;
	string out;

	if(mpfr_inf_p(mp))
	{ 
		if(mpfr_sgn(mp)>0) return "+@Inf@";
		else			   return "-@Inf@";
	}

	if(mpfr_zero_p(mp)) return "0";
	if(mpfr_nan_p(mp))  return "@NaN@";

		
	s  = mpfr_get_str(NULL,&exp,b,0,mp,mode);
	ns = mpfr_get_str(NULL,&exp,b,n,mp,mode);

	if(s!=NULL && ns!=NULL)
	{
		slen  = strlen(s);
		nslen = strlen(ns);
		if(nslen<=slen) 
		{
			mpfr_free_str(s);
			s = ns;
			slen = nslen;
		}
		else {
			mpfr_free_str(ns);
		}

		// Make human eye-friendly formatting if possible
		if (exp>0 && static_cast<size_t>(exp)<slen)
		{
			if(s[0]=='-')
			{
				// Remove zeros starting from right end
				char* ptr = s+slen-1;
				while (*ptr=='0' && ptr>s+exp) ptr--; 

				if(ptr==s+exp) out = string(s,exp+1);
				else		   out = string(s,exp+1)+'.'+string(s+exp+1,ptr-(s+exp+1)+1);

				//out = string(s,exp+1)+'.'+string(s+exp+1);
			}
			else
			{
				// Remove zeros starting from right end
				char* ptr = s+slen-1;
				while (*ptr=='0' && ptr>s+exp-1) ptr--; 

				if(ptr==s+exp-1) out = string(s,exp);
				else		     out = string(s,exp)+'.'+string(s+exp,ptr-(s+exp)+1);

				//out = string(s,exp)+'.'+string(s+exp);
			}

		}else{ // exp<0 || exp>slen
			if(s[0]=='-')
			{
				// Remove zeros starting from right end
				char* ptr = s+slen-1;
				while (*ptr=='0' && ptr>s+1) ptr--; 

				if(ptr==s+1) out = string(s,2);
				else		 out = string(s,2)+'.'+string(s+2,ptr-(s+2)+1);

				//out = string(s,2)+'.'+string(s+2);
			}
			else
			{
				// Remove zeros starting from right end
				char* ptr = s+slen-1;
				while (*ptr=='0' && ptr>s) ptr--; 

				if(ptr==s) out = string(s,1);
				else	   out = string(s,1)+'.'+string(s+1,ptr-(s+1)+1);

				//out = string(s,1)+'.'+string(s+1);
			}

			// Make final string
			if(--exp)
			{
				if(exp>0) out += "e+"+mpfr::to_string<mp_exp_t>(exp,std::dec);
				else 	  out += "e"+mpfr::to_string<mp_exp_t>(exp,std::dec);
			}
		}

		mpfr_free_str(s);
		return out;
	}else{
		return "conversion error!";
	}
}

//////////////////////////////////////////////////////////////////////////
// I/O
ostream& operator<<(ostream& os, const mpreal& v)
{
	return os<<v.to_string(os.precision());
}

istream& operator>>(istream &is, mpreal& v)
{
	char c;	
	string s = "";
	mpfr_t t;

	if(is.good())
	{
		is>>ws;
		while ((c = is.get())!=EOF)
		{
			if(c ==' ' || c == '\t' || c == '\n' || c == '\r')
			{
				is.putback(c);
				break;
			}
			s += c;
		}

		if(s.size() != 0)
		{
			// Protect current value from alternation in case of input error
			// so some error handling(roll back) procedure can be used 
			if(0==mpfr_init_set_str(t,s.c_str(),mpreal::default_base,mpreal::default_rnd))
			{
				mpfr_set(v.mp,t,mpreal::default_rnd);				
				mpfr_clear(t);

			}else{
				mpfr_clear(t);
				cerr<<"error reading from istream"<<endl;
				// throw an exception
			}
		}
	}
	return is;
}
}
