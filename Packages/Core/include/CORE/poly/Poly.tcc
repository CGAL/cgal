/******************************************************************
 * Core Library Version 1.6, June 2003
 * Copyright (c) 1995-2003 Exact Computation Project

 File: Poly.tcc
 Purpose: 
	Template implementations of the functions
	of the Polynomial<NT> class (found in Poly.h)

 OVERVIEW:  
	Each polynomial has a nominal "degree" (this
		is an upper bound on the true degree, which
		is determined by the first non-zero coefficient).
	Coefficients are parametrized by some number type "NT".
	Coefficients are stored in the "coeff" array of
		length "degree + 1".  
	IMPORTANT CONVENTION: the zero polynomial has degree -1
		while nonzero constant polynomials have degree 0.
	
 Author: Chee Yap, Sylvain Pion and Vikram Sharma
 Date:   May 28, 2002
 Core Library
 $Id$
 ************************************** */ 

template <class NT>
  int Polynomial<NT>::COEFF_PER_LINE  = 3;           // pretty print parameters
template <class NT>
  const char* Polynomial<NT>::INDENT_SPACE ="   ";  // pretty print parameters

// ==================================================
// Polynomial Constructors
// ==================================================

template <class NT>
  Polynomial<NT>::Polynomial(void) {
    degree = -1; 	// this is the zero polynomial!
    coeff = NULL;
  }

template <class NT>
  Polynomial<NT>::Polynomial(int n) {
    assert(n>= -1);
    degree = n;
    if (n == -1) return;	// return the zero polynomial!
    if (n>=0) coeff = new NT[n+1];
    coeff[0]=1;			// otherwise, return the unity polynomial
    for (int i=1; i<=n; i++)
      coeff[i]=0;
  }

template <class NT>
  Polynomial<NT>::Polynomial(int n, NT * c) {
    assert("array c has n+1 elements");
    degree = n; 
    if (n >= 0) {
	coeff = new NT[n+1];
	for (int i = 0; i <= n; i++) 
	  coeff[i] = c[i];
    }
  }

//
// Constructor from a vector of NT's
///////////////////////////////////////
template <class NT>
  Polynomial<NT>::Polynomial(const VecNT & vN) {
    degree = vN.size()-1; 
    if (degree >= 0) {
	coeff = new NT[degree+1];
	for (int i = 0; i <= degree; i++) 
	  coeff[i] = vN[i];
    }
  }

template <class NT>
  Polynomial<NT>::Polynomial(const Polynomial<NT> & p) {
    coeff = NULL;
    *this = p;	// reduce to assignment operator=
  }

// Constructor from a Character string of coefficients
///////////////////////////////////////
template <class NT>
  Polynomial<NT>::Polynomial(int n, const char * s[]) {
    assert("array s has n+1 elements");
    degree = n;
    if (n >= 0) {
        coeff = new NT[n+1];
        for (int i = 0; i <= n; i++)
          coeff[i] = s[i];
    }
  }

template <class NT>
  Polynomial<NT>::~Polynomial() {
     if (degree >= 0)
	delete[] coeff;
  }

  // assignment:
template <class NT>
  Polynomial<NT> & Polynomial<NT>::operator=(const Polynomial<NT>& p) {
    if (this == &p) return *this;	// self-assignment
    degree = p.getDegree();
    delete[] coeff;
    coeff = new NT[degree +1];
    for (int i = 0; i <= degree; i++)
	coeff[i] = p.coeff[i];
    return *this;
  }

  // getTrueDegree
template <class NT>
  int Polynomial<NT>::getTrueDegree() const {
    for (int i=degree; i>=0; i--)
	if (coeff[i] != 0) return i;
    return -1;	// Zero polynomial
  }

// ==================================================
// polynomial arithmetic
// ==================================================

  // Expands the nominal degree to n;  
  //	Returns n if nominal degree is changed to n
  //	Else returns -2
template <class NT>
  int Polynomial<NT>::expand(int n) { 
    if ((n <= degree)||(n < 0)) return -2;
    int i;
    NT * c = coeff;
    coeff = new NT[n+1];
    for (i = 0; i<= degree; i++) 
      coeff[i] = c[i];
    for (i = degree+1; i<=n; i++)
      coeff[i] = 0;
    delete[] c;
    degree = n;
    return n;
  }

  // contract() gets rid of leading zero coefficients
  //	and returns the new (true) degree;
  //	It returns -2 if this is a no-op
template <class NT>
  int Polynomial<NT>::contract() { 
    int d = getTrueDegree();
    if (d == degree) return (-2);  // nothing to do
    else degree = d;
    NT * c = coeff;
    coeff = new NT[d+1];
    for (int i = 0; i<= d; i++)
	coeff[i] = c[i];
    delete[] c;
    return d;
  }

template <class NT>
  Polynomial<NT> & Polynomial<NT>::operator+=(const Polynomial<NT>& p) { // +=
    int d = p.getDegree();
    if (d > degree) expand(d);
    for (int i = 0; i<=d; i++)
	  coeff[i] += p.coeff[i];
    return *this;
  }
template <class NT>
  Polynomial<NT> & Polynomial<NT>::operator-=(const Polynomial<NT>& p) { // -=
    int d = p.getDegree();
    if (d > degree) expand(d);
    for (int i = 0; i<=d; i++)
	  coeff[i] -= p.coeff[i];
    return *this;
  }

  // SELF-MULTIPLICATION
  // This is quadratic time multiplication!
template <class NT>
  Polynomial<NT> & Polynomial<NT>::operator*=(const Polynomial<NT>& p) { // *=
    int d = degree + p.getDegree();
    NT * c = new NT[d+1];
    for (int i = 0; i<=d; i++)
	c[i] = 0;

    for (int i = 0; i<=p.getDegree(); i++)
      for (int j = 0; j<=degree; j++) {
	  c[i+j] += p.coeff[i] * coeff[j];
      }
    degree = d;
    delete[] coeff;
    coeff = c;
    return *this;
  }

  // Multiply by a scalar
template <class NT>
  Polynomial<NT> & Polynomial<NT>::mulScalar( const NT & c) {
    for (int i = 0; i<=degree ; i++)
	coeff[i] *= c;
    return *this;
  }

  // mulXpower: Multiply by X^i (COULD be a divide if i<0)
template <class NT>
  Polynomial<NT> & Polynomial<NT>::mulXpower(int s) {
		// if s >= 0, then this is equivalent to 
		// multiplying by X^s;  if s < 0, to dividing by X^s
    if (s==0) return *this;
    int d = s+getTrueDegree();
    if (d < 0) {
	degree = -1; delete[] coeff;
        coeff = NULL;
	return *this;
    }
    NT * c = new NT[d+1];
    if (s>0) 
	for (int j=0;  j <= d; j++) {
	  if (j <= degree)
		c[d-j] = coeff[d-s-j];
	  else
		c[d-j] = 0;
	}
    if (s<0) {
	for (int j=0; j <= d; j++)
	  c[d-j] = coeff[d-s-j];  // since s<0, (d-s-j) > (d-j)
    }
    delete[] coeff;
    coeff = c; degree = d;
    return *this;
  }//mulXpower

  // REDUCE STEP (helper for PSEUDO-REMAINDER function)
  // Let THIS=(*this) be the current polynomial, and P be the input
  //	argument for reduceStep.  Let R be returned polynomial.
  //	R has the special form as a binomial,
  //		R = C + M
  //	where C is a constant and M a monomial (= coeff * power of X).
  //	Moreover, THIS is transsformed to a new polynomial THAT which
  //	is given by
  // 		THAT =  (C * THIS) - M * P
  //	MOREOVER: deg(THAT) < deg(THIS) unless deg(P)>deg(Q).  
  //	Basically, C is a power of the leading coefficient of P.
template <class NT>
  Polynomial<NT> Polynomial<NT>::reduceStep (
		Polynomial<NT>& p) {
	p.contract(); contract();	// first contract both polynomials
	Polynomial<NT> q(p);		// q is initially a copy of p
					//	but is eventually M*P
	int pDeg  = q.degree; int myDeg = degree;
	if (pDeg == -1)  return *(new Polynomial());  // Zero Polynomial
			// NOTE: pDeg=0 is really an error condition!
	if (myDeg < pDeg) return *(new Polynomial(0));  // Unity Polynomial
			// i.e., C=1, M=0.
	// Now (myDeg >= pDeg).  Start to form the Return Polynomial R=C+M
	Polynomial<NT> R(myDeg - pDeg + 1);  // deg(M)= myDeg - pDeg
	q.mulXpower(myDeg - pDeg);  	 // q is going to become M*P

	NT myLC = coeff[myDeg];	  // LC means "leading coefficient"
	NT qLC = q.coeff[myDeg];  // p also has degree "myDeg" (qLC non-zero)
	NT LC;
	if ((myLC % qLC) == 0) { // myLC is divisible by qLC 
	    LC = myLC / qLC;	 // exact division, LC not zero
	    R.setCoeff(0, 1);  		  //  C = 1,
	    if (LC != 1) 		  // Just a little optimization
	        R.setCoeff(R.degree, LC); //  M = LC * X^(myDeg-pDeg)
		q.mulScalar(LC); 	  //  q = M*P.
	}
	if ((qLC % myLC) == 0) {  // qLC is divisible by myLC 
	    LC = qLC / myLC;	  // exact division, LC not zero
	    if ((LC != 1) && (LC != -1)) { // IMPORTANT: unlike the previous
				  // case, we need an additional condition
				  // that LC != -1.  THis is because 
				  // if (LC = -1), then we have qLC and
				  // myLC are mutually divisible, and
				  // we would be updating R twice! 
		R.setCoeff(0, LC); 	   // C = LC, M = X^(myDeg-pDeg)
		mulScalar(LC); 	   	   // And q = M*P.
	    }
	}
	if (((qLC % myLC) != 0) && ((myLC % qLC)!= 0)) {  //
	  // NT g = gcd(qLC, myLC); 	// This ASSUMES gcd is defined in NT !!
	  NT g = 1;  			// In case no gcd is available
	  if (g == 1) {
	  	R.setCoeff(0, qLC);	  	// C = qLC
	  	R.setCoeff(R.degree, myLC);	 // M = (myLC) * X^{myDeg-pDeg}
	  	mulScalar(qLC);	 		// forming  C * THIS 
	  	q.mulScalar(myLC);		// forming  M * P
	  } else {
	  	R.setCoeff(0, qLC/g);	  	// C = qLC/g
	  	R.setCoeff(R.degree, myLC/g);	// M = (myLC/g) * X^{myDeg-pDeg}
	  	mulScalar(qLC/g);	 	// forming  C * THIS 
	  	q.mulScalar(myLC/g);		// forming  M * P
	  }
	} 
	(*this) -= q;		// THAT is finally C*THIS - M*P
	contract();		// Leading term is eliminated!
	return R;		// Returns R = C + M
  }// reduceStep

  // PSEUDO-REMAINDER (and PSEUDO-QUOTIENT)
  // Let THIS = (*this) be input polynomial and it is eventually
  // transformed to THAT.  Let R be the returned polynomial,
  // and P be the input polynomial.  Then we have
  // 		THAT =  (C * THIS) - R * P
  // where deg(THAT) < deg(P).
  // Moreover, C is uniquely determined (we won't spell it out)
  // except to note that
  //	C divides D = (LC)^{deg(THIS)-deg(P)+1} 
  //	where LC is the leading coefficient of P.
  // In this program, C is stored as the static variable ccc_ (HACK)
  // NOTE: Normally, Pseudo-Remainder is defined so that C = D

template <class NT>
  Polynomial<NT> Polynomial<NT>::pseudoRemainder (
		Polynomial<NT>& p) {
    contract(); p.contract();
    Polynomial<NT> q(p);
    if (q.degree == -1)  {
	std::cout << "ERROR in Polynomial<NT>::pseudoRemainder :\n" <<
		"    -- divide by zero polynomial" << std::endl;
	return Polynomial(0);  // Unit Polynomial (arbitrary!)
    }
    if (q.degree > degree) {
	return Polynomial(); // Zero Polynomial
		// CHECK: THAT = 1 * THIS - 0 * P,  deg(THAT) < deg(P)
    }
    Polynomial<NT> R;  	// accumulating the return polynomial R
    Polynomial<NT> tmpR;
    int sgn = -1;
//    ccc_ = *(new NT(1));  // ccc_ is a static NT variable (hack)
    while (degree >= p.degree) {
        q = p;
	tmpR = reduceStep(q);		// tmpR = C + M
        if (tmpR.coeff[0] < 0)
          sgn = -sgn;
	R.mulScalar(tmpR.coeff[0]);	// R *= C
//	ccc_ *= tmpR.coeff[0];		// hack!!
	tmpR.setCoeff(0,0);		// remove the zeroth coeff
	R += tmpR.mulXpower(-1);	// R += M   (so, R -> C*R + M)
    }

    if (sgn<0) *this = -*this;
    return R;	// R is the pseudo-quotient
  }//pseudoRemainder

template <class NT>
  Polynomial<NT> & Polynomial<NT>::operator-() {	// unary minus
    for (int i=0; i<=degree; i++)
	coeff[i] *= -1;
    return *this;
  }
template <class NT>
  Polynomial<NT> & Polynomial<NT>::power(unsigned int n) {	// self-power
    if (n == 0) {
	degree = 0;
	delete [] coeff;
	coeff = new NT[1]; coeff[0] = 1;
    } else {
      Polynomial<NT> p = *this;
      for (unsigned int i=1; i<n; i++)
	*this *= p;
    }
    return *this;
  }

  // evaluation of Expr value
template <class NT>
  Expr Polynomial<NT>::eval(const Expr& e) const {		// evaluation
    if (degree == -1) return Expr(0);
    if (degree == 0) return Expr(coeff[0]);
    Expr val(0);
    for (int i=degree; i>=0; i--) {
	val *= e;
	val += coeff[i];
    }
    return val;
  }//eval

  // evaluation of BigFloat value
template <class NT>
  BigFloat Polynomial<NT>::eval(const BigFloat& f) const {		// evaluation
    if (degree == -1) return BigFloat(0);
    if (degree == 0) return BigFloat(coeff[0]);
    BigFloat val(0);
    for (int i=degree; i>=0; i--) {
	val *= f;
	val.makeExact();	// in future, this should NOT be necessary
	val += coeff[i];
	val.makeExact();	// in future, this should NOT be necessary
    }
    return val;
  }//eval

  //============================================================
  // Bounds
  //============================================================

  // Cauchy Upper Bound on Roots
template < class NT >
  BigFloat Polynomial<NT>::CauchyUpperBound() const {
    if (zeroP(*this))
      return BigFloat(0);
    Expr mx = 0;
    int deg = getTrueDegree();
    for (int i = 0; i < deg; ++i) {
      mx = core_max(mx, Expr(core_abs(coeff[i])));
    }
    mx /= core_abs(coeff[deg]);
    mx.approx(CORE_INFTY, 2); 
    	// get an absolute approximate value with error < 1/4
    return (mx.BigFloatValue().makeExact() + 2);
}
  // Cauchy Lower Bound on Roots
template < class NT >
  BigFloat Polynomial<NT>::CauchyLowerBound() const {
    if (zeroP(*this))
      return BigFloat(0);
    Expr mx = 0;
    int deg = getTrueDegree();
    for (int i = 1; i <= deg; ++i) {
      mx = core_max(mx, Expr(core_abs(coeff[i])));
    }
    mx = core_abs(coeff[0])/(core_abs(coeff[0]) + mx) ;
    mx.approx(2, CORE_INFTY); 
    	// get an relative approximate value with error < 1/4
    return (mx.BigFloatValue().makeExact().div2());
}

  // Separation bound for polynomials that may have multiple roots
  // We use the Rump-Schwartz bound
template < class NT >
  BigFloat Polynomial<NT>::sepBound() const {
    BigInt d;
    BigFloat e;
    int deg = getTrueDegree();
    
    CORE::power(d, BigInt(deg), ((deg)+5)/2);
    e = CORE::power(height()+1, deg);
    return 1/e*2*d;
  }

  // height function
template < class NT >
  BigFloat Polynomial<NT>::height() const {
    if (zeroP(*this))
	return NT(0);
    int deg = getTrueDegree();
    NT ht = 0;
    for (int i = 0; i< deg; i++) 
	ht = std::max(ht, core_abs(coeff[i]) );
    return BigFloat(ht); 
  }

  // length function
template < class NT >
  BigFloat Polynomial<NT>::length() const {
   if (zeroP(*this))
	return NT(0);
    int deg = getTrueDegree();
    NT length = 0;
    for (int i = 0; i< deg; i++) 
	length += core_abs(coeff[i]*coeff[i]);
    return sqrt(BigFloat(length));
  }

  //============================================================
  // differentiation
  //============================================================

template <class NT>
  Polynomial<NT> & Polynomial<NT>::differentiate() {	// self-differentiation
    if (degree >= 0) {
      NT * c = new NT[degree];
      for (int i=1; i<=degree; i++) 
	c[i-1] = coeff[i] * i;
      degree--;
      delete[] coeff; 
      coeff = c;
    }
    return *this;
  }// differentiation

  // multi-differentiate
template <class NT>
  Polynomial<NT> & Polynomial<NT>::differentiate(int n) {
    assert(n >= 0);
    for (int i=1; i<=n; i++)
	this->differentiate();
    return *this;
  } // multi-differentiate


// ==================================================
// Square-free and Primitive Parts
// ==================================================

// squareFree not yet implemented


// Primitive Part:  *this is changed
//
template <class NT>
  Polynomial<NT> & Polynomial<NT>::primPart(){
  	// NOTE: GCD must be provided by NT
    int d = getTrueDegree();
    assert (d >= 0);
    if (d == 0) {
	if (coeff[0] > 0)
	   coeff[0] = 1;
	else
	   coeff[0] = -1;
	return *this;
    }

    NT g = coeff[d];
    NT c;
    for (int i=d-1; i>=0; i--) {
    	c = coeff[i];
        g = gcd(g, c);
    }
    if (g==1) return *this;

    for (int i=0; i<=d; i++){
    	coeff[i] =  coeff[i]/g;
    }
    return *this;
  }// primPart


// ==================================================
// Useful member functions
// ==================================================

  // reverse:
  // 	reverses the list of coefficients
template <class NT>
  void Polynomial<NT>::reverse() { 
    NT tmp;
    for (int i=0; i<= degree/2; i++) {
  	tmp = coeff[i];
  	coeff[i] =   coeff[degree-i];
  	coeff[degree-i] = tmp;
    }
  }//reverse

  // dump: prints out polynomial
template <class NT>
  void Polynomial<NT>::dump() const {
    std::cerr << "Polynomial: deg = " << degree << std::endl;
    int termsInLine = 1;
    int i=0;
    for (; i<= getTrueDegree();  ++i) 
	if (coeff[i] != 0) break;
    
    std::cerr << Polynomial<NT>::INDENT_SPACE << "(" << coeff[i] << ")";
    if (i>0)
      std::cerr << " * x ^ " << i;


    for (i++ ; i<= getTrueDegree(); ++i) {
      if (coeff[i] == 0) continue;
      termsInLine++;
      if (termsInLine % Polynomial<NT>::COEFF_PER_LINE == 0)
	  std::cerr << std::endl << Polynomial<NT>::INDENT_SPACE ;
      std::cerr << " + (" << coeff[i] << ")";
      std::cerr << " * x ^ " << i;
    }
    std::cerr << std::endl;
  }//dump

  // dump: with optional message string
template <class NT>
  void Polynomial<NT>::dump(std::string m1) const{
    std::cerr << m1;
    dump();
  }

  // Dump of Maple Code for Polynomial
template <class NT>
  void Polynomial<NT>::mapleDump() const {
    if (zeroP(*this)) {
      std::cerr << 0 << std::endl;
      return;
    }
    std::cerr << coeff[0];
    for (int i = 1; i<= getTrueDegree(); ++i) {
      std::cerr << " + (" << coeff[i] << ")";
      std::cerr << "*x^" << i;
    }
    std::cerr << std::endl;
  }//mapleDump

// ==================================================
// Useful friend functions for Polynomial<NT> class
// ==================================================

  // friend differentiation
template <class NT>
  Polynomial<NT> differentiate(const Polynomial<NT> & p) {	  // differentiate
    Polynomial<NT> q(p);
    return q.differentiate();
  }

  // friend multi-differentiation
template <class NT>
  Polynomial<NT> differentiate(const Polynomial<NT> & p, int n){//multi-differentiate
    Polynomial<NT> q(p);
    assert(n >= 0);
    for (int i=1; i<=n; i++)
	q.differentiate();
    return q;
  }

  // friend equality comparison
template <class NT>
  bool operator==(const Polynomial<NT>& p, const Polynomial<NT>& q) {	// ==
    int d, D;
    Polynomial<NT> P(p); P.contract();
    Polynomial<NT> Q(q); Q.contract();
    if (P.degree < Q.degree) {
	d = P.degree; D = Q.degree;
	for (int i = d+1; i<=D; i++)
	  if (Q.coeff[i] != 0) return false;	// return false
    } else {
	D = P.degree; d = Q.degree; 
	for (int i = d+1; i<=D; i++)
	  if (P.coeff[i] != 0) return false;	// return false
    }
    for (int i = 0; i <= d; i++)
	if (P.coeff[i] != Q.coeff[i]) return false;	// return false
    return true; 	// return true
  }

  // friend non-equality comparison
template <class NT>
  bool operator!=(const Polynomial<NT>& p, const Polynomial<NT>& q) {	// !=
    return (!(p == q));
  }

  // stream i/o
template <class NT>
  std::ostream& operator<<(std::ostream& o, const Polynomial<NT>& p) {
    o <<   "Polynomial<NT> ( deg = " << p.degree ;
    if (p.degree >= 0) {
      o << "," << std::endl;
      o << ">  coeff c0,c1,... = " << p.coeff[0];
      for (int i=1; i<= p.degree; i++)
	o << ", " <<  p.coeff[i] ;
    }
    o << ")" << std::endl;
    return o;
  }

  // fragile version...
template <class NT>
  std::istream& operator>>(std::istream& is, Polynomial<NT>& p) {
    is >> p.degree; 
    p.coeff = new NT[p.degree+1];
    for (int i=0; i<= p.degree; i++)
      is >> p.coeff[i];
    return is;
  }

// ==================================================
// Simple test of poly
// ==================================================

template <class NT>
bool testPoly() {
  int c[] = {1, 2, 3};
  Polynomial<NT> p(2, c);
  std::cout << p;

  Polynomial<NT> zeroP;
  std::cout << "zeroP  : " << zeroP << std::endl;

  Polynomial<NT> P5(5);
  std::cout << "Poly 5 : " << P5 << std::endl;

  return 0;
}

// ==================================================
// End of Polynomial<NT>
// ==================================================
