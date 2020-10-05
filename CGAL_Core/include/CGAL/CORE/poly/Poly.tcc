/*
 * Core Library Version 1.7, August 2004
 * Copyright (c) 1995-2004 Exact Computation Project
 * All rights reserved.
 *
 * This file is part of CGAL (www.cgal.org).
 *
 * File: Poly.tcc
 * Purpose:
 *        Template implementations of the functions
 *        of the Polynomial<NT> class (found in Poly.h)
 *
 * OVERVIEW:
 *        Each polynomial has a nominal "degree" (this
 *                is an upper bound on the true degree, which
 *                is determined by the first non-zero coefficient).
 *        Coefficients are parametrized by some number type "NT".
 *        Coefficients are stored in the "coeff" array of
 *                length "degree + 1".
 *        IMPORTANT CONVENTION: the zero polynomial has degree -1
 *                while nonzero constant polynomials have degree 0.
 *
 * Bugs:
 *        Currently, coefficient number type NT only accept
 *                        NT=BigInt and NT=int
 *
 *        To see where NT=Expr will give trouble,
 *                        look for NOTE_EXPR below.
 *
 * Author: Chee Yap, Sylvain Pion and Vikram Sharma
 * Date:   May 28, 2002
 *
 * WWW URL: http://cs.nyu.edu/exact/
 * Email: exact@cs.nyu.edu
 *
 * $URL$
 * $Id$
 * SPDX-License-Identifier: LGPL-3.0-or-later
 ***************************************************************************/


template <class NT>
const char Polynomial<NT>::INDENT_SPACE[3] = { ' ', ' ', ' ' };  // pretty print parameters

// ==================================================
// Polynomial Constructors
// ==================================================

template <class NT>
Polynomial<NT>::Polynomial(void) {
  degree = -1;         // this is the zero polynomial!
  coeff = NULL;
}

//Creates a polynomial with nominal degree n
template <class NT>
Polynomial<NT>::Polynomial(int n) {
  CGAL_assertion(n>= -1);
  degree = n;
  if (n == -1)
    return;        // return the zero polynomial!
  if (n>=0)
    coeff = new NT[n+1];
  coeff[0]=1;                        // otherwise, return the unity polynomial
  for (int i=1; i<=n; i++)
    coeff[i]=0;
}

template <class NT>
Polynomial<NT>::Polynomial(int n, const NT * c) {
  //CGAL_assertion("array c has n+1 elements");
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
  degree = static_cast<int>(vN.size())-1;
  if (degree >= 0) {
    coeff = new NT[degree+1];
    for (int i = 0; i <= degree; i++)
      coeff[i] = vN[i];
  }
}

template <class NT>
Polynomial<NT>::Polynomial(const Polynomial<NT> & p):degree(-1) {
  //degree must be initialized to -1 otherwise delete is called on coeff in operator=
  coeff = NULL;//WHY?
  *this = p;        // reduce to assignment operator=
}

// Constructor from a Character string of coefficients
///////////////////////////////////////
template <class NT>
Polynomial<NT>::Polynomial(int n, const char * s[]) {
  //CGAL_assertion("array s has n+1 elements");
  degree = n;
  if (n >= 0) {
    coeff = new NT[n+1];
    for (int i = 0; i <= n; i++)
      coeff[i] = s[i];
  }
}

//The BNF syntax is the following:-
//    [poly] -> [term]| [term] '+/-' [poly] |
//                    '-' [term] | '-' [term] '+/-' [poly]
//    [term] -> [basic term] | [basic term] [term] | [basic term]*[term]
//    [basic term] -> [number] | 'x' | [basic term] '^' [number]
//                    | '(' [poly] ')'
//COMMENT:
//  [number] is assumed to be a BigInt; in the future, we probably
//  want to generalize this to BigFloat, etc.
//
template <class NT>
Polynomial<NT>::Polynomial(const std::string & s, char myX) {
   std::string ss(s);
   constructFromString(ss, myX);
}
template <class NT>
Polynomial<NT>::Polynomial(const char * s, char myX) {
   std::string ss(s);
   constructFromString(ss, myX);
}
template <class NT>
void Polynomial<NT>::constructFromString(std::string & s, char myX) {
  if(myX != 'x' || myX != 'X'){
    //Replace myX with 'x'.
    std::string::size_type loc = s.find(myX, 0);
    while(loc != std::string::npos){
      s.replace(loc,1,1,'x');
      loc = s.find(myX, loc+1);
    }
  }

  coeff = NULL;//Did this to ape the constructor from polynomial above
  *this = getpoly(s);
}

template <class NT>
Polynomial<NT>::~Polynomial() {
  if(degree >= 0)
    delete[] coeff;
}

  ////////////////////////////////////////////////////////
  // METHODS USED BY STRING CONSTRUCTOR ABOVE
  ////////////////////////////////////////////////////////


//Sets the input Polynomial to X^n
template <class NT>
void Polynomial<NT>::constructX(int n, Polynomial<NT>& P){
  Polynomial<NT> q(n);//Nominal degree n
  q.setCoeff(n,NT(1));
  if (n>0) q.setCoeff(0,NT(0));
  P = q;
}

//Returns in P the coeffecient starting from start
template <class NT>
int Polynomial<NT>::getnumber(const char* c, int start, unsigned int len,
                          Polynomial<NT> & P){
  int j=0;
  char *temp = new char[len];
  while(isint(c[j+start])){
    temp[j]=c[j+start];j++;
  }
  temp[j] = '\0';
  NT cf = NT(temp);
  Polynomial<NT> q(0);
  q.setCoeff(0, cf);
  P = q;
  delete[] temp;
  return (j-1+start);//Exactly the length of the number
}

template <class NT>
bool Polynomial<NT>::isint(char c){
  if(c == '0' || c == '1' || c == '2' || c == '3' || c == '4' ||
     c == '5' || c == '6' || c == '7' || c == '8' || c == '9')
    return true;
  else
    return false;
}


//Returns as integer the number starting from start in c
template <class NT>
int Polynomial<NT>::getint(const char* c, int start, unsigned int len,
                          int & n){
  int j=0;
  char *temp = new char[len];
  while(isint(c[j+start])){
    temp[j]=c[j+start];j++;
  }
  temp[j] = '\n';
  n = atoi(temp);
  delete[] temp;
  return (j-1+start);//Exactly the length of the number
}

//Given a string starting with an open parentheses returns the place
// which marks the end of the corresponding closing parentheses.
//Strings of the form (A).
template <class NT>
int Polynomial<NT>::matchparen(const char* cstr, int start){
  int count = 0;
  int j=start;

  do{
    if(cstr[j] == '('){
      count++;
    }
    if(cstr[j] == ')'){
      count--;
    }
    j++;
  }while(count != 0 );//j is one more than the matching ')'

  return j-1;
}


template <class NT>
int Polynomial<NT>::getbasicterm(std::string & s, Polynomial<NT> & P){
  const char * cstr = s.c_str();
  unsigned int len = s.length();
  int i=0;
  //Polynomial<NT> * temp = new Polynomial<NT>();

  if(isint(cstr[i])){
    i = getnumber(cstr, i, len, P);
  }else if(cstr[i] == 'x'||cstr[i] == 'X'){
    constructX(1, P);
  }else if(cstr[i] =='('){
    int oldi = i;
    i = matchparen(cstr, i);
    std::string t = s.substr(oldi+1, i -oldi -1);
    P = getpoly(t);
  }else{
#ifdef CGAL_CORE_TRACE
    std::cout <<"ERROR IN PARSING BASIC TERM" << std::endl;
#endif
  }
  //i+1 points to the beginning of next syntactic object in the string.
 if(cstr[i+1] == '^'){
    int n;
    i = getint(cstr, i+2, len, n);
    P.power(n);
  }
  return i;
}


template <class NT>
int Polynomial<NT>::getterm(std::string & s, Polynomial<NT> & P){
  unsigned int len = s.length();
  if(len == 0){// Zero Polynomial
    P=Polynomial<NT>();
    return 0;
  }
  unsigned int ind, oind;
  const char* cstr =s.c_str();
  std::string t;
  //P will be used to accumulate the product of basic terms.
  ind = getbasicterm(s, P);
  while(ind != len-1 && cstr[ind + 1]!='+' && cstr[ind + 1]!='-' ){
    //Recursively get the basic terms till we reach the end or see
    // a '+' or '-' sign.
    if(cstr[ind + 1] == '*'){
      t = s.substr(ind + 2, len - ind -2);
      oind = ind + 2;
    }else{
      t = s.substr(ind + 1, len -ind -1);
      oind = ind + 1;
    }

    Polynomial<NT> R;
    ind = oind + getbasicterm(t, R);//Because the second term is the offset in
                                     //t
    P *= R;
  }

  return ind;
}

template <class NT>
Polynomial<NT> Polynomial<NT>::getpoly(std::string & s){

    //Remove white spaces from the string
    std::string::size_type cnt=s.find(' ',0);
    while(cnt != std::string::npos){
      s.erase(cnt, 1);
      cnt = s.find(' ', cnt);
    }

    unsigned int len = s.length();
    if(len <= 0){//Zero Polynomial
      return Polynomial<NT>();
    }

    //To handle the case when there is one '=' sign
    //Suppose s is of the form s1 = s2. Then we assign s to
    //s1 + (-1)(s2) and reset len
    std::string::size_type loc;
    if((loc=s.find('=',0)) != std::string::npos){
      s.replace(loc,1,1,'+');
      std::string s3 = "(-1)(";
      s.insert(loc+1, s3);
      len = s.length();
      s.insert(len, 1, ')');
    }
    len = s.length();

    const char *cstr = s.c_str();
    std::string t;
    Polynomial<NT> P;
    // P will be the polynomial in which we accumulate the
    //sum and difference of the different terms.
    unsigned int ind;
    if(cstr[0] == '-'){
      t = s.substr(1, len);
      ind = getterm(t,P) + 1;
      P.negate();
    }else{
      ind = getterm(s, P);
    }
    unsigned int oind =0;//the string between oind and ind is a term
    while(ind != len -1){
      Polynomial<NT> R;
      t = s.substr(ind + 2, len -ind -2);
      oind = ind;
      ind = oind + 2 + getterm(t, R);
      if(cstr[oind + 1] == '+')
                P += R;
      else if(cstr[oind + 1] == '-')
                P -= R;
      else{
#ifdef CGAL_CORE_TRACE
        std::cout << "ERROR IN PARSING POLY! " << std::endl;
#endif
      }
    }

    return (P);
}

  ////////////////////////////////////////////////////////
  // METHODS
  ////////////////////////////////////////////////////////

// assignment:
template <class NT>
Polynomial<NT> & Polynomial<NT>::operator=(const Polynomial<NT>& p) {
  if (this == &p)
    return *this;        // self-assignment
  if (degree >=0)  delete[] coeff;
  degree = p.getDegree();
  if (degree < 0) return *this;
  coeff = new NT[degree +1];
  for (int i = 0; i <= degree; i++)
    coeff[i] = p.coeff[i];
  return *this;
}

// getTrueDegree
template <class NT>
int Polynomial<NT>::getTrueDegree() const {
  for (int i=degree; i>=0; i--) {
    if (sign(coeff[i]) != 0)
      return i;
  }
  return -1;        // Zero polynomial
}

//get i'th Coeff. We check whether i is not greater than the
// true degree and if not then return coeff[i] o/w 0
template <class NT>
NT Polynomial<NT>::getCoeffi(int i) const {
  int deg = getTrueDegree();
  if(i > deg)
      return NT(0);
  return coeff[i];
}

// ==================================================
// polynomial arithmetic
// ==================================================

// Expands the nominal degree to n;
//        Returns n if nominal degree is changed to n
//        Else returns -2
template <class NT>
int Polynomial<NT>::expand(int n) {
  if ((n <= degree)||(n < 0))
    return -2;
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
//        and returns the new (true) degree;
//        It returns -2 if this is a no-op
template <class NT>
int Polynomial<NT>::contract() {
  int d = getTrueDegree();
  if (d == degree)
    return (-2);  // nothing to do
  else
    degree = d;
  NT * c = coeff;
  if (degree !=-1){
    coeff = new NT[d+1];
    for (int i = 0; i<= d; i++)
      coeff[i] = c[i];
  }
  delete[] c;
  return d;
}

template <class NT>
Polynomial<NT> & Polynomial<NT>::operator+=(const Polynomial<NT>& p) { // +=
  int d = p.getDegree();
  if (d > degree)
    expand(d);
  for (int i = 0; i<=d; i++)
    coeff[i] += p.coeff[i];

  return *this;
}
template <class NT>
Polynomial<NT> & Polynomial<NT>::operator-=(const Polynomial<NT>& p) { // -=
  int d = p.getDegree();
  if (d > degree)
    expand(d);
  for (int i = 0; i<=d; i++)
    coeff[i] -= p.coeff[i];
  return *this;
}

// SELF-MULTIPLICATION
// This is quadratic time multiplication!
template <class NT>
Polynomial<NT> & Polynomial<NT>::operator*=(const Polynomial<NT>& p) { // *=
  if (degree==-1) return *this;
  if (p.getDegree()==-1){
    degree=-1;
    delete[] coeff;
    return *this;
  }
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
  if (s==0)
    return *this;
  int d = s+getTrueDegree();
  if (d < 0) {
    degree = -1;
    delete[] coeff;
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
      c[d-j] = coeff[d-s-j];  // since s<0, (d-s-j) > (d-j) !!
  }
  delete[] coeff;
  coeff = c;
  degree = d;
  return *this;
}//mulXpower

// REDUCE STEP (helper for PSEUDO-REMAINDER function)
// Let THIS=(*this) be the current polynomial, and P be the input
//        argument for reduceStep.  Let R be returned polynomial.
//        R has the special form as a binomial,
//                R = C + X*M
//        where C is a constant and M a monomial (= coeff * some power of X).
//        Moreover, THIS is transformed to a new polynomial, THAT, which
//        is given by
//                 (C * THIS) = M * P  + THAT
//        MOREOVER: deg(THAT) < deg(THIS) unless deg(P)>deg(Q).
//        Basically, C is a power of the leading coefficient of P.
//        REMARK: R is NOT defined as C+M, because in case M is
//        a constant, then we cannot separate C from M.
//        Furthermore, R.mulScalar(-1) gives us M.
template <class NT>
Polynomial<NT> Polynomial<NT>::reduceStep (
  const Polynomial<NT>& p) {
  //         Chee: Omit the next 2 contractions as unnecessary
  //         since reduceStep() is only called by pseudoRemainder().
  //         Also, reduceStep() already does a contraction before returning.
  // p.contract();
  // contract();        // first contract both polynomials
  Polynomial<NT> q(p);                // q is initially a copy of p
  //        but is eventually M*P
  int pDeg  = q.degree;
  int myDeg = degree;
  if (pDeg == -1)
    return *(new Polynomial());  // Zero Polynomial
  // NOTE: pDeg=-1 (i.e., p=0) is really an error condition!
  if (myDeg < pDeg)
    return *(new Polynomial(0));  // Unity Polynomial
  // i.e., C=1, M=0.
  // Now (myDeg >= pDeg).  Start to form the Return Polynomial R=C+X*M
  Polynomial<NT> R(myDeg - pDeg + 1);  // deg(M)= myDeg - pDeg
  q.mulXpower(myDeg - pDeg);           // q is going to become M*P

  NT myLC = coeff[myDeg];          // LC means "leading coefficient"
  NT qLC = q.coeff[myDeg];  // p also has degree "myDeg" (qLC non-zero)
  NT LC;

  //  NT must support
  //  isDivisible(x,y), gcd(x,y), div_exact(x,y) in the following:
  //  ============================================================
  if (isDivisible(myLC, qLC)) { // myLC is divisible by qLC
    LC = div_exact(myLC, qLC);
    R.setCoeff(0, 1);                   //  C = 1,

    R.setCoeff(R.degree, LC); //  M = LC * X^(myDeg-pDeg)
    q.mulScalar(LC);           //  q = M*P.
  }
  else if (isDivisible(qLC, myLC)) { // qLC is divisible by myLC
    LC = div_exact(qLC, myLC);         //
    if ((LC != 1) && (LC != -1)) { // IMPORTANT: unlike the previous
      // case, we need an additional condition
      // that LC != -1.  THis is because
      // if (LC = -1), then we have qLC and
      // myLC are mutually divisible, and
      // we would be updating R twice!
      R.setCoeff(0, LC);            // C = LC,
      R.setCoeff(R.degree, 1);     // M = X^(myDeg-pDeg)(THIS WAS NOT DONE)
      mulScalar(LC);                       // THIS => THIS * LC

    }
  } else {                          // myLC and qLC are not mutually divisible
    NT g = gcd(qLC, myLC);         // This ASSUMES gcd is defined in NT !!
    //NT g = 1;                          // In case no gcd is available
    if (g == 1) {
      R.setCoeff(0, qLC);                  // C = qLC
      R.setCoeff(R.degree, myLC);         // M = (myLC) * X^{myDeg-pDeg}
      mulScalar(qLC);                         // forming  C * THIS
      q.mulScalar(myLC);                // forming  M * P
    } else {
      NT qLC1= div_exact(qLC,g);
      NT myLC1= div_exact(myLC,g);
      R.setCoeff(0, qLC1);                  // C = qLC/g
      R.setCoeff(R.degree, myLC1);        // M = (myLC/g) * X^{myDeg-pDeg}
      mulScalar(qLC1);                 // forming  C * THIS
      q.mulScalar(myLC1);                // forming  M * P
    }
  }
  (*this) -= q;                // THAT = (C * THIS) - (M * P)

  contract();

  return R;                // Returns R = C + X*M
}// reduceStep

// For internal use only:
// Checks that c*A = B*m + AA
//         where A=(*oldthis) and AA=(*newthis)
template <class NT>
Polynomial<NT> Polynomial<NT>::testReduceStep(const Polynomial<NT>& A,
        const Polynomial<NT>& B) {
std::cout << "+++++++++++++++++++++TEST REDUCE STEP+++++++++++++++++++++\n";
  Polynomial<NT> cA(A);
  Polynomial<NT> AA(A);
  Polynomial<NT> quo;
  quo = AA.reduceStep(B);                // quo = c + X*m  (m is monomial, c const)
                                // where c*A = B*m + (*newthis)
std::cout << "A = " << A << std::endl;
std::cout << "B = " << B << std::endl;
  cA.mulScalar(quo.coeff[0]);    // A -> c*A
  Polynomial<NT> m(quo);
  m.mulXpower(-1);            // m's value is now m
std::cout << "c = " << quo.coeff[0] << std::endl;
std::cout << "c + xm = " << quo << std::endl;
std::cout << "c*A = " << cA << std::endl;
std::cout << "AA = " << AA << std::endl;
std::cout << "B*m = " << B*m << std::endl;
std::cout << "B*m + AA = " << B*m + AA << std::endl;
  if (cA == (B*m + AA))
          std::cout << "CORRECT inside testReduceStep" << std::endl;
  else
          std::cout << "ERROR inside testReduceStep" << std::endl;
std::cout << "+++++++++++++++++END TEST REDUCE STEP+++++++++++++++++++++\n";
  return quo;
}

// PSEUDO-REMAINDER and PSEUDO-QUOTIENT:
// Let A = (*this) and B be the argument polynomial.
// Let Quo be the returned polynomial,
// and let the final value of (*this) be Rem.
// Also, C is the constant that we maintain.
// We are computing A divided by B.  The relation we guarantee is
//                 (C * A) = (Quo * B)  + Rem
// where deg(Rem) < deg(B).  So Rem is the Pseudo-Remainder
// and Quo is the Pseudo-Quotient.
// Moreover, C is uniquely determined (we won't spell it out)
// except to note that
//        C divides D = (LC)^{deg(A)-deg(B)+1}
//        where LC is the leading coefficient of B.
// NOTE: 1. Normally, Pseudo-Remainder is defined so that C = D
//          So be careful when using our algorithm.
//          2. We provide a version of pseudoRemainder which does not
//          require C as argument.  [For efficiency, we should provide this
//         version directly, instead of calling the other version!]

template <class NT>
Polynomial<NT> Polynomial<NT>::pseudoRemainder (
  const Polynomial<NT>& B) {
        NT temp;        // dummy argument to be discarded
        return pseudoRemainder(B, temp);
}//pseudoRemainder

template <class NT>
Polynomial<NT> Polynomial<NT>::pseudoRemainder (
  const Polynomial<NT>& B, NT & C) {
  contract();         // Let A = (*this).  Contract A.
  Polynomial<NT> tmpB(B);
  tmpB.contract();    // local copy of B
  int bTrueDegree = tmpB.degree;
  C = NT(1);  // Initialized to C=1.
  if (bTrueDegree == -1)  {
    core_error("ERROR in Polynomial<NT>::pseudoRemainder :\n    -- divide by zero polynomial", __FILE__, __LINE__, false);
    return Polynomial(0);  // Unit Polynomial (arbitrary!)
  }
  if (bTrueDegree > degree) {
    return Polynomial(); // Zero Polynomial
    // CHECK: 1*THIS = 0*B + THAT,  deg(THAT) < deg(B)
  }

  Polynomial<NT> Quo;  // accumulate the return polynomial, Quo
  Polynomial<NT> tmpQuo;
  while (degree >= bTrueDegree) {  // INVARIANT: C*A = B*Quo + (*this)
    tmpQuo = reduceStep(tmpB);  // Let (*this) be (*oldthis), which
                                // is transformed into (*newthis). Then,
                                //     c*(*oldthis) = B*m + (*newthis)
                                // where tmpQuo = c + X*m
    // Hence,   C*A =   B*Quo +   (*oldthis)      -- the old invariant
    //        c*C*A = c*B*Quo + c*(*oldthis)
    //              = c*B*Quo + (B*m + (*newthis))
    //              = B*(c*Quo + m)  + (*newthis)
    // i.e, to update invariant, we do C->c*C,  Quo --> c*Quo + m.
    C *= tmpQuo.coeff[0];            // C = c*C
    Quo.mulScalar(tmpQuo.coeff[0]); // Quo -> c*Quo
    tmpQuo.mulXpower(-1);           // tmpQuo is now equal to m
    Quo += tmpQuo;                  // Quo -> Quo + m
  }

  return Quo;        // Quo is the pseudo-quotient
}//pseudoRemainder

// Returns the negative of the pseudo-remainder
//         (self-modification)
template <class NT>
Polynomial<NT> & Polynomial<NT>::negPseudoRemainder (
  const Polynomial<NT>& B) {
        NT temp;        // dummy argument to be discarded
        pseudoRemainder(B, temp);
        if (temp < 0) return (*this);
        return negate();
}

template <class NT>
Polynomial<NT> & Polynomial<NT>::operator-() {        // unary minus
  for (int i=0; i<=degree; i++)
    coeff[i] *= -1;
  return *this;
}
template <class NT>
Polynomial<NT> & Polynomial<NT>::power(unsigned int n) {        // self-power
  if (n == 0) {
    degree = 0;
    delete [] coeff;
    coeff = new NT[1];
    coeff[0] = 1;
  } else {
    Polynomial<NT> p = *this;
    for (unsigned int i=1; i<n; i++)
      *this *= p;                // Would a binary power algorithm be better?
  }
  return *this;
}

// evaluation of BigFloat value
//   NOTE: we think of the polynomial as an analytic function in this setting
//      If the coefficients are more general than BigFloats,
//      then we may get inexact outputs, EVEN if the input value f is exact.
//      This is because we would have to convert these coefficients into
//      BigFloats, and this conversion precision is controlled by the
//      global variables defRelPrec and defAbsPrec.
//
/*
template <class NT>
BigFloat Polynomial<NT>::eval(const BigFloat& f) const {        // evaluation
  if (degree == -1)
    return BigFloat(0);
  if (degree == 0)
    return BigFloat(coeff[0]);
  BigFloat val(0);
  for (int i=degree; i>=0; i--) {
    val *= f;
    val += BigFloat(coeff[i]);
  }
  return val;
}//eval
*/

/// Evaluation Function (generic version, always returns the exact value)
///
///  This evaluation function is easy to use, but may not be efficient
///  when you have BigRat or Expr values.
///
/// User must be aware that the return type of eval is Max of Types NT and T.
///
/// E.g., If NT is BigRat, and T is Expr then Max(NT,T)=Expr.
///
/// REMARK: If NT is BigFloat, it is assumed that the BigFloat is error-free.

template <class NT>
template <class T>
MAX_TYPE(NT, T) Polynomial<NT>::eval(const T& f) const {        // evaluation
  typedef MAX_TYPE(NT, T) ResultT;
  if (degree == -1)
    return ResultT(0);
  if (degree == 0)
    return ResultT(coeff[0]);
  ResultT val(0);
  for (int i=degree; i>=0; i--) {
    val *= ResultT(f);
    val += ResultT(coeff[i]);
  }
  return val;
}//eval


/// Approximate Evaluation of Polynomials
///         the coefficients of the polynomial are approximated to some
///        specified composite precision (r,a).
/// @param  f evaluation point
/// @param  r relative precision to which the coefficients are evaluated
/// @param  a absolute precision to which the coefficients are evaluated
/// @return a BigFloat with error containing value of the polynomial.
///     If zero is in this BigFloat interval, then we don't know the sign.
//
//         ASSERT: NT = BigRat or Expr
//
template <class NT>
BigFloat Polynomial<NT>::evalApprox(const BigFloat& f,
        const extLong& r, const extLong& a) const {        // evaluation
  if (degree == -1)
    return BigFloat(0);
  if (degree == 0)
    return BigFloat(coeff[0], r, a);

  BigFloat val(0), c;
  for (int i=degree; i>=0; i--) {
    c = BigFloat(coeff[i], r, a);
    val *= f;
    val += c;
  }
  return val;
}//evalApprox

// This BigInt version of evalApprox should never be called...
template <>
CORE_INLINE
BigFloat Polynomial<BigInt>::evalApprox(const BigFloat& /*f*/,
        const extLong& /*r*/, const extLong& /*a*/) const {        // evaluation
  CGAL_assertion(0);
  return BigFloat(0);
}



/**
 * Evaluation at a BigFloat value
 * using "filter" only when NT is BigRat or Expr.
 * Using filter means we approximate the polynomials
 * coefficients using BigFloats.  If this does not give us
 * the correct sign, we will resort to an "exact" evaluation
 * using Expr.
 *
 * If NT <= BigFloat, we just call eval().
 *
   We use the following heuristic estimates of precision for coefficients:

      r = 1 + lg(|P|_\infty) + lg(d+1)                  if f <= 1
      r = 1 + lg(|P|_\infty) + lg(d+1) + d*lg|f|         if f > 1

   if the filter fails, then we use Expr to do evaluation.

   This function is mainly called by Newton iteration (which
   has some estimate for oldMSB from previous iteration).

   @param p polynomial to be evaluated
   @param val the evaluation point
   @param oldMSB an rough estimate of the lg|p(val)|
   @return bigFloat interval contain p(val), with the correct sign

 ***************************************************/
template <class NT>
BigFloat Polynomial<NT>::evalExactSign(const BigFloat& val,
         const extLong& oldMSB) const {
    CGAL_assertion(val.isExact());
    if (getTrueDegree() == -1)
      return BigFloat(0);

    extLong r;
    r = 1 + BigFloat(height()).uMSB() + clLg(long(getTrueDegree()+1));
    if (val > 1)
      r += getTrueDegree() * val.uMSB();
    r += core_max(extLong(0), -oldMSB);

    if (hasExactDivision<NT>::check()) { // currently, only to detect NT=Expr and NT=BigRat
        BigFloat rVal = evalApprox(val, r);
        if (rVal.isZeroIn()) {
          Expr eVal = eval(Expr(val));        // eval gives exact value
          eVal.approx(54,CORE_INFTY);  // if value is 0, we get exact 0
          return eVal.BigFloatValue();
        } else
          return rVal;
    } else
        return BigFloat(eval(val));

   //return 0; // unreachable
  }//evalExactSign


//============================================================
// Bounds
//============================================================

// Cauchy Upper Bound on Roots
// -- ASSERTION: NT is an integer type
template < class NT >
BigFloat Polynomial<NT>::CauchyUpperBound() const {
  if (zeroP(*this))
    return BigFloat(0);
  NT mx = 0;
  int deg = getTrueDegree();
  for (int i = 0; i < deg; ++i) {
    mx = core_max(mx, abs(coeff[i]));
  }
  Expr e = mx;
  e /= Expr(abs(coeff[deg]));
  e.approx(CORE_INFTY, 2);
  // get an absolute approximate value with error < 1/4
  return (e.BigFloatValue().makeExact() + 2);
}

//============================================================
// An iterative version of computing Cauchy Bound from Erich Kaltofen.
// See the writeup under collab/sep/.
//============================================================
template < class NT >
BigInt Polynomial<NT>::CauchyBound() const {
  int deg = getTrueDegree();
  BigInt B(1);
  BigFloat lhs(0), rhs(1);
  while (true) {
    /* compute \sum_{i=0}{deg-1}{|a_i|B^i} */
    lhs = 0;
    for (int i=deg-1; i>=0; i--) {
      lhs *= B;
      lhs += abs(coeff[i]);
    }
    lhs /= abs(coeff[deg]);
    lhs.makeFloorExact();
    /* compute B^{deg} */
    if (rhs <= lhs) {
      B <<= 1;
      rhs *= (BigInt(1)<<deg);
    } else
      break;
  }
  return B;
}

//====================================================================
//Another iterative bound which is at least as good as the above bound
//by Erich Kaltofen.
//====================================================================
template < class NT >
BigInt Polynomial<NT>::UpperBound() const {
  int deg = getTrueDegree();

  BigInt B(1);
  BigFloat lhsPos(0), lhsNeg(0), rhs(1);
  while (true) {
    /* compute \sum_{i=0}{deg-1}{|a_i|B^i} */
    lhsPos = lhsNeg = 0;
    for (int i=deg-1; i>=0; i--) {
      if (coeff[i]>0) {
              lhsPos = lhsPos * B + coeff[i];
              lhsNeg = lhsNeg * B;
      } else {
              lhsNeg = lhsNeg * B - coeff[i];
              lhsPos = lhsPos * B;
      }
    }
    lhsNeg /= abs(coeff[deg]);
    lhsPos /= abs(coeff[deg]);
    lhsPos.makeCeilExact();
    lhsNeg.makeCeilExact();

    /* compute B^{deg} */
    if (rhs <= (std::max)(lhsPos,lhsNeg)) {
      B <<= 1;
      rhs *= (BigInt(1)<<deg);
    } else
      break;
  }
  return B;
}

// Cauchy Lower Bound on Roots
// -- ASSERTION: NT is an integer type
template < class NT >
BigFloat Polynomial<NT>::CauchyLowerBound() const {
  if ((zeroP(*this)) || coeff[0] == 0)
    return BigFloat(0);
  NT mx = 0;
  int deg = getTrueDegree();
  for (int i = 1; i <= deg; ++i) {
    mx = core_max(mx, abs(coeff[i]));
  }
  Expr e = Expr(abs(coeff[0]))/ Expr(abs(coeff[0]) + mx);
  e.approx(2, CORE_INFTY);
  // get an relative approximate value with error < 1/4
  return (e.BigFloatValue().makeExact().div2());
}

// Separation bound for polynomials that may have multiple roots.
// We use the Rump-Schwartz bound.
//
//    ASSERT(the return value is an exact BigFloat and a Lower Bound)
//
template < class NT >
BigFloat Polynomial<NT>::sepBound() const {
  BigInt d;
  BigFloat e;
  int deg = getTrueDegree();

  CORE::power(d, BigInt(deg), ((deg)+4)/2);
  e = CORE::power(height()+1, deg);
  e.makeCeilExact(); // see NOTE below
  return (1/(e*2*d)).makeFloorExact();
        // BUG fix: ``return 1/e*2*d'' was wrong
        // NOTE: the relative error in this division (1/(e*2*d))
        //   is defBFdivRelPrec (a global variable), but
        //   since this is always positive, we are OK.
        //   But to ensure that defBFdivRelPrec is used,
        //   we must make sure that e and d are exact.
        // Finally, using "makeFloorExact()" is OK because
        //   the mantissa minus error (i.e., m-err) will remain positive
        //   as long as the relative error (defBFdivRelPrec) is >1.
}

/// height function
/// @return a BigFloat with error
template < class NT >
BigFloat Polynomial<NT>::height() const {
  if (zeroP(*this))
    return BigFloat(0);
  int deg = getTrueDegree();
  NT ht = 0;
  for (int i = 0; i< deg; i++)
    if (ht < abs(coeff[i]))
      ht = abs(coeff[i]);
  return BigFloat(ht);
}

/// length function
/// @return a BigFloat with error
template < class NT >
BigFloat Polynomial<NT>::length() const {
  if (zeroP(*this))
    return BigFloat(0);
  int deg = getTrueDegree();
  NT length = 0;
  for (int i = 0; i< deg; i++)
    length += abs(coeff[i]*coeff[i]);
  return sqrt(BigFloat(length));
}

//============================================================
// differentiation
//============================================================

template <class NT>
Polynomial<NT> & Polynomial<NT>::differentiate() {        // self-differentiation
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
  CGAL_assertion(n >= 0);
  for (int i=1; i<=n; i++)
    this->differentiate();
  return *this;
} // multi-differentiate


// ==================================================
// GCD, content, primitive and square-free parts
// ==================================================

/// divisibility predicate for polynomials
// isDivisible(P,Q) returns true iff Q divides P
// To FIX: the predicate name is consistent with Expr.h but not with BigInt.h
template <class NT>
bool isDivisible(Polynomial<NT> p, Polynomial<NT> q) {
  if(zeroP(p))
    return true;
  if(zeroP(q))
    return false;  // We should really return error!
  if(p.getTrueDegree() < q.getTrueDegree())
    return false;
  p.pseudoRemainder(q);
  if(zeroP(p))
    return true;
  else
    return false;
}//isDivisible

//Content of a polynomial P
//      -- content(P) is just the gcd of all the coefficients
//      -- REMARK: by definition, content(P) is non-negative
//                 We rely on the fact that NT::gcd always
//                 return a non-negative value!!!
template <class NT>
NT content(const Polynomial<NT>& p) {
  if(zeroP(p))
    return 0;
  int d = p.getTrueDegree();
  if(d == 0){
    if(p.coeff[0] > 0)
      return p.coeff[0];
    else
      return -p.coeff[0];
  }

  NT content = p.coeff[d];
  for (int i=d-1; i>=0; i--) {
    content = gcd(content, p.coeff[i]);
    if(content == 1) break;   // remark: gcd is non-negative, by definition
  }
  //if (p.coeff[0] < 0) return -content;(BUG!)
  return content;
}//content

// Primitive Part:  (*this) is transformed to primPart and returned
//        -- primPart(P) is just P/content(P)
//        -- Should we return content(P) instead? [SHOULD IMPLEMENT THIS]
// IMPORTANT: we require that content(P)>0, hence
//         the coefficients of primPart(P) does
//         not change sign; this is vital for use in Sturm sequences
template <class NT>
Polynomial<NT> & Polynomial<NT>::primPart() {
  // ASSERT: GCD must be provided by NT
  int d = getTrueDegree();
  CGAL_assertion (d >= 0);
  if (d == 0) {
    if (coeff[0] > 0) coeff[0] = 1;
    else coeff[0] = -1;
    return *this;
  }

  NT g = content(*this);
  if (g == 1 && coeff[d] > 0)
     return (*this);
  for (int i=0; i<=d; i++) {
     coeff[i] =  div_exact(coeff[i], g);
  }
  return *this;
}// primPart

//GCD of two polynomials.
//   --Assumes that the coeffecient ring NT has a gcd function
//   --Returns the gcd with a positive leading coefficient (*)
//     otherwise division by gcd causes a change of sign affecting
//     Sturm sequences.
//   --To Check: would a non-recursive version be much better?
template <class NT>
Polynomial<NT> gcd(const Polynomial<NT>& p, const Polynomial<NT>& q) {

  // If the first polynomial has a smaller degree then the second,
  // then change the order of calling
  if(p.getTrueDegree() < q.getTrueDegree())
    return gcd(q,p);

  // If any polynomial is zero then the gcd is the other polynomial
  if(zeroP(q)){
    if(zeroP(p))
       return p;
    else{
       if(p.getCoeffi(p.getTrueDegree()) < 0){
         return Polynomial<NT>(p).negate();
       }else
         return p;        // If q<>0, then we know p<>0
   }
  }
  Polynomial<NT> temp0(p);
  Polynomial<NT> temp1(q);

  // We want to compute:
  // gcd(p,q) = gcd(content(p),content(q)) * gcd(primPart(p), primPart(q))

  NT cont0 = content(p);        // why is this temporary needed?
  NT cont1 = content(q);
  NT cont = gcd(cont0,cont1);
  temp0.primPart();
  temp1.primPart();

  temp0.pseudoRemainder(temp1);
  return (gcd(temp1, temp0).mulScalar(cont));
}//gcd

// sqFreePart()
//         -- this is a self-modifying operator!
//         -- Let P =(*this) and Q=square-free part of P.
//         -- (*this) is transformed into P, and gcd(P,P') is returned
// NOTE: The square-free part of P is defined to be P/gcd(P, P')
template <class NT>
Polynomial<NT>  Polynomial<NT>::sqFreePart() {

  int d = getTrueDegree();
  if(d <= 1) // linear polynomials and constants are square-free
    return *this;

  Polynomial<NT> temp(*this);
  Polynomial<NT> R = gcd(*this, temp.differentiate()); // R = gcd(P, P')

  // If P and P' have a constant gcd, then P is squarefree
  if(R.getTrueDegree() == 0)
    return (Polynomial<NT>(0)); // returns the unit polynomial as gcd

  (*this)=pseudoRemainder(R); // (*this) is transformed to P/R, the sqfree part
  //Note: This is up to multiplication by units
  return (R); // return the gcd
}//sqFreePart()


// ==================================================
// Useful member functions
// ==================================================

// reverse:
//         reverses the list of coefficients
template <class NT>
void Polynomial<NT>::reverse() {
  NT tmp;
  for (int i=0; i<= degree/2; i++) {
    tmp = coeff[i];
    coeff[i] =   coeff[degree-i];
    coeff[degree-i] = tmp;
  }
}//reverse

// negate:
//         multiplies the polynomial by -1
//         Chee: 4/29/04 -- added negate() to support negPseudoRemainder(B)
template <class NT>
Polynomial<NT> & Polynomial<NT>::negate() {
  for (int i=0; i<= degree; i++)
    coeff[i] *= -1;          // every NT must be able to construct from -1
  return *this;
}//negate

// makeTailCoeffNonzero
//         Divide (*this) by X^k, so that the tail coeff becomes non-zero.
//        The value k is returned.  In case (*this) is 0, return k=-1.
//        Otherwise, if (*this) is unchanged, return k=0.
template <class NT>
int Polynomial<NT>::makeTailCoeffNonzero() {
  int k=-1;
  for (int i=0; i <= degree; i++) {
    if (coeff[i] != 0) {
      k=i;
      break;
    }
  }
  if (k <= 0)
    return k;        // return either 0 or -1
  degree -=k;                // new (lowered) degree
  NT * c = new NT[1+degree];
  for (int i=0; i<=degree; i++)
    c[i] = coeff[i+k];
  delete[] coeff;
  coeff = c;
  return k;
}//

// filedump(string msg, ostream os, string com, string com2):
//      Dumps polynomial to output stream os
//      msg is any message
//      NOTE: Default is com="#", which is placed at start of each
//            output line.
template <class NT>
void Polynomial<NT>::filedump(std::ostream & os,
                          std::string msg,
                          std::string commentString,
                          std::string commentString2) const {
  int d= getTrueDegree();
  if (msg != "") os << commentString << msg << std::endl;
  int i=0;
  if (d == -1) { // if zero polynomial
    os << commentString << "0";
    return;
  }
  for (; i<=d;  ++i)        // get to the first non-zero coeff
    if (coeff[i] != 0)
      break;
  int termsInLine = 1;

  // OUTPUT the first nonzero term
  os << commentString;
  if (coeff[i] == 1) {                        // special cases when coeff[i] is
    if (i>1) os << "x^" << i;        // either 1 or -1
    else if (i==1) os << "x" ;
    else os << "1";
  } else if (coeff[i] == -1) {
    if (i>1) os << "-x^" << i;
    else if (i==1) os << "-x" ;
    else os << "-1";
  } else {                                // non-zero coeff
    os << coeff[i];
    if (i>1) os << "*x^" << i;
    else if (i==1) os << "x" ;
  }
  // OUTPUT the remaining nonzero terms
  for (i++ ; i<= getTrueDegree(); ++i) {
    if (coeff[i] == 0)
      continue;
    termsInLine++;
    if (termsInLine % Polynomial<NT>::COEFF_PER_LINE == 0) {
      os << std::endl;
      os << commentString2;
    }
    if (coeff[i] == 1) {                // special when coeff[i] = 1
      if (i==1) os << " + x";
      else os << " + x^" << i;
    } else if (coeff[i] == -1) {        // special when coeff[i] = -1
      if (i==1) os << " - x";
      else os << " -x^" << i;
    } else {                                // general case
      if(coeff[i] > 0){
        os << " + ";
        os << coeff[i];
      }else
        os << coeff[i];

      if (i==1) os << "*x";
      else os << "*x^" << i;
    }
  }
}//filedump

// dump(message, ofstream, commentString) -- dump to file
template <class NT>
void Polynomial<NT>::dump(std::ofstream & ofs,
                std::string msg,
                std::string commentString,
                std::string commentString2) const {
  filedump(ofs, msg, commentString, commentString2);
}

// dump(message)         -- to std output
template <class NT>
void Polynomial<NT>::dump(std::string msg, std::string com,
                std::string com2) const {
  filedump(std::cout, msg, com, com2);
}

// Dump of Maple Code for Polynomial
template <class NT>
void Polynomial<NT>::mapleDump() const {
  if (zeroP(*this)) {
    std::cout << 0 << std::endl;
    return;
  }
  std::cout << coeff[0];
  for (int i = 1; i<= getTrueDegree(); ++i) {
    std::cout << " + (" << coeff[i] << ")";
    std::cout << "*x^" << i;
  }
  std::cout << std::endl;
}//mapleDump

// ==================================================
// Useful friend functions for Polynomial<NT> class
// ==================================================

// friend differentiation
template <class NT>
Polynomial<NT> differentiate(const Polynomial<NT> & p) {          // differentiate
  Polynomial<NT> q(p);
  return q.differentiate();
}

// friend multi-differentiation
template <class NT>
Polynomial<NT> differentiate(const Polynomial<NT> & p, int n) {//multi-differentiate
  Polynomial<NT> q(p);
  CGAL_assertion(n >= 0);
  for (int i=1; i<=n; i++)
    q.differentiate();
  return q;
}

// friend equality comparison
template <class NT>
bool operator==(const Polynomial<NT>& p, const Polynomial<NT>& q) {        // ==
  int d, D;
  Polynomial<NT> P(p);
  P.contract();
  Polynomial<NT> Q(q);
  Q.contract();
  if (P.degree < Q.degree) {
    d = P.degree;
    D = Q.degree;
    for (int i = d+1; i<=D; i++)
      if (Q.coeff[i] != 0)
        return false;        // return false
  } else {
    D = P.degree;
    d = Q.degree;
    for (int i = d+1; i<=D; i++)
      if (P.coeff[i] != 0)
        return false;        // return false
  }
  for (int i = 0; i <= d; i++)
    if (P.coeff[i] != Q.coeff[i])
      return false;        // return false
  return true;         // return true
}

// friend non-equality comparison
template <class NT>
bool operator!=(const Polynomial<NT>& p, const Polynomial<NT>& q) {        // !=
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
  // Don't you need to first do  "delete[] p.coeff;"  ??
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


//Resultant of two polynomials.
//Since we use pseudoRemainder we have to modify the original algorithm.
//If C * P = Q * R + S, where C is a constant and S = prem(P,Q), m=deg(P),
// n=deg(Q) and l = deg(S).
//Then res(P,Q) = (-1)^(mn) b^(m-l) res(Q,S)/C^(n)
//The base case being res(P,Q) = Q^(deg(P)) if Q is a constant, zero otherwise
template <class NT>
NT res(Polynomial<NT> p, Polynomial<NT> q) {

  int m, n;
  m = p.getTrueDegree();
  n = q.getTrueDegree();

  if(m == -1 || n == -1) return 0;  // this definition is not certified
  if(m == 0 && n == 0) return 1;    // this definition is at variance from
                                    // Yap's book (but is OK for purposes of
                                    // checking the vanishing of resultants
  if(n > m) return (res(q, p));

  NT b = q.coeff[n];//Leading coefficient of Q
  NT lc = p.coeff[m], C;

  p.pseudoRemainder(q, C);

  if(zeroP(p) && n ==0)
     return (pow(q.coeff[0], m));

  int l = p.getTrueDegree();

  return(pow(NT(-1), m*n)*pow(b,m-l)*res(q,p)/pow(C,n));
}

//i^th Principal Subresultant Coefficient (psc) of two polynomials.
template <class NT>
NT psc(int i,Polynomial<NT> p, Polynomial<NT> q) {

  CGAL_assertion(i >= 0);
  if(i == 0)
     return res(p,q);

  int m = p.getTrueDegree();
  int n = q.getTrueDegree();

  if(m == -1 || n == -1) return 0;
  if(m < n) return psc(i, q, p);

  if ( i == n) //This case occurs when i > degree of gcd
     return pow(q.coeff[n], m - n);

  if(n < i && i <= m)
     return 0;

  NT b = q.coeff[n];//Leading coefficient of Q
  NT lc = p.coeff[m], C;


  p.pseudoRemainder(q, C);

  if(zeroP(p) && i < n)//This case occurs when i < deg(gcd)
     return 0;

  if(zeroP(p) && i == n)//This case occurs when i=deg(gcd) might be constant
     return pow(q.coeff[n], m - n);

  int l = p.getTrueDegree();
  return pow(NT(-1),(m-i)*(n-i))*pow(b,m-l)*psc(i,q,p)/pow(C, n-i);

}

//factorI(p,m)
//    computes the polynomial q which containing all roots
//    of multiplicity m of a given polynomial P
//    Used to determine the nature of intersection of two algebraic curves
//    The algorithm is given in Wolperts Thesis
template <class NT>
Polynomial<NT> factorI(Polynomial<NT> p, int m){
  int d=p.getTrueDegree();
  Polynomial<NT> *w = new Polynomial<NT>[d+1];
  w[0] = p;
  Polynomial<NT> temp;

  for(int i = 1; i <=d ; i ++){
     temp = w[i-1];
     w[i] = gcd(w[i-1],temp.differentiate());
  }

  Polynomial<NT> *u = new Polynomial<NT>[d+1];
  u[d] = w[d-1];
  for(int i = d-1; i >=m; i--){
        temp = power(u[i+1],2);
     for(int j=i+2; j <=d; j++){
        temp *= power(u[j],j-i+1);
     }
     u[i] = w[i-1].pseudoRemainder(temp);//i-1 since the array starts at 0
  }

        delete[] w;
  return u[m];
}

// ==================================================
// End of Polynomial<NT>
// ==================================================
