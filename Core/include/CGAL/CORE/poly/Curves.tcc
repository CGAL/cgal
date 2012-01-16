/****************************************************************************
 * Core Library Version 1.7, August 2004
 * Copyright (c) 1995-2004 Exact Computation Project
 * All rights reserved.
 *
 * This file is part of CORE (http://cs.nyu.edu/exact/core/).
 * You can redistribute it and/or modify it under the terms of the GNU
 * General Public License as published by the Free Software Foundation,
 * either version 3 of the License, or (at your option) any later version.
 *
 * Licensees holding a valid commercial license may use this file in
 * accordance with the commercial license agreement provided with the
 * software.
 *
 * This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
 * WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 *
 *
 *  File: Curves.tcc
 *
 *  Description: 
 *  	This file contains the implementations of
 *		functions defined by Curves.h
 *	It is included through Curves.h
 *  	Please see Curves.h for more details.
 *
 *  Author:  Vikram Sharma and Chee Yap
 *  Date:    April 12, 2004
 *
 * WWW URL: http://cs.nyu.edu/exact/
 * Email: exact@cs.nyu.edu
 *
 * $URL$
 * $Id$
 ***************************************************************************/


//CONSTRUCTORS FOR THE BIPOLY CLASS
  ////////////////////////////////////////////////////////
  //Constructors
  ////////////////////////////////////////////////////////

template <class NT>
BiPoly<NT>::BiPoly(){ // zero polynomial
    ydeg = -1;
  }
  
  //BiPoly(n)

template <class NT>
BiPoly<NT>::BiPoly(int n){// creates a BiPoly with nominal y-degree equal to n.

     // To support this constructor, you need functions
     // that are equivalent to "setCoeff(i, PolyNT q)" in the Poly Class
     // You also want "getCoeff(i)" functions.
     // E.g. BiPoly<NT> p(10);
     //      PolyNT q(3, {1,2,3,4});
     //      p.setCoeff(10, q);
     //      p.setCoeff(0, PolyNT());
     //
    ydeg = n;
    if (n<0) return; // coeffX not needed
    for(int i=0; i<= ydeg; i++){
      Polynomial<NT> temp;
      coeffX.push_back(temp);      
    }
  }

  //BiPoly(vp)
template <class NT>
BiPoly<NT>::BiPoly(std::vector<Polynomial<NT> > vp){ 
    // From vector of Polynomials
    ydeg = vp.size() - 1;
    if(ydeg >=0){
      coeffX = vp;
    }
  }
  
  //BiPoly(p, flag):
  //	if true, it converts polynomial p(X) into P(Y)
  // 	if false, it creates the polynomial Y-p(X)

template <class NT>
BiPoly<NT>::BiPoly(Polynomial<NT> p, bool flag){
    if (flag){
      ydeg = p.getTrueDegree();
      if(ydeg >=0){
        for(int i=0; i<=ydeg; i++){
	  Polynomial<NT> temp(0);
	  temp.setCoeff(0, p.getCoeffi(i));
	  coeffX.push_back(temp);	// does STL make a copy of temp?
        }//for
      }//if
    } else {
      ydeg = 1;
      coeffX.push_back(p);
      coeffX.push_back(Polynomial<NT>::polyUnity());
    }//else
  }//BiPoly(p)


  //BiPoly(deg, d[], C[]):
  //	Takes in a list of list of coefficients.
  //	Each cofficient list represents a polynomial in X
  //
  //  deg - ydeg of the bipoly
  //  d[] - array containing the degrees of each coefficient (i.e., X poly)
  //  C[] - list of coefficients, we use array d to select the
  //      coefficients

template <class NT>
BiPoly<NT>::BiPoly(int deg, int *d, NT *C){

    ydeg = deg;
    Polynomial<NT> temp;
    int max=0;
    for(int i=0; i <=deg; i++)
      max = core_max(d[i],max);

    NT *c = new NT[max];

    for(int i=0; i<= deg; i++){
      for(int j=0; j <=d[i]; j++)
	c[j] = C[i+j];
      temp = Polynomial<NT>(d[i],c);
      coeffX.push_back(temp);      
    }
    delete[] c;
  }//BiPoly(deg,d[],C[])


//The BNF syntax is the following:-
//    [bipoly] -> [term]| [term] +/- [bipoly]
//    [term] -> [basic term]|[basic term] [term]|[basic term]*[term]
//    [basic term] -> [number]|'x'|'y'|[basic term]'^'[number]
//                    | '(' [bipoly]')'|'-'
//Unary minus is treated as a basic term
template <class NT>
BiPoly<NT>::BiPoly(const char * s, char myX, char myY){
	string ss(s);
	constructFromString(ss, myX, myY);
}
template <class NT>
BiPoly<NT>::BiPoly(const string & s, char myX, char myY){
	string ss(s);
	constructFromString(ss, myX, myY);
}
template <class NT>
void BiPoly<NT>::constructFromString(string & s, char myX, char myY){
  if((myX != 'x' || myX != 'X') && (myY != 'y' || myY != 'Y')){
    //Replace myX with 'x' and myY with 'y' in s.
    unsigned int loc = s.find(myX, 0);
    while(loc != string::npos){
      s.replace(loc,1,1,'x');
      loc = s.find(myX, loc+1);
    }
    loc = s.find(myY, 0);
    while(loc != string::npos){
      s.replace(loc,1,1,'y');
      loc = s.find(myY, loc+1);
    }
  }
  (*this) =  (*this).getbipoly(s);

}

//Constructor from another BiPoly
template <class NT>
BiPoly<NT>::BiPoly(const BiPoly<NT> &P) : ydeg(P.ydeg), coeffX(P.coeffX) {
}

  //Destructor
template <class NT>
void BiPoly<NT>::deleteCoeffX(){
  coeffX.clear();
}

template <class NT>
BiPoly<NT>::~BiPoly(){
  if (ydeg >= 0)
    deleteCoeffX();
}




  ////////////////////////////////////////////////////////
  // METHODS USED BY STRING CONSTRUCTOR ABOVE
  ////////////////////////////////////////////////////////


//Sets the input BiPoly to X^n
template <class NT>
void BiPoly<NT>::constructX(int n, BiPoly<NT>& P){
  assert(n>= -1);
  P.deleteCoeffX();//Clear the present coeffecients
  Polynomial<NT> q(n);//Nominal degree n
  q.setCoeff(n,NT(1)); 
  if (n>0) q.setCoeff(0,NT(0));
  P.coeffX.push_back(q);
  P.ydeg = 0;
}

//Sets the input BiPoly to Y^n
template <class NT>
void BiPoly<NT>::constructY(int n, BiPoly<NT>& P){
  P.ydeg = 0;
  P.deleteCoeffX();//Clear the present coeffecients
  NT cs[] = {1};
  Polynomial<NT> q(0, cs);
  P.coeffX.push_back(q);//P is the unit bipoly. Now we just multiply Y^n
  P.mulYpower(n);
}

//Returns in P the coeffecient starting from start
template <class NT>
int BiPoly<NT>::getnumber(const char* c, int start, unsigned int len,
			  BiPoly<NT> & P){
  int j=0;
  char *temp = new char[len];
  while(isint(c[j+start])){
    temp[j]=c[j+start];j++;
  }
  temp[j] = '\0';
  NT cf = NT(temp);
  P.deleteCoeffX();
  Polynomial<NT> q(0);
  q.setCoeff(0, cf);
  P.coeffX.push_back(q);
  P.ydeg = 0;
  delete[] temp;

  return (j-1+start);//Exactly the length of the number
}

template <class NT>
bool BiPoly<NT>::isint(char c){
  if(c == '0' || c == '1' || c == '2' || c == '3' || c == '4' ||
     c == '5' || c == '6' || c == '7' || c == '8' || c == '9')
    return true;
  else
    return false;
}


//Returns as integer the number starting from start in c
template <class NT>
int BiPoly<NT>::getint(const char* c, int start, unsigned int len,
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
int BiPoly<NT>::matchparen(const char* cstr, int start){
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
int BiPoly<NT>::getbasicterm(string s, BiPoly<NT> & P){
  const char * cstr = s.c_str();
  unsigned int len = s.length();
  int i=0;
  //BiPoly<NT> * temp = new BiPoly<NT>();

  if(isint(cstr[i])){
    i = getnumber(cstr, i, len, P);
  }else if(cstr[i] == 'x'||cstr[i] == 'X'){
    constructX(1, P);
  }else if(cstr[i] == 'y'||cstr[i] == 'Y'){
    constructY(1, P);
  }else if(cstr[i] =='('){
    int oldi = i;
    i = matchparen(cstr, i);
    string t = s.substr(oldi+1, i -oldi -1);
    P = getbipoly(t);
  }else if(cstr[i] == '-'){//Unary Minus
    P.ydeg = 0;
    P.deleteCoeffX();//Clear the present coeffecients
    NT cs[] = {-1};
    Polynomial<NT> q(0, cs);
    P.coeffX.push_back(q);
  }else{
    std::cout <<"ERROR IN PARSING BASIC TERM" << std::endl;
  }
  //i+1 points to the beginning of next syntaxtic object in the string.
 if(cstr[i+1] == '^'){
    int n;
    i = getint(cstr, i+2, len, n);
    P.pow(n);
  }
  return i;
}


template <class NT>
int BiPoly<NT>::getterm(string s, BiPoly<NT> & P){
  unsigned int len = s.length();
  if(len == 0){// Zero BiPoly
    P = BiPoly<NT>();
    return 0;
  }
  unsigned int ind, oind;
  const char* cstr =s.c_str();
  string t;
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

    BiPoly<NT> R;
    ind = oind + getbasicterm(t, R);//Because the second term is the offset in
                                     //t
    P *= R;
  }

  return ind;
}

template <class NT>
BiPoly<NT> BiPoly<NT>::getbipoly(string s){

    //Remove white spaces from the string
    unsigned int cnt=s.find(' ',0);
    while(cnt != string::npos){
      s.erase(cnt, 1);
      cnt = s.find(' ', cnt);
    }

    unsigned int len = s.length();
    if(len <= 0){//Zero Polynomial
      return BiPoly<NT>();
    }

    //To handle the case when there is one '=' sign
    //Suppose s is of the form s1 = s2. Then we assign s to
    //s1 + (-1)(s2) and reset len
    unsigned int loc;
    if((loc=s.find('=',0)) != string::npos){
      s.replace(loc,1,1,'+');
      string s3 = "(-1)(";
      s.insert(loc+1, s3);
      len = s.length();
      s.insert(len, 1, ')');
    }
    len = s.length();

    const char *cstr = s.c_str();
    string t;
    BiPoly<NT> P;
    // P will be the polynomial in which we accumulate the
    //sum and difference of the different terms.
    unsigned int ind;

    if(cstr[0] == '-'){
      t = s.substr(1, len);
      ind = getterm(t,P) + 1;
      NT negone(-1);
      P.mulScalar(negone);
    }else{
      ind = getterm(s, P);
    }

    unsigned int oind =0;//the string between oind and ind is a term
    while(ind != len -1){
      BiPoly<NT> R;
      t = s.substr(ind + 2, len -ind -2);
      oind = ind;
      ind = oind + 2 + getterm(t, R);
      if(cstr[oind + 1] == '+')
	P += R;
      else if(cstr[oind + 1] == '-')
	P -= R;
      else
	std::cout << "ERROR IN PARSING BIPOLY! " << std::endl;
    }

    return P;
}

  ////////////////////////////////////////////////////////
  // METHODS
  ////////////////////////////////////////////////////////
  
  // filedump (msg, ofs, com, com2)
  // 	where msg and com are strings.
  // 	msg is an message and com is the character preceding each line
  // 	(e.g., msg=""  and com=com2="# ")
  // This is called by the other dump functions
template <class NT>
void BiPoly<NT>::dump(std::ostream & os, std::string msg,
    std::string com, std::string com2) const {
    if (msg != "")
      os << msg << std::endl;
    if(ydeg == -1) {
      os << com << " Zero Polynomial" << std::endl;
      return;
    }
    bool first = true;
    os << com;
    for (int i=0; i <= ydeg; i++){
      if (!zeroP(coeffX[i])){
	if (i % 3 == 0) os << std::endl << com2 ;  // output 3 coefficients per line
	if (first) first = false;
	else os << " + ";
	os << "[";
	coeffX[i].filedump(os,"","",com2); // Note: first comment string is ""
	os << "]";
        if (i > 0) os <<" * y^"<< i ;
      }
    }
  }//dump

  // dump(ofs, msg, com, com2) -- dump to file
  //
/*
template <class NT>
void BiPoly<NT>::dump(std::ofstream & ofs, std::string msg,
              std::string com, std::string com2) const {
    dump(ofs, msg, com, com2);
  }
*/
  // dump(msg, com) -- dump to std output
template <class NT>
void BiPoly<NT>::dump(std::string msg, std::string com,
    std::string com2) const {
    dump(std::cout, msg, com, com2);
  }

  /* ***********************************************************
    We want to substitute X or Y with a general Expression or BigFloat
    into a BiPoly.  E.g., the following produces a Y-polynomial
    after replacing X by x:

  BiPoly<NT> yPolynomial(const Expr & x) {

    VecExpr vE;
    for (int i=0; i<= ydeg; i++) {
        vE.push_back(coeffX[i].eval(x));
    }
    return Polynomial<Expr>(vE);
  }//yPolynomial
    
    But this has many problems.
    Solution below is to have special yExprPolynomial(x).
    *********************************************************** */
  
template <class NT>
Polynomial<NT> BiPoly<NT>::yPolynomial(const NT & x) {
    NT coeffVec[ydeg+1];
    int d=-1;
    for(int i=ydeg; i >= 0 ; i--){
      coeffVec[i] = coeffX[i].eval(x);
      if ((d < 0) && (coeffVec[i] != 0))
	      d = i;
    }
    return Polynomial<NT>(d, coeffVec);
  }

  // Specialized version of yPolynomial for Expressions
template <class NT>
  Polynomial<Expr> BiPoly<NT>::yExprPolynomial(const Expr & x) {
    Expr coeffVec[ydeg+1];
    int d=-1;
    for(int i=ydeg; i >= 0 ; i--){
      coeffVec[i] = coeffX[i].eval(x);
      if ((d < 0) && (coeffVec[i] != 0))
	      d = i;
    }
    return Polynomial<Expr>(d, coeffVec);
  }

  // Specialized version of yPolynomial for BigFloat
  // ASSERTION: BigFloat x is exact
template <class NT>
  Polynomial<BigFloat> BiPoly<NT>::yBFPolynomial(const BigFloat & x) {
    BigFloat coeffVec[ydeg+1];
    int d=-1;
    for(int i=ydeg; i >= 0 ; i--){
      coeffVec[i] = coeffX[i].eval(x);
      if ((d < 0) && (coeffVec[i] != 0))
	      d = i;
    }
    return Polynomial<BigFloat>(d, coeffVec);
  }

  // xPolynomial(y) 
  //   returns the polynomial (in X) when we substitute Y=y
  //   
  //   N.B. May need the
  //   		Polynomial<Expr> xExprPolynomial(Expr y)
  //   version too...
  //

template <class NT>
Polynomial<NT> BiPoly<NT>::xPolynomial(const NT & y) {
    Polynomial<NT> res = coeffX[0];
    NT yPower(y);
    for(int i=1; i <= ydeg ; i++){
      res += coeffX[i].mulScalar(yPower);
      yPower *= y;
    }
    return res;
  }//xPolynomial
  

  // getYdegree()
template <class NT>
int BiPoly<NT>::getYdegree() const {
    return ydeg;
  }
  
  // getXdegree()
template <class NT>
int BiPoly<NT>::getXdegree(){
    int deg=-1;
    for(int i=0; i <=ydeg; i++)
      deg = max(deg, coeffX[i].getTrueDegree());
    return deg;
  }

  // getTrueYdegree
template <class NT>
int BiPoly<NT>::getTrueYdegree(){
    for (int i=ydeg; i>=0; i--){
      coeffX[i].contract();
      if (!zeroP(coeffX[i]))
	return i;
    }
    return -1;	// Zero polynomial
  }//getTrueYdegree


  //eval(x,y)
template <class NT>
Expr BiPoly<NT>::eval(Expr x, Expr y){//Evaluate the polynomial at (x,y)
    Expr e = 0;

    for(int i=0; i <=ydeg; i++){
      e += coeffX[i].eval(x)*pow(y,i);
    }
    return e;
  }//eval

  ////////////////////////////////////////////////////////
  // Polynomial arithmetic (these are all self-modifying)
  ////////////////////////////////////////////////////////
  
  // Expands the nominal y-degree to n;
  //	Returns n if nominal y-degree is changed to n
  //	Else returns -2

template <class NT>
int BiPoly<NT>::expand(int n) {
    if ((n <= ydeg)||(n < 0))
      return -2;
    
    for(int i=ydeg+1; i <=n ;i++)
      coeffX.push_back(Polynomial<NT>::polyZero());
    
    ydeg = n;
    return n;
  }//expand

  // contract() gets rid of leading zero polynomials
  //	and returns the new (true) y-degree;
  //	It returns -2 if this is a no-op

template <class NT>
int BiPoly<NT>::contract() {
    int d = getTrueYdegree();
    if (d == ydeg)
      return (-2);  // nothing to do
    else{
      for (int i = ydeg; i> d; i--)
	coeffX.pop_back();
      ydeg = d;
    }
    return d;
  }//contract

  // Self-assignment
template <class NT>
BiPoly<NT> & BiPoly<NT>::operator=( const BiPoly<NT>& P) {
  if (this == &P)
    return *this;	// self-assignment
  ydeg = P.getYdegree();
  coeffX = P.coeffX;
  return *this;
  }//operator=


  // Self-addition
template <class NT>
BiPoly<NT> & BiPoly<NT>::operator+=( BiPoly<NT>& P) { // +=

    int d = P.getYdegree();
    if (d > ydeg)
      expand(d);
    for (int i = 0; i<=d; i++)
      coeffX[i] += P.coeffX[i];

  return *this;
  }//operator+=
   
  // Self-subtraction
template <class NT>
BiPoly<NT> & BiPoly<NT>::operator-=( BiPoly<NT>& P) { // -=
    int d = P.getYdegree();
    if (d > ydeg)
      expand(d);
  for (int i = 0; i<=d; i++)
    coeffX[i] -= P.coeffX[i];

  return *this;
  }//operator-=

  // Self-multiplication
template <class NT>
BiPoly<NT> & BiPoly<NT>::operator*=( BiPoly<NT>& P) { // *=
    int d = ydeg + P.getYdegree();

    std::vector<Polynomial<NT> > vP;

    Polynomial<NT>* c = new Polynomial<NT> [d + 1];
    for(int i=0; i <=d; i++)
      c[i] = Polynomial<NT>();
    
    for (int i = 0; i<=P.getYdegree(); i++)
      for (int j = 0; j<=ydeg; j++) {
	if(!zeroP(P.coeffX[i]) && !zeroP(coeffX[j]))
	  c[i+j] += P.coeffX[i] * coeffX[j];
      }

    for(int i=0; i <= d; i++)
      vP.push_back(c[i]);

    delete[] c;
    coeffX.clear();
    coeffX = vP;
    ydeg = d;
    return *this;
  }//operator*=
  
  // Multiply by a polynomial in X
template <class NT>
BiPoly<NT> & BiPoly<NT>::mulXpoly( Polynomial<NT> & p) {
    contract();
    if (ydeg == -1)
      return *this;

    for (int i = 0; i<=ydeg ; i++)
      coeffX[i] *= p;
    return *this;
  }//mulXpoly

  //Multiply by a constant
template <class NT>
BiPoly<NT> & BiPoly<NT>::mulScalar( NT & c) {
    for (int i = 0; i<=ydeg ; i++)
      coeffX[i].mulScalar(c);
    return *this;
  }//mulScalar

  // mulYpower: Multiply by Y^i (COULD be a divide if i<0)
template <class NT>
BiPoly<NT> & BiPoly<NT>::mulYpower(int s) {
  // if s >= 0, then this is equivalent to
  // multiplying by Y^s;  if s < 0, to dividing by Y^s
  if (s==0)
    return *this;
  int d = s+getTrueYdegree();
  if (d < 0) {
    ydeg = -1;
    deleteCoeffX();
    return *this;
  }

  std::vector<Polynomial<NT> > vP;
  if(s > 0){
    for(int i=0; i < s; i ++)
      vP.push_back(Polynomial<NT>::polyZero());
    for(int i=s; i<=d; i++)
      vP.push_back(coeffX[i-s]);
  }

  if (s < 0) {
    for (int j=-1*s; j <= d-s; j++)
      vP.push_back(coeffX[j]);
  }

  coeffX = vP;
  ydeg= d;
  return *this;
  }//mulYpower


  
  // Divide by a polynomial in X.
  // We replace the coeffX[i] by the pseudoQuotient(coeffX[i], P)
template <class NT>
BiPoly<NT> & BiPoly<NT>::divXpoly( Polynomial<NT> & p) {
    contract();
    if (ydeg == -1)
      return *this;

    for (int i = 0; i<=ydeg ; i++)
      coeffX[i] = coeffX[i].pseudoRemainder(p);

    return *this;
  }// divXpoly


  
  //Using the standard definition of pseudRemainder operation.
  //	--No optimization!
template <class NT>
BiPoly<NT> BiPoly<NT>::pseudoRemainderY (BiPoly<NT> & Q){
    contract();
    Q.contract();
    int qdeg = Q.getYdegree();
    Polynomial<NT> LC(Q.coeffX[qdeg]), temp1(LC), temp;
    BiPoly<NT> P(Q);

    std::vector<Polynomial<NT> > S;//The quotient in reverse order
    if(ydeg < qdeg){
      return (*this);
    }
    (*this).mulXpoly(temp1.power(ydeg - qdeg +1));

    while(ydeg >= qdeg){
      temp1 = coeffX[ydeg];
      temp = temp1.pseudoRemainder(LC);
      S.push_back(temp);
      P.mulXpoly(temp);
      P.mulYpower(ydeg - qdeg);
      *this -= P;
      contract();
      P = Q;//P is used as a temporary since mulXpower and mulYpower are
            // self modifying
    }

    //Correct the order of S
    std::vector<Polynomial<NT> > vP;

    for(int i= S.size()-1; i>=0; i--)
      vP.push_back(S[i]);


    return BiPoly<NT>(vP);
    
  }//pseudoRemainder

  //Partial Differentiation
  //Partial Differentiation wrt Y
template <class NT>
BiPoly<NT> & BiPoly<NT>::differentiateY() {	
  if (ydeg >= 0) {
    for (int i=1; i<=ydeg; i++)
      coeffX[i-1] = coeffX[i].mulScalar(i);
    ydeg--;
  }
  return *this;
  }// differentiation wrt Y

template <class NT>
BiPoly<NT> & BiPoly<NT>::differentiateX() {	
    if (ydeg >= 0)
      for (int i=0; i<=ydeg; i++)
	coeffX[i].differentiate();
    
    return *this;
  }// differentiation wrt X

template <class NT>
BiPoly<NT> & BiPoly<NT>::differentiateXY(int m, int n) {//m times wrt X and n times wrt Y
    assert(m >=0); assert(n >=0);
    for(int i=1; i <=m; i++)
      (*this).differentiateX();
    for(int i=1; i <=n; i++)
      (*this).differentiateY();
    
    return *this;
  }

  //Represents the bivariate polynomial in (R[X])[Y] as a member
  //of (R[Y])[X].
  //But since our polynomials in X can only have NT coeff's thus
  // to represent the above polynomial we switch X and Y once
  // the conversion has been done.
  //NOTE: This is different from replacing X by Y which was
  //      done in the case of the constructor from a polynomial in X
  //Need to calculate resultant wrt X.
template <class NT>
BiPoly<NT> & BiPoly<NT>::convertXpoly(){
    getTrueYdegree();
    if(ydeg == -1) return (*this);
    NT *cs = new NT[ydeg +1];
    int xdeg = getXdegree();
    std::vector<Polynomial<NT> > vP;

    for(int i=0; i<=xdeg; i++){
      for(int j=0; j<=ydeg; j++){
	cs[j] = coeffX[j].getCoeffi(i);
      }
      
      vP.push_back(Polynomial<NT>(ydeg, cs));
    }
    delete[] cs;
      
    ydeg = xdeg;
    coeffX = vP;
    return (*this);
  }

  //Set Coeffecient to the polynomial passed as a parameter
template <class NT>
bool BiPoly<NT>::setCoeff(int i, Polynomial<NT> p){
    if(i < 0 || i > ydeg)
      return false;
    coeffX[i] = p;
    return true;
  }
  
template <class NT>
void BiPoly<NT>::reverse() {
    Polynomial<NT> tmp;
    for (int i=0; i<= ydeg/2; i++) {
      tmp = coeffX[i];
      coeffX[i] =   coeffX[ydeg-i];
      coeffX[ydeg-i] = tmp;
    }
  }//reverse

  //replaceYwithX()
  //   used used when the coeffX in BiPoly are constants,
  //   to get the corresponding univariate poly in X
  //   E.g., Y^2+2*Y+9 will be converted to X^2+2*X+9
template <class NT>
Polynomial<NT>  BiPoly<NT>::replaceYwithX(){
    int m = getTrueYdegree();
    NT *cs = new NT[m+1];
    for(int i=0; i <= m ; i++)
      cs[i]=coeffX[i].getCoeffi(0);
    delete[] cs;

    return Polynomial<NT>(m,cs);
  }//replaceYwithX
  
template <class NT>
BiPoly<NT>&  BiPoly<NT>::pow(unsigned int n){

  if (n == 0) {
    ydeg = 0;
    deleteCoeffX();
    Polynomial<NT> unity(0);
    coeffX.push_back(unity);
  }else if(n == 1){

  }else{
    BiPoly<NT> x(*this);

    while ((n % 2) == 0) { // n is even
      x *= x;
      n >>= 1;
    }
    BiPoly<NT> u(x);
    while (true) {
      n >>= 1;
      if (n == 0){
	(*this) = u;
	return (*this);
      }
      x *= x;
      if ((n % 2) == 1) // n is odd
        u *= x;
    }
    (*this) = u;
  }
  return (*this);
}//pow



  ////////////////////////////////////////////////////////
  // Helper Functions
  ////////////////////////////////////////////////////////

// isZeroPinY(P)
//	checks whether a Bi-polynomial is a zero Polynomial
template <class NT>
bool isZeroPinY(BiPoly<NT> & P){
    if(P.getTrueYdegree() == -1)
      return true;
    return false;
}

// gcd(P,Q)
//   This gcd is based upon the subresultant PRS to avoid
//   exponential coeffecient growth and gcd computations, both of which 
//   are expensive since the coefficients are polynomials

template <class NT>
BiPoly<NT> gcd( BiPoly<NT>& P ,BiPoly<NT>& Q){
    int m = P.getTrueYdegree();
    int n = Q.getTrueYdegree();

    //If both Bi-Polys are zero then return zero
    if( m==-1 && n==-1){
      return BiPoly<NT>();//Had to do this for the type of
                          //return value o/w causes problem in
    }                    //assignment operation

    if(m < n) {
      return gcd(Q, P);
    }

    //If one of the Bi-Polys are zero then return the other
    if(n==-1)
      return P;

    //When the two BiPolys are just univariate Polynomials in X
    if(m==0 && n==0){
      std::vector<Polynomial<NT> > vP;
      vP.push_back(gcd(P.coeffX[0], Q.coeffX[0]));
      return BiPoly<NT>(vP);//Had to do this for the type of
                            //return value
      }
    
    int delta = m - n;
    Polynomial<NT> a(Q.coeffX[n]);
    std::vector<NT> vN;
    vN.push_back(pow(BigInt(-1),delta + 1));
    Polynomial<NT> beta(vN);
    Polynomial<NT> phi(power(a, delta)), t;
    m = n;
    BiPoly<NT> temp;
    P.pseudoRemainderY(Q);

    while(!isZeroPinY(P)){
      P.divXpoly(beta);
      n = P.getTrueYdegree();
      delta = m - n;
      beta = power(phi, delta)*a.mulScalar(pow(BigInt(-1),delta+1));
      a = P.coeffX[n];
      t = power(phi, delta -1);
      phi = power(a,delta).pseudoRemainder(t);
      m = n;
      temp = Q;//Swapping the pseudo-remainder for the next step
      Q = P;
      P = temp;
      P.pseudoRemainderY(Q);
  }
    return Q;
}//gcd

// resY(P,Q):
//      Resultant of Bi-Polys P and Q w.r.t. Y.
//      So the resultant is a polynomial in X
template <class NT>
Polynomial<NT>  resY( BiPoly<NT>& P ,BiPoly<NT>& Q){

  int m = P.getTrueYdegree();
  int n = Q.getTrueYdegree();

  //If either polynomial is zero, return zero
  if( m == -1 || n == -1) return Polynomial<NT>();//::polyZero();

  if(n > m)
    return resY(Q, P);

  Polynomial<NT> b(Q.coeffX[n]);
  Polynomial<NT> lc(P.coeffX[m]), C, temp;
  BiPoly<NT> r;
  
  r = P.pseudoRemainderY(Q);
  C = b * r.coeffX[r.getTrueYdegree()];
  C = C.pseudoRemainder(lc);

  if(isZeroPinY(P) && n > 0)
    return Polynomial<NT>();//::polyZero();

  if(Q.getTrueYdegree() == 0 && isZeroPinY(P))
    return power(Q.coeffX[0], m);

  int l = P.getTrueYdegree();

  temp = power(b, m-l).mulScalar(pow(BigInt(-1),m*n))*resY(Q,P);
  temp = temp.pseudoRemainder(power(C,n));
  return temp;

}//resY

// resX(P,Q):
//      Resultant of Bi-Polys P and Q w.r.t. X.
//      So the resultant is a polynomial in Y
//	We first convert P, Q to polynomials in X. Then 
// 	call resY and then turns it back into a polynomial in Y
//	QUESTION: is this last switch really necessary???
template <class NT>
BiPoly<NT>  resX( BiPoly<NT>& P ,BiPoly<NT>& Q){
  P.convertXpoly();
  Q.convertXpoly();

  // Polynomial<NT> p(resY(P,Q));
  // return BiPoly<NT>(p);
  return(resY(P,Q)); // avoid the switch back
}//resX


//Equality operator for BiPoly
template <class NT>
bool operator==(const BiPoly<NT>& P, const BiPoly<NT>& Q) {	// ==
  BiPoly<NT> P1(P);
  BiPoly<NT> Q1(Q);
  P1.contract();
  Q1.contract();
  if(P1.getYdegree() != Q1.getYdegree())
    return false;
  else{
    for(int i=0; i <= P1.getYdegree() ; i++){
      if(P1.coeffX[i] != Q1.coeffX[i])
	return false;
    }
  }
  return true;
}

  // Addition P + Q
template <class NT>
BiPoly<NT> operator+(const BiPoly<NT>& P, const BiPoly<NT>& Q ) { // +
  return BiPoly<NT>(P) += Q;
}
  // Subtraction P - Q
template <class NT>
BiPoly<NT> operator-(const  BiPoly<NT>& P,const  BiPoly<NT>& Q ) { // +
  return BiPoly<NT>(P) -= Q;
}

  // Multiplication P * Q
template <class NT>
BiPoly<NT> operator*(const  BiPoly<NT>& P, const BiPoly<NT>& Q ) { // +
  return BiPoly<NT>(P) *= Q;
}


  ////////////////////////////////////////////////////////
  //Curve Class
  //  	extends BiPoly Class
  ////////////////////////////////////////////////////////

template < class NT >
Curve<NT>::Curve(){} // zero polynomial
  
  //Curve(vp):
  //    construct from a vector of polynomials
template < class NT >
Curve<NT>::Curve(std::vector<Polynomial<NT> > vp)
	  : BiPoly<NT>(vp){
  }
  
  //Curve(p):
  //	Converts a polynomial p(X) to a BiPoly in one of two ways:
  // 	    (1) if flag is false, the result is Y-p(X) 
  // 	    (2) if flag is true, the result is p(Y) 
  //    The default is (1) because we usually want to plot the
  //        graph of the polynomial p(X)
template < class NT >
Curve<NT>::Curve(Polynomial<NT> p, bool flag)
	  : BiPoly<NT>(p, flag){
  }

  //Curve(deg, d[], C[]):
  //	Takes in a list of list of coefficients.
  //	Each cofficient list represents a polynomial in X
  //
  //  deg - ydeg of the bipoly
  //  d[] - array containing the degrees of each coefficient (i.e., X poly)
  //  C[] - list of coefficients, we use array d to select the
  //      coefficients
template < class NT >
Curve<NT>::Curve(int deg, int *d, NT *C)
	  : BiPoly<NT>(deg, d, C){
  }

template < class NT >
Curve<NT>::Curve(const BiPoly<NT> &P)
	  : BiPoly<NT>(P){
  }

  //Curve(n)
template < class NT >
Curve<NT>::Curve(int n)
	  : BiPoly<NT>(n){// creates a Curve with nominal y-degree equal to n
  }

  //Creates a curve from a string (no parentheses, no *)
template < class NT >
Curve<NT>::Curve(const string & s, char myX, char myY)
	  : BiPoly<NT>(s, myX, myY){
  }
template < class NT >
Curve<NT>::Curve(const char * s, char myX, char myY)
	  : BiPoly<NT>(s, myX, myY){
  }


  // verticalIntersections(x, vecI, aprec=0):
  //
  //    The list vecI is passed an isolating intervals for y's such that (x,y)
  //    lies on the curve.
  //    If aprec is non-zero (!), the intervals have with < 2^{-aprec}.
  //    Return is -2 if curve equation does not depend on Y
  //    	-1 if infinitely roots at x,
  //    	0 if no roots at x
  //    	1 otherwise
  //
  //    ASSERTION: x is an exact BigFloat

template < class NT >
int Curve<NT>::verticalIntersections(const BigFloat & x, BFVecInterval & vI,
	int aprec) {
    int d= Curve<NT>::getTrueYdegree();
    if(d <= 0) return(-2);
    	   // This returns a NULL vI, which should be caught by caller
	
    Polynomial<Expr> PY = this->yExprPolynomial(x); // should be replaced
    // assert(x.isExact());
    // Polynomial<BigFloat> PY = yBFPolynomial(x); // unstable still

    d = PY.getTrueDegree();
    if(d <= 0) return(d);

    Sturm<Expr> Ss(PY); // should be replaced by BigFloat version
    // Sturm<BigFloat> Ss(PY); // unstable still
    Ss.isolateRoots(vI);

    int s = vI.size();
    if ((aprec != 0) && (s>0))
	Ss.newtonRefineAllRoots(vI, aprec);
    
    return s;
  }
  
  // plot(eps, x1, y1, x2, y2)
  //
  // 	All parameters have defaults
  //
  //    Gives the points on the curve at resolution "eps".  Currently,
  //    eps is viewed as delta-x step size (but it could change).
  //    The display is done in the rectangale 
  //    defined by [(x1, y1), (x2, y2)].
  //    The output is written into a file in the format specified
  //    by our drawcurve function (see COREPATH/ext/graphics).
  //
  //    Heuristic: the open polygonal lines end when number of roots
  //    changes...
  //
  //    TO DO:
  //       (1) need to automatically readjust x- and y-scales
  //              independently
  //       (2) reorder parameters in this order: x1, x2, y1, y2
  //            [Then y1, y2 can be automatically determined]
  //       (3) should allow the ability to look for interesting
  //             features
  //
  //    ASSERTION: all input BigFloats are exact
  //
template < class NT >
int Curve<NT>::plot( BigFloat eps, BigFloat x1,
	BigFloat y1, BigFloat x2, BigFloat y2, int fileNo){

  const char* filename[] = {"data/input", "data/plot", "data/plot2"};

  assert(eps.isExact()); // important for plotting...
  assert(x1.isExact());
  assert(y1.isExact());

  ofstream outFile;
  outFile.open(filename[fileNo]); // ought to check if open is successful!
  outFile << "########################################\n";
  outFile << "# Curve equation: \n";
  this->dump(outFile,"", "# ");
  outFile << std::endl;
  outFile << "# Plot parameters:  step size (eps) = " << eps << std::endl;
  outFile << "#                   x1 = " << x1 << ",\t y1 = " << y1 << std::endl;
  outFile << "#                   x2 = " << x2 << ",\t y2 = " << y2 << std::endl;
  outFile << "########################################\n";
  outFile << "# X-axis " << std::endl;
  outFile << "o 2" << std::endl;
  outFile << "0.99 \t 0.99 \t 0.99" << std::endl;
  outFile << x1 << "\t" << 0 << std::endl;
  outFile << x2 << "\t" << 0 << std::endl;
  outFile << "########################################\n";
  outFile << "# Y-axis " << std::endl;
  outFile << "o 2" << std::endl;
  outFile << "0.99 \t 0.99 \t 0.99" << std::endl;
  outFile << 0 << "\t" << y1 << std::endl;
  outFile << 0 << "\t" << y2 << std::endl;
  // assert(eps>0)
  int aprec = -(eps.lMSB()).toLong(); // By definition, eps.lMSB() <= lg(eps)

  BigFloat xCurr = x1;
  BigFloat xLast = x1;
  unsigned int numRoots=0;
  unsigned int numPoints=0;
  BFVecInterval vI;

cout <<"Current value of x " << xCurr << endl;
  //===================================================================
  // FIRST locate the first x-value where there are roots for plotting
  //===================================================================
  do {
    vI.clear();
    if (verticalIntersections(xCurr, vI, aprec) > 0) {
      numRoots = vI.size();
cout <<"Number of roots at " << xCurr << " are " << numRoots<<endl;
    }
    xCurr += eps;

  } while ((numRoots == 0) && (xCurr < x2));//numRoots <= ydeg

  if (numRoots == 0 && x2 <= xCurr) {//if numRoots > 0 then there exists a
                                     //valid point for plotting
	  return -1; // nothing to plot!
  }

  int limit = ((x2 - xCurr + eps)/eps).intValue()+1;
  //std::cout << "Limit = " << limit << std::endl;
  machine_double plotCurves[this->getTrueYdegree()][limit];//plot buffer 
  machine_double yval;

  for (unsigned int i=0; i< numRoots; i++) {
     yval = (vI[i].first + vI[i].second).doubleValue()/2;
     if (yval < y1) plotCurves[i][numPoints] = y1.doubleValue();
     else if (yval > y2) plotCurves[i][numPoints] = y2.doubleValue();
     else plotCurves[i][numPoints] = yval;
  }

  vI.clear();
  xLast = xCurr - eps;  // -eps to compensate for the forward value of xCurr
  numPoints = 1;

  //===================================================================
  // Get all the curves in a main loop
  // 	-- dump the curves when an x-interval is discovered
  // 	We define an "x-interval" to be a maximal interval
  // 	where the number of curves is constant (this is a heuristic!)
  // 	Note that this includes the special case where number is 0.
  //===================================================================

  BigFloat tmp; // used to step from xLast to xCurr in loops

  while (xCurr < x2) { //main loop
    //std::cout << "Doing verticalintersec at " << xCurr << std::endl;
    verticalIntersections(xCurr, vI, aprec);
    if (vI.size() != numRoots) { // an x-interval discovered!
	// write previous x-interval to output file
        outFile << "########################################\n";
        outFile << "# New x-interval with " << numRoots << " roots\n";
	for (unsigned int i=0; i< numRoots; i++) {
          outFile << "#=======================================\n";
          outFile << "# Curve No. " << i+1 << std::endl;
	  outFile << "o " << numPoints << std::endl;
	  outFile << red_comp(i) << "\t"
		  << green_comp(i)  << "\t"
		  << blue_comp(i) << std::endl;
	  tmp = xLast;
	  for (unsigned int j=0; j< numPoints; j++) {

		  outFile << tmp << "	\t\t"
			  << plotCurves[i][j] << "\n";
		  tmp += eps;
	  }//for j
	}//for i
	numPoints = 0;          // reset
	numRoots = vI.size();   // reset
	xLast = xCurr;		// reset
    }//if vI.size() !=
    if (numRoots>0){ // record curr. vertical intersections if numRoots>0
      for (unsigned int i=0; i< numRoots; i++) {
         yval = (vI[i].first + vI[i].second).doubleValue()/2;
	 // HERE SHOULD BE A LOOP TO OUTPUT MORE POINTS IN CASE THE slope IS LARGE
	 // Idea: let previous value of yval be yval-old.  
	 //
	 // Two cases:
	 // (1) When i=0:
	 // 	you need to search backwards and forwards
	 // 	(yval-old is not defined in this case)
	 //
	 // (2) When i>0:
	 // 	If  |yval-old - yval| > 2eps, we must plot using horizontalIntersection()
	 // 	We must start from
	 // 		y = yval-old until
	 // 			EITHER (y = (vI[i+1].first + vI[i+1].second)/2)
	 // 			OR (sign of slope changes)
	 // 			OR (hit the ymax or ymin boundary)
	 //
         if (yval < y1) plotCurves[i][numPoints] = y1.doubleValue();
         else if (yval > y2) plotCurves[i][numPoints] = y2.doubleValue();
         else plotCurves[i][numPoints] = yval;
      }//for i
      vI.clear();  //Should clear the intersection points
      numPoints++;
    }//if
    xCurr += eps;
if (!xCurr.isExact()) std::cout<<"xCurr has error! xCurr=" << xCurr << std::endl;
   }//main while loop

   // Need to flush out the final x-interval:
   if ((numRoots>0) && (numPoints >0)) { 
	// write to output file
        outFile << "########################################\n";
        outFile << "# New x-interval with " << numRoots << " roots\n";
	for (unsigned int i=0; i< numRoots; i++) {
          outFile << "#=======================================\n";
          outFile << "# Curve No. " << i+1 << std::endl;
	  outFile << "o " << numPoints << std::endl;
	  outFile << red_comp(i) << "\t"
		  << green_comp(i)  << "\t"
		  << blue_comp(i) << std::endl;
	  tmp = xLast;
	  for (unsigned int j=0; j< numPoints; j++) {

		  outFile << tmp << "	\t\t"
			  << plotCurves[i][j] << "\n";
		  tmp += eps;
	  }//for j
	}//for i
    }//if

    // Put out the final frame (this hides the artificial cut-off of curves
    outFile << "########################################\n";
    outFile << "# FRAME around our figure " << std::endl;
    outFile << "p 4" << std::endl;
    outFile << "0.99 \t 0.99 \t 0.99" << std::endl;
    outFile << x1 << "\t" << y1 << std::endl;
    outFile << x2 << "\t" << y1 << std::endl;
    outFile << x2 << "\t" << y2 << std::endl;
    outFile << x1 << "\t" << y2 << std::endl;
    outFile << "######### End of File ##################\n";

    outFile.close();
    return 0;
  }//plot

// selfIntersections():
//   this should be another member function that lists
//   all the self-intersections of a curve
//
//  template <class NT>
//  void selfIntersections(BFVecInterval &vI){
//  ...
//  }



  ////////////////////////////////////////////////////////
  // Curve helper functions
  ////////////////////////////////////////////////////////


//Xintersections(C, D, vI):
//  returns the list vI of x-ccordinates of possible intersection points.
//  Assumes that C & D are quasi-monic.(or generally aligned)
template <class NT>
void  Xintersections( Curve<NT>& P ,Curve<NT>& Q, BFVecInterval &vI){
  Sturm<NT> Ss(resY(P, Q));
  Ss.isolateRoots(vI);
}

//Yintersections(C, D, vI):
//	similar to Xintersections
template <class NT>
void  Yintersections( Curve<NT>& P ,Curve<NT>& Q, BFVecInterval &vI){
  Sturm<NT> Ss(resX(P, Q));
  Ss.isolateRoots(vI);
}

// Display Intervals
// 
template <class NT>//DO I NEED THIS OVERHERE AS WELL?
void showIntervals(char* s, BFVecInterval &vI) {
   std::cout << s;
   for (unsigned int i=0; i< vI.size(); i++) {
   	std::cout << "[ " << vI[i].first << ", " 
   		<< vI[i].second << " ],  " ;
   }
   std::cout << std::endl;
}


/*************************************************************************** */
// END
/*************************************************************************** */

