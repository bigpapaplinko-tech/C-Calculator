
#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <numeric>
#include <complex>
#include <map>
#include <functional>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <cassert>

using namespace std;
using cd = complex<double>;

// ─────────────────────────────────────────────────────────────
//  GLOBAL CONSTANTS (Put them all here!)
// ─────────────────────────────────────────────────────────────
const double PI       = acos(-1.0);
const double E_NUM    = exp(1.0);
const double PHI      = (1.0 + sqrt(5.0)) / 2.0;
const double INF      = numeric_limits<double>::infinity();

// Physics Constants for Choice 21 & 22
const double C_LIGHT  = 299792458.0;      // Speed of light
const double H_PLANCK = 6.62607015e-34;   // Planck's constant
const double H_BAR    = H_PLANCK / (2.0 * PI);

// ─────────────────────────────────────────────────────────────
//  EXPRESSION PARSER  (handles full math expressions as strings)
// ─────────────────────────────────────────────────────────────
struct Parser {
    string src;
    size_t pos;

    Parser(const string& s) : src(s), pos(0) {}

    void skipSpaces() { while (pos < src.size() && src[pos] == ' ') pos++; }

    double parseAtom() {
        skipSpaces();
        if (pos < src.size() && src[pos] == '(') {
            pos++; // consume '('
            double v = parseExpr();
            skipSpaces();
            if (pos < src.size() && src[pos] == ')') pos++;
            return v;
        }
        if (pos < src.size() && isalpha(src[pos])) return parseFnOrConst();
        // number
        bool neg = false;
        if (pos < src.size() && src[pos] == '-') { neg = true; pos++; skipSpaces(); }
        if (pos < src.size() && src[pos] == '+') { pos++; skipSpaces(); }
        if (pos < src.size() && src[pos] == '(') {
            pos++;
            double v = parseExpr();
            skipSpaces();
            if (pos < src.size() && src[pos] == ')') pos++;
            return neg ? -v : v;
        }
        double intPart = 0, fracPart = 0, fracDiv = 1;
        while (pos < src.size() && isdigit(src[pos])) intPart = intPart * 10 + (src[pos++] - '0');
        if (pos < src.size() && src[pos] == '.') {
            pos++;
            while (pos < src.size() && isdigit(src[pos])) { fracPart = fracPart * 10 + (src[pos++] - '0'); fracDiv *= 10; }
        }
        double val = intPart + fracPart / fracDiv;
        // scientific notation
        if (pos < src.size() && (src[pos] == 'e' || src[pos] == 'E') && (pos+1 < src.size()) && (isdigit(src[pos+1]) || src[pos+1]=='-' || src[pos+1]=='+')) {
            pos++;
            double expSign = 1;
            if (src[pos] == '-') { expSign = -1; pos++; }
            else if (src[pos] == '+') pos++;
            double expVal = 0;
            while (pos < src.size() && isdigit(src[pos])) expVal = expVal * 10 + (src[pos++] - '0');
            val *= pow(10, expSign * expVal);
        }
        return neg ? -val : val;
    }

    double parseFnOrConst() {
        // read identifier
        string name;
        while (pos < src.size() && (isalpha(src[pos]) || isdigit(src[pos]) || src[pos] == '_')) name += src[pos++];
        skipSpaces();

        // constants
        if (name == "pi" || name == "PI") return PI;
        if (name == "e"  || name == "E")  return E_NUM;
        if (name == "phi" || name == "PHI") return PHI;
        if (name == "inf" || name == "INF") return INF;

        // functions — expect '('
        if (pos < src.size() && src[pos] == '(') {
            pos++; // consume '('
            double a = parseExpr();
            // check for second argument
            double b = 0; bool hasB = false;
            skipSpaces();
            if (pos < src.size() && src[pos] == ',') { pos++; b = parseExpr(); hasB = true; }
            skipSpaces();
            if (pos < src.size() && src[pos] == ')') pos++;

            if (name=="sin")    return sin(a);
            if (name=="cos")    return cos(a);
            if (name=="tan")    return tan(a);
            if (name=="asin")   return asin(a);
            if (name=="acos")   return acos(a);
            if (name=="atan")   return hasB ? atan2(a,b) : atan(a);
            if (name=="atan2")  return atan2(a,b);
            if (name=="sinh")   return sinh(a);
            if (name=="cosh")   return cosh(a);
            if (name=="tanh")   return tanh(a);
            if (name=="asinh")  return asinh(a);
            if (name=="acosh")  return acosh(a);
            if (name=="atanh")  return atanh(a);
            if (name=="sqrt")   return sqrt(a);
            if (name=="cbrt")   return cbrt(a);
            if (name=="exp")    return exp(a);
            if (name=="log"||name=="ln") return log(a);
            if (name=="log2")   return log2(a);
            if (name=="log10")  return log10(a);
            if (name=="abs")    return fabs(a);
            if (name=="ceil")   return ceil(a);
            if (name=="floor")  return floor(a);
            if (name=="round")  return round(a);
            if (name=="sign")   return (a>0)-(a<0);
            if (name=="pow")    return pow(a,b);
            if (name=="max")    return hasB ? max(a,b) : a;
            if (name=="min")    return hasB ? min(a,b) : a;
            if (name=="mod")    return fmod(a,b);
            if (name=="gcd")    return __gcd((long long)a,(long long)b);
            if (name=="lcm")    return (long long)a / __gcd((long long)a,(long long)b) * (long long)b;
            if (name=="deg")    return a * 180.0 / PI;   // rad->deg
            if (name=="rad")    return a * PI / 180.0;   // deg->rad
            if (name=="fact"||name=="factorial") {
                long long n = (long long)round(a); long long r=1;
                for(long long i=2;i<=n;i++) r*=i; return (double)r;
            }
            if (name=="C"||name=="nCr") { // combinations
                long long n=(long long)a,k=(long long)b;
                if(k>n) return 0;
                if(k==0||k==n) return 1;
                double r=1; for(long long i=0;i<k;i++) r=r*(n-i)/(i+1);
                return round(r);
            }
            if (name=="P"||name=="nPr") { // permutations
                long long n=(long long)a,k=(long long)b;
                double r=1; for(long long i=0;i<k;i++) r*=(n-i);
                return r;
            }
            throw runtime_error("Unknown function: " + name);
        }
        throw runtime_error("Unknown identifier: " + name);
    }

    double parsePow() {
        double base = parseAtom();
        skipSpaces();
        if (pos < src.size() && src[pos] == '^') { pos++; double exp_ = parsePow(); return pow(base, exp_); } // right-assoc
        return base;
    }

    double parseTerm() {
        double left = parsePow();
        while (true) {
            skipSpaces();
            if (pos < src.size() && src[pos] == '*') { pos++; left *= parsePow(); }
            else if (pos < src.size() && src[pos] == '/') { pos++; double r=parsePow(); if(r==0) throw runtime_error("Division by zero"); left /= r; }
            else if (pos < src.size() && src[pos] == '%') { pos++; left = fmod(left, parsePow()); }
            else break;
        }
        return left;
    }

    double parseExpr() {
        double left = parseTerm();
        while (true) {
            skipSpaces();
            if (pos < src.size() && src[pos] == '+') { pos++; left += parseTerm(); }
            else if (pos < src.size() && src[pos] == '-') { pos++; left -= parseTerm(); }
            else break;
        }
        return left;
    }

    double parse() { double v = parseExpr(); skipSpaces(); return v; }
};

double evalExpr(const string& s) {
    Parser p(s);
    return p.parse();
}

// ─────────────────────────────────────────────────────────────
//  NUMERICAL CALCULUS
// ─────────────────────────────────────────────────────────────
// Numerical derivative (5-point stencil, very accurate)
double derivative(function<double(double)> f, double x, int order=1) {
    const double h = 1e-5;
    if (order == 1)
        return (-f(x+2*h) + 8*f(x+h) - 8*f(x-h) + f(x-2*h)) / (12*h);
    if (order == 2)
        return (-f(x+2*h) + 16*f(x+h) - 30*f(x) + 16*f(x-h) - f(x-2*h)) / (12*h*h);
    if (order == 3)
        return (f(x+2*h) - 2*f(x+h) + 2*f(x-h) - f(x-2*h)) / (2*h*h*h);
    if (order == 4)
        return (f(x+2*h) - 4*f(x+h) + 6*f(x) - 4*f(x-h) + f(x-2*h)) / (h*h*h*h);
    return 0;
}

// Adaptive Simpson integration (very accurate)
double simpsonStep(function<double(double)> f, double a, double b) {
    double m = (a+b)/2;
    return (b-a)/6.0 * (f(a) + 4*f(m) + f(b));
}

double adaptiveSimpson(function<double(double)> f, double a, double b, double eps, double whole, int depth) {
    double m = (a+b)/2;
    double left  = simpsonStep(f,a,m);
    double right = simpsonStep(f,m,b);
    if (depth >= 50) return left+right;
    if (fabs(left+right-whole) <= 15*eps) return left+right + (left+right-whole)/15.0;
    return adaptiveSimpson(f,a,m,eps/2,left,depth+1) + adaptiveSimpson(f,m,b,eps/2,right,depth+1);
}

double integrate(function<double(double)> f, double a, double b, double eps=1e-9) {
    return adaptiveSimpson(f, a, b, eps, simpsonStep(f,a,b), 0);
}

// Limit (numerical, from both sides)
struct LimitResult { double value; bool exists; bool leftOnly; bool rightOnly; };
LimitResult computeLimit(function<double(double)> f, double x0) {
    vector<double> eps = {1e-3,1e-4,1e-5,1e-6,1e-7,1e-8,1e-9,1e-10};
    double lv = f(x0 - eps.back()), rv = f(x0 + eps.back());
    bool okL = isfinite(lv), okR = isfinite(rv);
    if (!okL && !okR) return {0,false,false,false};
    if (!okL) return {rv,false,false,true};
    if (!okR) return {lv,false,true,false};
    if (fabs(lv-rv) < 1e-6) return {(lv+rv)/2.0, true, false, false};
    return {0, false, false, false};
}

// Newton-Raphson root finding
double findRoot(function<double(double)> f, double x, double eps=1e-12, int maxIter=1000) {
    for (int i=0;i<maxIter;i++) {
        double fx = f(x);
        double dfx = derivative(f,x);
        if (fabs(dfx) < 1e-15) break;
        double xn = x - fx/dfx;
        if (fabs(xn-x) < eps) return xn;
        x = xn;
    }
    return x;
}

// Taylor series expansion (prints first N terms)
void taylorSeries(function<double(double)> f, double a, int terms) {
    cout << "\nTaylor series around x = " << a << ":\n";
    cout << "  f(x) ≈ ";
    double factorial = 1;
    bool first = true;
    for (int n=0; n<terms; n++) {
        if (n>0) factorial *= n;
        // nth derivative at a (finite difference)
        double coeff;
        if (n==0) coeff = f(a);
        else {
            // use high-order finite differences for nth derivative
            // Stirling coefficients
            double h=1e-3, sum=0;
            for (int k=0; k<=n; k++) {
                double binom=1;
                for(int j=0;j<k;j++) binom=binom*(n-j)/(j+1);
                sum += (k%2==0?1:-1) * binom * f(a + (n-2.0*k)*h/2.0);
            }
            coeff = sum / pow(h,n);
        }
        double c = coeff / factorial;
        if (fabs(c) < 1e-12) continue;
        if (!first) cout << (c>=0?" + ":" - ");
        if (!first) c = fabs(c);
        first = false;
        cout << fixed << setprecision(6) << c;
        if (n==1) cout << "(x-" << a << ")";
        else if (n>1) cout << "(x-" << a << ")^" << n;
    }
    cout << " + ...\n";
}

// ─────────────────────────────────────────────────────────────
//  LINEAR ALGEBRA
// ─────────────────────────────────────────────────────────────
using Mat = vector<vector<double>>;

Mat matMul(const Mat& A, const Mat& B) {
    int n=A.size(), m=B[0].size(), k=B.size();
    Mat C(n, vector<double>(m,0));
    for(int i=0;i<n;i++) for(int j=0;j<m;j++) for(int l=0;l<k;l++) C[i][j]+=A[i][l]*B[l][j];
    return C;
}

double determinant(Mat A) {
    int n=A.size(); double det=1;
    for(int col=0;col<n;col++){
        int pivot=col;
        for(int row=col+1;row<n;row++) if(fabs(A[row][col])>fabs(A[pivot][col])) pivot=row;
        if(pivot!=col){ swap(A[pivot],A[col]); det*=-1; }
        if(fabs(A[col][col])<1e-12) return 0;
        det*=A[col][col];
        for(int row=col+1;row<n;row++){
            double f=A[row][col]/A[col][col];
            for(int j=col;j<n;j++) A[row][j]-=f*A[col][j];
        }
    }
    return det;
}

Mat inverse(Mat A) {
    int n=A.size();
    Mat aug(n, vector<double>(2*n,0));
    for(int i=0;i<n;i++){ for(int j=0;j<n;j++) aug[i][j]=A[i][j]; aug[i][n+i]=1; }
    for(int col=0;col<n;col++){
        int pivot=col;
        for(int r=col+1;r<n;r++) if(fabs(aug[r][col])>fabs(aug[pivot][col])) pivot=r;
        swap(aug[pivot],aug[col]);
        double d=aug[col][col];
        if(fabs(d)<1e-12) throw runtime_error("Matrix not invertible");
        for(int j=0;j<2*n;j++) aug[col][j]/=d;
        for(int r=0;r<n;r++) if(r!=col){ double f=aug[r][col]; for(int j=0;j<2*n;j++) aug[r][j]-=f*aug[col][j]; }
    }
    Mat inv(n, vector<double>(n));
    for(int i=0;i<n;i++) for(int j=0;j<n;j++) inv[i][j]=aug[i][n+j];
    return inv;
}

// Gauss-Jordan to solve Ax=b
vector<double> solveLinear(Mat A, vector<double> b) {
    int n=A.size();
    Mat aug(n, vector<double>(n+1));
    for(int i=0;i<n;i++){ for(int j=0;j<n;j++) aug[i][j]=A[i][j]; aug[i][n]=b[i]; }
    for(int col=0;col<n;col++){
        int pivot=col;
        for(int r=col+1;r<n;r++) if(fabs(aug[r][col])>fabs(aug[pivot][col])) pivot=r;
        swap(aug[pivot],aug[col]);
        double d=aug[col][col];
        if(fabs(d)<1e-12) throw runtime_error("No unique solution");
        for(int j=col;j<=n;j++) aug[col][j]/=d;
        for(int r=0;r<n;r++) if(r!=col){ double f=aug[r][col]; for(int j=col;j<=n;j++) aug[r][j]-=f*aug[col][j]; }
    }
    vector<double> x(n);
    for(int i=0;i<n;i++) x[i]=aug[i][n];
    return x;
}

void printMatrix(const Mat& M, const string& label="") {
    if(!label.empty()) cout << "\n" << label << ":\n";
    for(auto& row : M){
        cout << "  [ ";
        for(double v : row) cout << setw(10) << fixed << setprecision(4) << v << " ";
        cout << "]\n";
    }
}

// ─────────────────────────────────────────────────────────────
//  STATISTICS
// ─────────────────────────────────────────────────────────────
struct Stats {
    double mean, median, mode, variance, stddev, min, max, range, q1, q3, iqr, skewness, kurtosis;
};

Stats computeStats(vector<double> data) {
    if(data.empty()) throw runtime_error("Empty dataset");
    sort(data.begin(), data.end());
    int n = data.size();
    Stats s;
    s.min = data.front(); s.max = data.back(); s.range = s.max - s.min;
    s.mean = accumulate(data.begin(),data.end(),0.0)/n;
    s.median = (n%2==0)?(data[n/2-1]+data[n/2])/2.0:data[n/2];
    // mode
    map<double,int> freq;
    for(double x:data) freq[x]++;
    s.mode = max_element(freq.begin(),freq.end(),[](auto&a,auto&b){return a.second<b.second;})->first;
    // variance (sample)
    double var=0; for(double x:data) var+=pow(x-s.mean,2); var/=(n-1);
    s.variance=var; s.stddev=sqrt(var);
    // quartiles
    s.q1=(n%4==0)?(data[n/4-1]+data[n/4])/2.0:data[n/4];
    s.q3=(n*3%4==0)?(data[3*n/4-1]+data[3*n/4])/2.0:data[3*n/4];
    s.iqr=s.q3-s.q1;
    // skewness (Fisher)
    double sk=0; for(double x:data) sk+=pow((x-s.mean)/s.stddev,3);
    s.skewness=sk*n/((n-1.0)*(n-2.0));
    // excess kurtosis
    double ku=0; for(double x:data) ku+=pow((x-s.mean)/s.stddev,4);
    s.kurtosis=ku*(double)n*(n+1)/((n-1.0)*(n-2.0)*(n-3.0)) - 3.0*(n-1.0)*(n-1.0)/((n-2.0)*(n-3.0));
    return s;
}

// Linear regression
struct Regression { double slope, intercept, r2; };
Regression linearRegression(const vector<double>& x, const vector<double>& y) {
    int n=x.size();
    double sx=0,sy=0,sxy=0,sx2=0,sy2=0;
    for(int i=0;i<n;i++){ sx+=x[i]; sy+=y[i]; sxy+=x[i]*y[i]; sx2+=x[i]*x[i]; sy2+=y[i]*y[i]; }
    double slope=(n*sxy-sx*sy)/(n*sx2-sx*sx);
    double intercept=(sy-slope*sx)/n;
    double sse=0,sst=0,ym=sy/n;
    for(int i=0;i<n;i++){ double pred=slope*x[i]+intercept; sse+=pow(y[i]-pred,2); sst+=pow(y[i]-ym,2); }
    double r2=1-sse/sst;
    return {slope,intercept,r2};
}

// ─────────────────────────────────────────────────────────────
//  NUMBER THEORY
// ─────────────────────────────────────────────────────────────
bool isPrime(long long n) {
    if (n < 2) return false;
    if (n==2||n==3) return true;
    if (n%2==0||n%3==0) return false;
    for (long long i=5; i*i<=n; i+=6) if(n%i==0||n%(i+2)==0) return false;
    return true;
}

vector<long long> primeFactors(long long n) {
    vector<long long> factors;
    for(long long d=2; d*d<=n; d++) while(n%d==0){ factors.push_back(d); n/=d; }
    if(n>1) factors.push_back(n);
    return factors;
}

vector<long long> primesUpTo(long long n) {
    vector<bool> sieve(n+1,true); sieve[0]=sieve[1]=false;
    for(long long i=2;i*i<=n;i++) if(sieve[i]) for(long long j=i*i;j<=n;j+=i) sieve[j]=false;
    vector<long long> primes;
    for(long long i=2;i<=n;i++) if(sieve[i]) primes.push_back(i);
    return primes;
}

long long eulerTotient(long long n) {
    long long result=n;
    for(long long p=2;p*p<=n;p++) if(n%p==0){ while(n%p==0) n/=p; result-=result/p; }
    if(n>1) result-=result/n;
    return result;
}

long long modPow(long long base, long long exp, long long mod) {
    long long result=1; base%=mod;
    while(exp>0){ if(exp%2==1) result=result*base%mod; exp/=2; base=base*base%mod; }
    return result;
}

// Extended GCD  -> ax + by = gcd(a,b)
long long extGCD(long long a, long long b, long long& x, long long& y) {
    if(b==0){ x=1;y=0;return a; }
    long long x1,y1,g=extGCD(b,a%b,x1,y1);
    x=y1; y=x1-(a/b)*y1; return g;
}

// ─────────────────────────────────────────────────────────────
//  COMPLEX NUMBERS
// ─────────────────────────────────────────────────────────────
void complexOps(cd a, cd b) {
    cout << fixed << setprecision(6);
    cout << "\n--- Complex Number Operations ---\n";
    cout << "  A = " << a.real() << (a.imag()>=0?"+":"") << a.imag() << "i\n";
    cout << "  B = " << b.real() << (b.imag()>=0?"+":"") << b.imag() << "i\n";
    cout << "  A + B = " << (a+b) << "\n";
    cout << "  A - B = " << (a-b) << "\n";
    cout << "  A * B = " << (a*b) << "\n";
    cout << "  A / B = " << (a/b) << "\n";
    cout << "  |A|   = " << abs(a) << "\n";
    cout << "  arg(A)= " << arg(a) << " rad  (" << arg(a)*180/PI << " deg)\n";
    cout << "  conj(A)= " << conj(a) << "\n";
    cout << "  A^2   = " << pow(a,2.0) << "\n";
    cout << "  sqrt(A)= " << sqrt(a) << "\n";
    cout << "  exp(A) = " << exp(a) << "\n";
    cout << "  ln(A)  = " << log(a) << "\n";
    cout << "  sin(A) = " << sin(a) << "\n";
    cout << "  cos(A) = " << cos(a) << "\n";
}

// ─────────────────────────────────────────────────────────────
//  POLYNOMIALS
// ─────────────────────────────────────────────────────────────
// poly stored as coefficients [a0, a1, a2, ...] => a0 + a1*x + a2*x^2 ...
using Poly = vector<double>;

double polyEval(const Poly& p, double x) {
    double result=0, xp=1;
    for(double c : p){ result+=c*xp; xp*=x; }
    return result;
}

Poly polyDerive(const Poly& p) {
    if(p.size()<=1) return {0};
    Poly d(p.size()-1);
    for(int i=1;i<(int)p.size();i++) d[i-1]=p[i]*i;
    return d;
}

Poly polyAdd(const Poly& a, const Poly& b) {
    Poly res(max(a.size(),b.size()),0);
    for(int i=0;i<(int)a.size();i++) res[i]+=a[i];
    for(int i=0;i<(int)b.size();i++) res[i]+=b[i];
    return res;
}

Poly polyMul(const Poly& a, const Poly& b) {
    Poly res(a.size()+b.size()-1,0);
    for(int i=0;i<(int)a.size();i++) for(int j=0;j<(int)b.size();j++) res[i+j]+=a[i]*b[j];
    return res;
}

void printPoly(const Poly& p, const string& var="x") {
    bool first=true;
    for(int i=p.size()-1;i>=0;i--){
        if(fabs(p[i])<1e-12) continue;
        if(!first) cout<<(p[i]>=0?" + ":" - ");
        if(!first) cout<<fabs(p[i]);
        else cout<<p[i];
        if(i==1) cout<<var;
        else if(i>1) cout<<var<<"^"<<i;
        first=false;
    }
    if(first) cout<<"0";
}

// Roots of quadratic ax^2+bx+c=0
void quadraticRoots(double a, double b, double c) {
    cout << "\nQuadratic: " << a << "x² + " << b << "x + " << c << " = 0\n";
    double disc = b*b - 4*a*c;
    cout << "  Discriminant = " << disc << "\n";
    if(disc > 0){
        double r1=(-b+sqrt(disc))/(2*a), r2=(-b-sqrt(disc))/(2*a);
        cout << "  Two real roots: x₁ = " << r1 << ", x₂ = " << r2 << "\n";
    } else if(disc==0){
        cout << "  One repeated root: x = " << -b/(2*a) << "\n";
    } else {
        cout << "  Two complex roots: x = " << -b/(2*a) << " ± " << sqrt(-disc)/(2*a) << "i\n";
    }
}

// Roots of cubic ax^3+bx^2+cx+d=0 (Cardano)
void cubicRoots(double a, double b, double c, double d) {
    cout << "\nCubic: " << a << "x³ + " << b << "x² + " << c << "x + " << d << " = 0\n";
    b/=a; c/=a; d/=a; // normalize
    double p = c - b*b/3;
    double q = 2*b*b*b/27 - b*c/3 + d;
    double disc = q*q/4 + p*p*p/27;
    double shift = b/3;
    if(disc > 0){
        double u=cbrt(-q/2+sqrt(disc)), v=cbrt(-q/2-sqrt(disc));
        double x1=(u+v)-shift;
        double re=-(u+v)/2-shift, im=sqrt(3)/2*(u-v);
        cout << "  x₁ = " << x1 << "\n";
        cout << "  x₂ = " << re << " + " << im << "i\n";
        cout << "  x₃ = " << re << " - " << im << "i\n";
    } else if(fabs(disc)<1e-12){
        double u=cbrt(-q/2);
        cout << "  x₁ = " << 2*u-shift << "\n";
        cout << "  x₂ = x₃ = " << -u-shift << "\n";
    } else {
        double r=sqrt(-p*p*p/27), theta=acos(-q/(2*r))/3;
        double mag=2*cbrt(r);
        for(int k=0;k<3;k++) cout << "  x" << k+1 << " = " << mag*cos(theta+2*PI*k/3)-shift << "\n";
    }
}

// ─────────────────────────────────────────────────────────────
//  SERIES & SEQUENCES
// ─────────────────────────────────────────────────────────────
void arithmeticSeries(double a, double d, int n) {
    cout << "\nArithmetic Series: a=" << a << ", d=" << d << ", n=" << n << "\n";
    cout << "  nth term    = " << a+(n-1)*d << "\n";
    cout << "  Sum(1..n)   = " << n/2.0*(2*a+(n-1)*d) << "\n";
    cout << "  First 10    = ";
    for(int i=0;i<min(10,n);i++) cout << a+i*d << " ";
    cout << "\n";
}

void geometricSeries(double a, double r, int n) {
    cout << "\nGeometric Series: a=" << a << ", r=" << r << ", n=" << n << "\n";
    cout << "  nth term    = " << a*pow(r,n-1) << "\n";
    if(fabs(r)!=1) cout << "  Sum(1..n)   = " << a*(pow(r,n)-1)/(r-1) << "\n";
    if(fabs(r)<1)  cout << "  Sum(inf)    = " << a/(1-r) << "\n";
    cout << "  First 10    = ";
    for(int i=0;i<min(10,n);i++) cout << a*pow(r,i) << " ";
    cout << "\n";
}
// ─────────────────────────────────────────────────────────────
//  POLAR MATH LOGIC (MUST BE ABOVE MAIN)
// ─────────────────────────────────────────────────────────────
struct PolarPoint { 
    double r, theta_rad, theta_deg; 
};

PolarPoint toPolar(double x, double y) {
    double r = sqrt(x*x + y*y);
    double theta = atan2(y, x); 
    return {r, theta, theta * 180.0 / PI};
}
// 20. Probability Logic
double normalPDF(double x, double mu, double sigma) {
    return (1.0 / (sigma * sqrt(2.0 * PI))) * exp(-0.5 * pow((x - mu) / sigma, 2.0));
}
double normalCDF(double x, double mu, double sigma) {
    return 0.5 * (1.0 + erf((x - mu) / (sigma * sqrt(2.0))));
}

// 21 & 22. Physics Logic
double energyFromMass(double m) { return m * C_LIGHT * C_LIGHT; }
double photonEnergy(double freq) { return H_PLANCK * freq; }
           

} else if(choice == 0) {
                cout << "  Goodbye!\n";
                break; 

            } else {
                cout << "  Invalid choice or not yet implemented.\n";
            }

        } catch(exception& ex) {
            cout << "  ERROR: " << ex.what() << "\n";
        }
        cout << "\n------------------------------------------------\n";
    }
    return 0;
}
}}
// ─────────────────────────────────────────────────────────────
//  DISPLAY MENU
// ─────────────────────────────────────────────────────────────
void printMenu() {
    cout << R"(
╔══════════════════════════════════════════════════════════════════╗
║              ULTIMATE C++ CALCULATOR — COMMAND MENU              ║
╠══════════════════════════════════════════════════════════════════╣
║  1.  EXPRESSION EVALUATOR     — eval any math expression         ║
║  2.  DERIVATIVE               — f'(x) at a point, any order      ║
║  3.  DEFINITE INTEGRAL        — ∫[a,b] f(x) dx                   ║
║  4.  LIMIT                    — lim(x→c) f(x)                    ║
║  5.  ROOT FINDER (Newton)     — find root near x₀                ║
║  6.  TAYLOR SERIES            — expand f(x) around a             ║
║  7.  STATISTICS               — mean, median, stddev, etc.       ║
║  8.  LINEAR REGRESSION        — best fit line for data           ║
║  9.  MATRIX OPERATIONS        — det, inverse, multiply, solve    ║
║  10. QUADRATIC ROOTS          — ax²+bx+c = 0                     ║
║  11. CUBIC ROOTS              — ax³+bx²+cx+d = 0                 ║
║  12. COMPLEX NUMBERS          — add/mul/div/pow/trig complex     ║
║  13. NUMBER THEORY            — prime, factors, gcd, totient     ║
║  14. ARITHMETIC SERIES        — sum & terms of A.P.              ║
║  15. GEOMETRIC SERIES         — sum & terms of G.P.              ║
║  16. COMBINATIONS & PERMUTATIONS                                 ║
║  17. UNIT CONVERTER           — length/mass/temp/angle           ║
║  18. TRIG FULL SOLVE                                             ║    
║  19. POLAR COORDINATES - — convert between polar/rectangular     ║
║  20. PROBABILITIY DISTRIBUTIONS                                  ║
║  0.  QUIT                                                        ║
╚══════════════════════════════════════════════════════════════════╝
)";
    cout << "Enter choice: ";
}

// ─────────────────────────────────────────────────────────────
//  INPUT HELPERS
// ─────────────────────────────────────────────────────────────
string getLine(const string& prompt="") {
    if(!prompt.empty()) cout << prompt;
    string s; getline(cin,s); return s;
}
double getDouble(const string& prompt="") {
    cout << prompt; double v; cin>>v; cin.ignore(); return v;
}
int getInt(const string& prompt="") {
    cout << prompt; int v; cin>>v; cin.ignore(); return v;
}

vector<double> getDataset() {
    cout << "  Enter values separated by spaces: ";
    string line; getline(cin,line);
    istringstream iss(line); vector<double> data;
    double x; while(iss>>x) data.push_back(x);
    return data;
}

// ─────────────────────────────────────────────────────────────
//  MAIN
// ─────────────────────────────────────────────────────────────
int main() {
    cout << fixed << setprecision(8);
    cout << "\n  ULTIMATE C++ CALCULATOR  |  All answers computed directly\n";
    cout << "  Tip: expressions support sin,cos,tan,log,sqrt,exp,abs,pi,e,^ etc.\n\n";

    int choice;
    while(true) {
        printMenu();
        if(!(cin>>choice)) break;
        cin.ignore();
        cout << "\n";

        try {
            // ─────────── 1. EXPRESSION ───────────
            if(choice==1){
                string expr = getLine("  Enter expression: ");
                double result = evalExpr(expr);
                cout << "  = " << result << "\n";
                // also show hex / scientific
                cout << "  (scientific) = " << scientific << result << "\n" << fixed;
                if(result==(long long)result)
                    cout << "  (integer)    = " << (long long)result << "\n";

            // ─────────── 2. DERIVATIVE ───────────
            } else if(choice==2){
                string fn  = getLine("  Enter f(x): ");
                double x   = evalExpr(getLine("  At x = "));
                int order  = getInt("  Derivative order (1-4): ");
                auto f = [&](double xv){ return evalExpr([&](){
                    // replace x symbol — inline substitution via parser variable
                    string s=fn; string xstr=to_string(xv);
                    // Simple: build expression with x replaced
                    string result2;
                    for(size_t i=0;i<s.size();i++){
                        if(s[i]=='x'&&(i==0||!isalpha(s[i-1]))&&(i+1>=s.size()||!isalpha(s[i+1])))
                            result2+="("+xstr+")";
                        else result2+=s[i];
                    }
                    return result2;
                }()); };
                double d = derivative(f, x, order);
                cout << "  f" << string(order,'\'') << "(" << x << ") = " << d << "\n";

            // ─────────── 3. INTEGRAL ───────────
            } else if(choice==3){
                string fn = getLine("  Enter f(x): ");
                double a  = evalExpr(getLine("  From a = "));
                double b  = evalExpr(getLine("  To   b = "));
                auto f = [&](double xv){
                    string s=fn, result2;
                    for(size_t i=0;i<s.size();i++){
                        if(s[i]=='x'&&(i==0||!isalpha(s[i-1]))&&(i+1>=s.size()||!isalpha(s[i+1])))
                            result2+="("+to_string(xv)+")";
                        else result2+=s[i];
                    }
                    return evalExpr(result2);
                };
                double result = integrate(f, a, b);
                cout << "  ∫[" << a << " to " << b << "] " << fn << " dx = " << result << "\n";

            // ─────────── 4. LIMIT ───────────
            } else if(choice==4){
                string fn = getLine("  Enter f(x): ");
                double x0 = evalExpr(getLine("  As x → "));
                auto f = [&](double xv){
                    string s=fn, r;
                    for(size_t i=0;i<s.size();i++){
                        if(s[i]=='x'&&(i==0||!isalpha(s[i-1]))&&(i+1>=s.size()||!isalpha(s[i+1])))
                            r+="("+to_string(xv)+")";
                        else r+=s[i];
                    }
                    try{ return evalExpr(r); } catch(...){ return numeric_limits<double>::quiet_NaN(); }
                };
                auto L = computeLimit(f, x0);
                if(L.exists) cout << "  lim(x→" << x0 << ") = " << L.value << "\n";
                else if(L.leftOnly)  cout << "  Left-hand limit  = " << L.value << "  (right undefined)\n";
                else if(L.rightOnly) cout << "  Right-hand limit = " << L.value << "  (left undefined)\n";
                else cout << "  Limit does NOT exist at x = " << x0 << "\n";

            // ─────────── 5. ROOT FINDER ───────────
            } else if(choice==5){
                string fn = getLine("  Enter f(x): ");
                double x0 = evalExpr(getLine("  Initial guess x₀ = "));
                auto f = [&](double xv){
                    string s=fn, r;
                    for(size_t i=0;i<s.size();i++){
                        if(s[i]=='x'&&(i==0||!isalpha(s[i-1]))&&(i+1>=s.size()||!isalpha(s[i+1])))
                            r+="("+to_string(xv)+")";
                        else r+=s[i];
                    }
                    return evalExpr(r);
                };
                double root = findRoot(f, x0);
                cout << "  Root ≈ " << root << "\n";
                cout << "  Verification: f(" << root << ") = " << f(root) << "\n";

            // ─────────── 6. TAYLOR ───────────
            } else if(choice==6){
                string fn = getLine("  Enter f(x): ");
                double a = evalExpr(getLine("  Expand around a = "));
                int terms = getInt("  Number of terms: ");
                auto f = [&](double xv){
                    string s=fn, r;
                    for(size_t i=0;i<s.size();i++){
                        if(s[i]=='x'&&(i==0||!isalpha(s[i-1]))&&(i+1>=s.size()||!isalpha(s[i+1])))
                            r+="("+to_string(xv)+")";
                        else r+=s[i];
                    }
                    return evalExpr(r);
                };
                taylorSeries(f, a, terms);

            // ─────────── 7. STATISTICS ───────────
            } else if(choice==7){
                cout << "  Enter dataset:\n";
                auto data = getDataset();
                auto s = computeStats(data);
                cout << "\n  n           = " << data.size() << "\n";
                cout << "  Mean        = " << s.mean << "\n";
                cout << "  Median      = " << s.median << "\n";
                cout << "  Mode        = " << s.mode << "\n";
                cout << "  Std Dev     = " << s.stddev << "\n";
                cout << "  Variance    = " << s.variance << "\n";
                cout << "  Min         = " << s.min << "\n";
                cout << "  Max         = " << s.max << "\n";
                cout << "  Range       = " << s.range << "\n";
                cout << "  Q1          = " << s.q1 << "\n";
                cout << "  Q3          = " << s.q3 << "\n";
                cout << "  IQR         = " << s.iqr << "\n";
                cout << "  Skewness    = " << s.skewness << "\n";
                cout << "  Kurtosis    = " << s.kurtosis << "\n";

            // ─────────── 8. LINEAR REGRESSION ───────────
            } else if(choice==8){
                cout << "  Enter x values:\n"; auto xv = getDataset();
                cout << "  Enter y values:\n"; auto yv = getDataset();
                if(xv.size()!=yv.size()) throw runtime_error("x and y must have same size");
                auto r = linearRegression(xv,yv);
                cout << "\n  y = " << r.slope << "x + " << r.intercept << "\n";
                cout << "  Slope     = " << r.slope << "\n";
                cout << "  Intercept = " << r.intercept << "\n";
                cout << "  R²        = " << r.r2 << "\n";
                cout << "  R (corr)  = " << sqrt(r.r2) * (r.slope>=0?1:-1) << "\n";

            // ─────────── 9. MATRIX ───────────
            } else if(choice==9){
                int n = getInt("  Matrix size n (for n×n): ");
                cout << "  Choose: (1) Determinant  (2) Inverse  (3) Multiply  (4) Solve Ax=b\n  > ";
                int sub; cin>>sub; cin.ignore();
                auto readMat = [&](const string& label) {
                    cout << "  Enter " << label << " row by row (" << n << " values each):\n";
                    Mat M(n, vector<double>(n));
                    for(int i=0;i<n;i++){
                        cout << "    Row " << i+1 << ": ";
                        for(int j=0;j<n;j++) cin>>M[i][j];
                    } cin.ignore(); return M;
                };
                if(sub==1){
                    auto A = readMat("A"); printMatrix(A,"A");
                    cout << "\n  det(A) = " << determinant(A) << "\n";
                } else if(sub==2){
                    auto A = readMat("A"); printMatrix(A,"A");
                    auto inv = inverse(A); printMatrix(inv,"A⁻¹");
                } else if(sub==3){
                    auto A = readMat("A"), B = readMat("B");
                    printMatrix(A,"A"); printMatrix(B,"B");
                    printMatrix(matMul(A,B),"A × B");
                } else if(sub==4){
                    auto A = readMat("A");
                    cout << "  Enter b vector (" << n << " values): ";
                    vector<double> b(n); for(auto& v:b) cin>>v; cin.ignore();
                    auto x = solveLinear(A, b);
                    cout << "\n  Solution x:\n";
                    for(int i=0;i<n;i++) cout << "    x[" << i+1 << "] = " << x[i] << "\n";
                }

            // ─────────── 10. QUADRATIC ───────────
            } else if(choice==10){
                double a=getDouble("  a = "), b=getDouble("  b = "), c=getDouble("  c = ");
                quadraticRoots(a,b,c);

            // ─────────── 11. CUBIC ───────────
            } else if(choice==11){
                double a=getDouble("  a = "),b=getDouble("  b = "),c=getDouble("  c = "),d=getDouble("  d = ");
                cubicRoots(a,b,c,d);

            // ─────────── 12. COMPLEX ───────────
            } else if(choice==12){
                double ar=getDouble("  A real part: "), ai=getDouble("  A imaginary part: ");
                double br=getDouble("  B real part: "), bi=getDouble("  B imaginary part: ");
                complexOps({ar,ai},{br,bi});

            // ─────────── 13. NUMBER THEORY ───────────
            } else if(choice==13){
                long long n = (long long)evalExpr(getLine("  Enter n: "));
                cout << "\n  n              = " << n << "\n";
                cout << "  isPrime        = " << (isPrime(n)?"YES":"NO") << "\n";
                cout << "  Prime factors  = ";
                for(auto p:primeFactors(n)) cout<<p<<" "; cout<<"\n";
                cout << "  Euler totient  = " << eulerTotient(n) << "\n";
                // primes up to n (if reasonable)
                if(n<=10000){
                    auto primes=primesUpTo(n);
                    cout << "  # primes ≤ n   = " << primes.size() << "\n";
                }
                long long m = (long long)evalExpr(getLine("  Enter m (for gcd/lcm/mod): "));
                long long g = __gcd(n,m);
                cout << "  gcd(" << n << "," << m << ") = " << g << "\n";
                cout << "  lcm(" << n << "," << m << ") = " << n/__gcd(n,m)*m << "\n";
                long long mod = (long long)evalExpr(getLine("  Modulus for modPow: "));
                long long exp2 = (long long)evalExpr(getLine("  Exponent for modPow: "));
                cout << "  " << n << "^" << exp2 << " mod " << mod << " = " << modPow(n,exp2,mod) << "\n";
                long long x2,y2; extGCD(n,m,x2,y2);
                cout << "  Bezout: " << n << "*" << x2 << " + " << m << "*" << y2 << " = " << g << "\n";

            // ─────────── 14. ARITHMETIC SERIES ───────────
            } else if(choice==14){
                double a=getDouble("  First term a = ");
                double d=getDouble("  Common diff d = ");
                int n=getInt("  # terms n = ");
                arithmeticSeries(a,d,n);

            // ─────────── 15. GEOMETRIC SERIES ───────────
            } else if(choice==15){
                double a=getDouble("  First term a = ");
                double r=getDouble("  Common ratio r = ");
                int n=getInt("  # terms n = ");
                geometricSeries(a,r,n);

            // ─────────── 16. COMBINATIONS ───────────
            } else if(choice==16){
                long long n=(long long)getDouble("  n = ");
                long long k=(long long)getDouble("  k = ");
                // C(n,k) using log to avoid overflow
                double logC=0; for(long long i=0;i<k;i++) logC+=log(n-i)-log(i+1);
                double Cnk=round(exp(logC));
                double Pnk=1; for(long long i=0;i<k;i++) Pnk*=(n-i);
                cout << "\n  C(" << n << "," << k << ") = " << (long long)Cnk << "\n";
                cout << "  P(" << n << "," << k << ") = " << (long long)Pnk << "\n";

            // ─────────── 17. UNIT CONVERTER ───────────
            } else if(choice==17){
                cout << "\n  Categories: (1)Length (2)Mass (3)Temperature (4)Angle (5)Area (6)Volume (7)Speed\n  > ";
                int cat; cin>>cat; cin.ignore();
                double val=getDouble("  Value: ");
                if(cat==1){
                    cout << "  " << val << " m =\n";
                    cout << "    " << val*100 << " cm\n    " << val*1000 << " mm\n";
                    cout << "    " << val*39.3701 << " inches\n    " << val*3.28084 << " feet\n";
                    cout << "    " << val*1.09361 << " yards\n    " << val/1609.34 << " miles\n";
                    cout << "    " << val/1000 << " km\n    " << val/1852 << " nautical miles\n";
                } else if(cat==2){
                    cout << "  " << val << " kg =\n";
                    cout << "    " << val*1000 << " g\n    " << val*1e6 << " mg\n";
                    cout << "    " << val*2.20462 << " lbs\n    " << val*35.274 << " oz\n";
                    cout << "    " << val/1000 << " tonnes\n    " << val*0.000984207 << " long tons\n";
                } else if(cat==3){
                    cout << "  " << val << "°C = " << val*9/5+32 << "°F = " << val+273.15 << " K\n";
                    cout << "  " << val << "°F = " << (val-32)*5/9 << "°C = " << (val-32)*5/9+273.15 << " K\n";
                    cout << "  " << val << " K  = " << val-273.15 << "°C = " << (val-273.15)*9/5+32 << "°F\n";
                } else if(cat==4){
                    cout << "  " << val << " deg = " << val*PI/180 << " rad\n";
                    cout << "  " << val << " rad = " << val*180/PI << " deg\n";
                    cout << "  " << val << " deg = " << val/360 << " rev = " << val*60 << " arcmin = " << val*3600 << " arcsec\n";
                } else if(cat==5){
                    cout << "  " << val << " m² =\n";
                    cout << "    " << val*1e4 << " cm²\n    " << val*10.7639 << " ft²\n";
                    cout << "    " << val*1.19599 << " yd²\n    " << val/4046.86 << " acres\n    " << val/1e6 << " km²\n";
                } else if(cat==6){
                    cout << "  " << val << " L =\n";
                    cout << "    " << val*1000 << " mL\n    " << val/1000 << " m³\n";
                    cout << "    " << val*0.264172 << " US gal\n    " << val*33.814 << " fl oz\n";
                    cout << "    " << val*1.75975 << " UK pints\n";
                } else if(cat==7){
                    cout << "  " << val << " m/s =\n";
                    cout << "    " << val*3.6 << " km/h\n    " << val*2.23694 << " mph\n";
                    cout << "    " << val*1.94384 << " knots\n    " << val/340.29 << " Mach\n";
                }
                // ─────────── 18. TRIG ───────────
            } else if(choice==18){
                double x = getDouble("  Enter angle in degrees: ");
                double r = x * PI / 180.0;
                cout << "\n  Angle: " << x << "° = " << r << " rad\n";
                cout << "  sin(" << x << "°) = " << sin(r) << "\n";
                cout << "  cos(" << x << "°) = " << cos(r) << "\n";
                if(fabs(cos(r))>1e-12) cout << "  tan(" << x << "°) = " << tan(r) << "\n"; else cout << "  tan = undefined\n";
                if(fabs(sin(r))>1e-12) cout << "  csc(" << x << "°) = " << 1/sin(r) << "\n"; else cout << "  csc = undefined\n";
                if(fabs(cos(r))>1e-12) cout << "  sec(" << x << "°) = " << 1/cos(r) << "\n"; else cout << "  sec = undefined\n";
                if(fabs(sin(r))>1e-12) cout << "  cot(" << x << "°) = " << cos(r)/sin(r) << "\n"; else cout << "  cot = undefined\n";
                
                cout << "\n  Inverse trig (input = value, output = degrees):\n";
                double v = getDouble("  Value for asin/acos/atan: ");
                if(v>=-1&&v<=1) cout << "  asin(" << v << ") = " << asin(v)*180/PI << "°\n";
                if(v>=-1&&v<=1) cout << "  acos(" << v << ") = " << acos(v)*180/PI << "°\n";
                cout << "  atan(" << v << ") = " << atan(v)*180/PI << "°\n";
            
            // ─────────── 19. POLAR ───────────
            } else if(choice==19){
                cout << "  1. Cartesian to Polar (x,y -> r,θ)\n";
                cout << "  2. Polar to Cartesian (r,θ -> x,y)\n";
                int sub = getInt("  Choice: ");

                if(sub == 1) {
                    double x = getDouble("  x: ");
                    double y = getDouble("  y: ");
                    
                    auto p = toPolar(x, y); 
                    cout << "  Radius (r): " << p.r << "\n";
                    cout << "  Angle (θ):  " << p.theta_deg << "° (" << p.theta_rad << " rad)\n";
                } 
                else if(sub == 2) {
                    double r = getDouble("  r: ");
                    double t = getDouble("  θ (degrees): ");
                    
                    double x = r * cos(t * PI / 180.0);
                    double y = r * sin(t * PI / 180.0);
                    cout << "  Cartesian Point (x, y): (" << x << ", " << y << ")\n";
                } 
                else {
                    cout << "  Invalid sub-choice.\n";
                }
                // ─────────── 19. POLAR ───────────
            } else if(choice == 19) {
                double x = getDouble("  Enter X: ");
                double y = getDouble("  Enter Y: ");
                PolarPoint p = toPolar(x, y);
                cout << "  r = " << p.r << ", theta = " << p.theta_rad << " rad\n";

            // ─────────── 20. PROBABILITY ───────────
            } else if(choice == 20) {
                double x = getDouble("  Enter x: ");
                double mu = getDouble("  Mean: ");
                double sigma = getDouble("  Std Dev: ");
                cout << "  PDF: " << normalPDF(x, mu, sigma) << "\n";
                cout << "  CDF: " << normalCDF(x, mu, sigma) << "\n";

            // ─────────── 21. RELATIVITY ───────────
            } else if(choice == 21) {
                double m = getDouble("  Enter mass (kg): ");
                cout << "  E = " << scientific << energyFromMass(m) << " J\n" << fixed;

            // ─────────── 22. QUANTUM ───────────
            } else if(choice == 22) {
                double f = getDouble("  Enter freq (Hz): ");
                cout << "  E = " << scientific << photonEnergy(f) << " J\n" << fixed;

            } else if(choice == 0) {
                break;
            }
        } catch(exception& e) {
            cout << "  Error: " << e.what() << "\n";
        }

                // ─────────── 0. EXIT & CATCH-ALL ───────────
            } else if(choice==0){
                cout << "  Goodbye!\n"; break;
            } else {
                cout << "  Invalid choice.\n";
            }

        } catch(exception& ex) {
            cout << "  ERROR: " << ex.what() << "\n";
        }

        cout << "\n";
    }
    return 0;
}