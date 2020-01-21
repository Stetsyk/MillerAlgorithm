
#include "utils.h"

using boost::multiprecision::cpp_int;

cpp_int readFromHexadecimal(std::string s){
    if(s.empty()) return cpp_int(0);
    cpp_int result(0);
    for(int i = 0; i < (int)s.size(); ++i){
        int digit = ('0' <= s[i] && s[i] <= '9') ? (s[i] - '0') :  (s[i] - 'A' + 10);
        result = cpp_int(16) * result + digit;
    }
    return result;
}

cpp_int binPow(cpp_int a, cpp_int deg, cpp_int modulo){
    cpp_int res(1);
    while(deg > 0){
        if(deg % 2 == 1){
            res = (res * a) % modulo;
        }
        a = (a * a) % modulo;
        deg = deg / 2;
    }
    return res;
}

cpp_int inverse(cpp_int a, cpp_int primeModulo){
    return binPow(a, primeModulo - 2, primeModulo);
}

cpp_int tonelli_shanks(cpp_int n, cpp_int p) {
  cpp_int s = 0;
  cpp_int q = p - 1;
  while ((q % 2) == 0) { q /= 2; s = s + 1; }
  if (s == 1) {
    cpp_int r = binPow(n, (p+1)/4, p);
    if ((r * r) % p == n) return r;
    return 0;
  }
  // Find the first quadratic non-residue z by brute-force search
  long long z = 1;
    while (binPow(++z, (p-1)/2, p) != p - 1){
        
    }
  cpp_int c = binPow(z, q, p);
  cpp_int r = binPow(n, (q+1)/2, p);
  cpp_int t = binPow(n, q, p);
  cpp_int m = s;
  while (t != 1) {
    cpp_int tt = t;
    cpp_int i = 0;
    while (tt != 1) {
      tt = (tt * tt) % p;
      i = i + 1;
      if (i == m) return 0;
    }
    cpp_int b = binPow(c, binPow(2, m-i-1, p-1), p);
    cpp_int b2 = (b * b) % p;
    r = (r * b) % p;
    t = (t * b2) % p;
    c = b2;
    m = i;
  }
  if ((r * r) % p == n) return r;
  return 0;
}
