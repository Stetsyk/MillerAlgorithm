#pragma once
#include <boost/multiprecision/cpp_int.hpp>
using boost::multiprecision::cpp_int;

cpp_int readFromHexadecimal(std::string s);

cpp_int binPow(cpp_int a, cpp_int deg, cpp_int modulo);

cpp_int inverse(cpp_int a, cpp_int primeModulo);


/* Takes as input an odd prime p and n < p and returns r
 * such that r * r = n [mod p]. */
cpp_int tonelli_shanks(cpp_int n, cpp_int p) ;
