//
//  ElipticCurveGeometry.hpp
//  Miller_algorithm
//
//  Created by Oleksii Stetsyk on 11/28/19.
//  Copyright © 2019 Oleksii Stetsyk. All rights reserved.
//

#pragma once

#include <stdio.h>
#include "utils.h"
#include <vector>
#include <boost/multiprecision/cpp_int.hpp>
using boost::multiprecision::cpp_int;

struct BigintOnField {
    cpp_int value;
    static cpp_int primeModulo;
    cpp_int toBigint(){
        return value;
    }
    BigintOnField(){}
    
    BigintOnField(long long v) {
        value = v;
    }

    BigintOnField(const std::string &s) {
        value = cpp_int(s);
    }
    
    BigintOnField(const cpp_int & v) {
        value = v;
    }
    
    void operator=(const BigintOnField &v) {
        value = v.value;
    }
};

std::ostream & operator << (std::ostream & stream, const BigintOnField & num);

BigintOnField operator + (const BigintOnField & lhs, const BigintOnField & rhs);
BigintOnField operator - (const BigintOnField & lhs, const BigintOnField & rhs);
BigintOnField operator * (const BigintOnField & lhs, const BigintOnField & rhs);
BigintOnField operator / (const BigintOnField & lhs, const BigintOnField & rhs);
bool operator ==(const BigintOnField & lhs, const BigintOnField & rhs);

BigintOnField binPow(BigintOnField a, BigintOnField deg, BigintOnField modulo);
BigintOnField inverse(BigintOnField a, BigintOnField primeModulo);

struct Point{
    BigintOnField x, y;
};

struct Line{
    BigintOnField a, b, c;
    BigintOnField value(Point p){
        return (a * p.x + b * p.y + c);
    }
};

std::ostream & operator << (std::ostream & stream, const Point & P);

std::ostream & operator << (std::ostream & stream, const Line & L);


bool operator==(const Point & lhs, const Point & rhs);
bool operator!=(const Point & lhs, const Point & rhs);

bool operator != (const BigintOnField & lhs, const BigintOnField & rhs);


class ElipticCurveGeometry{
public:
    BigintOnField pow(BigintOnField a, BigintOnField deg);
    BigintOnField inv(BigintOnField a);
    BigintOnField sum(const BigintOnField & lhs, const BigintOnField & rhs) ;
    BigintOnField sub(const BigintOnField & lhs, const BigintOnField & rhs);
    BigintOnField mul(const BigintOnField & lhs, const BigintOnField & rhs);
    BigintOnField div(const BigintOnField & lhs, const BigintOnField & rhs);

    Point sum(const Point & P, const Point & Q);
    Point sub(const Point & P, const Point & Q);
    Point neg(const Point & P);
    Point mul(Point P, cpp_int deg);
    
    Line lineTroughPoints(const Point & P, const Point & Q);
    Line verticalLine(const Point & P);
    BigintOnField slopeOfLine(const Point & P, const Point & Q){
        if(P.x == Q.x && P.y != Q.y){
            return primeModulo;
        }
        else if(P.x == Q.x && P.y == Q.y){
            return (3 * P.x * P.x + coefA) / (2 * P.y);
        }
        else{
            return (P.y - Q.y) / (P.x - Q.x);
        }
    }

    Point gen();
    
    bool onCurve(Point P){
        BigintOnField first = (P.x * P.x * P.x + getCoefA() * P.x + getCoefB());
        BigintOnField second = (P.y * P.y);
        return first == second;
    }
    
    void setPrimeModulo(BigintOnField primeModulo);
    BigintOnField getPrimeModulo();
    void setCoefA(BigintOnField coefA);
    BigintOnField getCoefA();
    void setCoefB(BigintOnField coefB);
    BigintOnField getCoefB();
    Point getOInf(){
        return OInf;
    }
    
    std::vector<Point> genPoints(int number) {
        Point P = gen();
        std::vector<Point> res;
        res.push_back(P);
        Point Q = P;
        for(int i = 2; i <= number; ++i){
            Q = sum(Q, P);
            res.push_back(Q);
        }
        return res;
    }
    
private:
    BigintOnField primeModulo;
    Point OInf = {0, primeModulo};
    BigintOnField coefA;
    BigintOnField coefB;
};
std::ostream & operator << (std::ostream & stream, ElipticCurveGeometry & el);

