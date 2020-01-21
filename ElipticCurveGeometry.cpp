//
//  ElipticCurveGeometry::.cpp
//  Miller_algorithm
//
//  Created by Oleksii Stetsyk on 11/28/19.
//  Copyright © 2019 Oleksii Stetsyk. All rights reserved.
//

#include "ElipticCurveGeometry.h"
#include <random>
using namespace std;
//BigintOnField operator +

std::ostream & operator << (std::ostream & stream, const BigintOnField & num){
    stream << num.value;
    return stream;
}
    
BigintOnField binPow(BigintOnField a, BigintOnField deg, BigintOnField modulo){
    return BigintOnField(binPow(a.value, deg.value, modulo.value));
}
BigintOnField inverse(BigintOnField a, BigintOnField primeModulo){
    return BigintOnField(inverse(a.value, primeModulo.value));
}

BigintOnField operator + (const BigintOnField & lhs, const BigintOnField & rhs){
    BigintOnField res;
    res.value = lhs.value + rhs.value;
    if(res.value >= res.primeModulo){
        res.value -= res.primeModulo;
    }
    return res;
}

BigintOnField operator - (const BigintOnField & lhs, const BigintOnField & rhs){
    BigintOnField res;
    res.value = lhs.value - rhs.value;
    if(res.value < 0){
        res.value += res.primeModulo;
    }
    return res;
}

BigintOnField operator * (const BigintOnField & lhs, const BigintOnField & rhs){
    BigintOnField res;
    res.value = lhs.value * rhs.value % res.primeModulo;
    return res;
}

BigintOnField operator / (const BigintOnField & lhs, const BigintOnField & rhs){
    BigintOnField res;
    res.value = lhs.value * inverse(rhs.value, res.primeModulo) % res.primeModulo;
    return res;
}

bool operator == (const BigintOnField & lhs, const BigintOnField & rhs){
    return lhs.value == rhs.value;
}
    
bool operator != (const BigintOnField & lhs, const BigintOnField & rhs){
    return lhs.value != rhs.value;
}

std::ostream & operator << (std::ostream & stream, const Point & P){
    stream << P.x << " " << P.y;
    return stream;
}
    
std::ostream & operator << (std::ostream & stream, const Line & l){
    stream << l.a << " " << l.b <<  " " << l.c;
    return stream;
}

bool operator==(const Point & lhs, const Point & rhs){
    return (lhs.x == rhs.x && lhs.y == rhs.y);
}
    
bool operator!=(const Point & lhs, const Point & rhs){
    return !(lhs == rhs);
}

BigintOnField ElipticCurveGeometry::pow(BigintOnField a, BigintOnField deg){
    return binPow(a, deg, primeModulo);
}

BigintOnField ElipticCurveGeometry::inv(BigintOnField a){
    return inverse(a, primeModulo);
}

BigintOnField ElipticCurveGeometry::sum(const BigintOnField & lhs, const BigintOnField & rhs){
    return lhs + rhs;
}

BigintOnField ElipticCurveGeometry::sub(const BigintOnField & lhs, const BigintOnField & rhs){
    return lhs - rhs;
}

BigintOnField ElipticCurveGeometry::mul(const BigintOnField & lhs, const BigintOnField & rhs){
    return lhs * rhs;
}

BigintOnField ElipticCurveGeometry::div(const BigintOnField & lhs, const BigintOnField & rhs){
    return lhs / rhs;
}

Point ElipticCurveGeometry::sum(const Point & P, const Point & Q){
    if(P == OInf || Q == OInf){
        return (P == OInf) ? Q : P;
    }
    if(P.x != Q.x){
        BigintOnField s = div(sub(P.y, Q.y), sub(P.x, Q.x));
        BigintOnField xRes = sub(sub(pow(s, 2), P.x), Q.x);
        BigintOnField yRes = sum(P.y, mul(s, sub(xRes, P.x)));
        yRes = (primeModulo - yRes);
        return {xRes, yRes};
    }
    else{
        if(P.y == Q.y){
            BigintOnField s = div(sum(coefA, mul(BigintOnField(3), pow(P.x, 2))), mul(BigintOnField(2), P.y));
            BigintOnField xRes = sub(pow(s, 2), mul(BigintOnField(2), P.x));
            BigintOnField yRes = sum(P.y, mul(s, sub(xRes, P.x)));
            yRes = sub(primeModulo, yRes);
            return {xRes, yRes};
        }
        else{
            return OInf; // it's O, special point
        }
    }
}
    
Point ElipticCurveGeometry::neg(const Point & P){
    Point res = P;
    res.y = sub(0, res.y);
    return res;
}
    
Point ElipticCurveGeometry::sub(const Point & P, const Point & Q){
    Point res = neg(Q);
    res = sum(res, P);
    return res;
}


Point ElipticCurveGeometry::mul(Point P, cpp_int deg){
    Point res = OInf;
    while(deg > 0){
        if(deg % 2 == 1){
            res = sum(res, P);
        }
        deg = deg / 2;
        P = sum(P, P);
    }
    return res;
}

Point ElipticCurveGeometry::gen(){
//    std::srand(time(0));
    random();
    for(BigintOnField i = random() + 2; ; i = i + 1){
        BigintOnField xRes = i;
        BigintOnField y2 = sum(coefB, sum(mul(coefA, xRes), pow(xRes, 3)));
        BigintOnField yRes = tonelli_shanks(y2.value, primeModulo.value);
        if(yRes != 0){
            if((yRes * yRes) != y2){
                throw std::logic_error("first mistake");
            }
            return {xRes, yRes};
        }
    }
}

void ElipticCurveGeometry::setPrimeModulo(BigintOnField primeModulo){
    BigintOnField::primeModulo = primeModulo.value;
    this->primeModulo = primeModulo;
    OInf = {0, primeModulo};
}
BigintOnField ElipticCurveGeometry::getPrimeModulo(){
    return primeModulo;
}
void ElipticCurveGeometry::setCoefA(BigintOnField coefA){this->coefA = coefA;}
BigintOnField ElipticCurveGeometry::getCoefA(){
    return coefA;
    
}
void ElipticCurveGeometry::setCoefB(BigintOnField coefB){this->coefB = coefB;}
BigintOnField ElipticCurveGeometry::getCoefB(){
    return coefB;
}
    

std::ostream & operator << (std::ostream & stream, ElipticCurveGeometry & el){
    std::cout << "Prime modulo = " << el.getPrimeModulo() << std::endl;
    std::cout << "Coef A " << el.getCoefA() << std::endl;
    std::cout << "Coef B " << el.getCoefB() << std::endl;
    return stream;
}
    
    Line ElipticCurveGeometry::verticalLine(const Point & first){
        Line l;
    l.a = 1;
    l.b = 0;
    l.c = sub(0, first.x);
    return l;
    }
    
Line ElipticCurveGeometry::lineTroughPoints(const Point & first, const Point & second){
    Line l;
    if(first == second){
        // points are equal
        l.b = 1;
        l.a = sub(0, div(sum(coefA, mul(BigintOnField(3), pow(first.x, 2))), mul(BigintOnField(2), first.y)));
        BigintOnField tmp1 = sub(0, pow(first.y, 2));
        BigintOnField tmp2 = mul(mul(BigintOnField(2), coefA), first.x);
        BigintOnField tmp3 = mul(BigintOnField(3), coefB);
        BigintOnField tmp4 = mul(BigintOnField(2), first.y);
        l.c = div(sum(tmp1, sum(tmp2, tmp3)), tmp4);
        l.c = sub(0, l.c);
        return l;
    }
    else{
        if(first.x == second.x){
            // vertical line
            l.a = 1;
            l.b = 0;
            l.c = sub(0, first.x);
            return l;
        }
        else {
            // just regular line
            l.a = sub(0, div(sub(second.y, first.y), sub(second.x, first.x)));
            l.b = 1;
            l.c = div(sub(mul(first.x, second.y), mul(first.y, second.x)), sub(second.x,first.x));
            return l;
        }
    }
}


