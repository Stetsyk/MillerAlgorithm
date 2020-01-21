//
//  main.cpp
//  Miller_algorithm
//
//  Created by Oleksii Stetsyk on 11/28/19.
//  Copyright © 2019 Oleksii Stetsyk. All rights reserved.
//

#include <iostream>
#include "ElipticCurveGeometry.h"
#include "test_runner.h"
#include "utils.h"
#include "profile.h"
#include <boost/multiprecision/cpp_int.hpp>
using namespace boost::multiprecision;
using namespace std;

BigintOnField h(ElipticCurveGeometry el, Point P, Point Q, Point R){
    if ((P == Q && P.y == 0) || (P != Q && P.x == Q.x)){
        return R.x - P.x;
    }
    BigintOnField L = el.slopeOfLine(P, Q);
    BigintOnField p = R.y - P.y - L * (R.x - P.x);
    BigintOnField q = R.x + P.x + Q.x - L * L;
    return p / q;
}

BigintOnField miller(ElipticCurveGeometry el, Point P, Point Q, cpp_int r){
    if(P == Q) return 1;
    Point T = P;
    BigintOnField f = 1;
    bool flag = false;
    for(int i = 1000; i >= 0; --i){
        if(!flag){
            if(r & (cpp_int(1) << i)){
                flag = true;
            }
            continue;
        }
        cpp_int rem = (r & ( cpp_int(1) << i ) );
        f = f * f * h(el, T, T, Q);
        T = el.sum(T, T);
        if (rem) {
            f = f * h(el, T, P, Q);
            T = el.sum(T, P);
        }
    }
    return f;
}

BigintOnField weil_pairing(ElipticCurveGeometry el, Point P, Point Q, cpp_int order, Point S = {0, 0}){
    if(S == Point({0, 0})){
        S = el.gen();
    }
    BigintOnField fpqs = miller(el, P, el.sum(Q, S), order);
    BigintOnField fps = miller(el, P, S, order);
    BigintOnField fqps = miller(el, Q, el.sub(P, S), order);
    BigintOnField fqs = miller(el, Q, el.neg(S), order);
    return fpqs * fqs / fps / fqps;2
}

ElipticCurveGeometry createElipticCurve(){
    ElipticCurveGeometry el;
    el.setPrimeModulo(readFromHexadecimal("E95E4A5F737059DC60DFC7AD95B3D8139515620F"));
    el.setCoefA(readFromHexadecimal("340E7BE2A280EB74E2BE61BADA745D97E8F7C300"));
    el.setCoefB(readFromHexadecimal("1E589A8595423412134FAA2DBDEC95C8D8675E58"));
    return el;
}

ElipticCurveGeometry createElipticCurve2(){
    ElipticCurveGeometry el;
    el.setPrimeModulo(readFromHexadecimal("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFEE37"));
    el.setCoefA(0);
    el.setCoefB(3);
    return el;
}

ElipticCurveGeometry createElipticCurve3(){
    ElipticCurveGeometry el;
    el.setPrimeModulo(cpp_int("631"));
    el.setCoefA(30);
    el.setCoefB(34);
    return el;
}

void TestPointGeneration(){
    auto el = createElipticCurve();
    Point P = el.gen();
    ASSERT(el.onCurve(P));
}

void TestPointAddition(){
    auto el = createElipticCurve();
    Point P = {BigintOnField("846930889"), BigintOnField("914895438423939880721154620442986359470249180327")};
    Point Q = {BigintOnField("1681692781"), BigintOnField("413392049233386513684462954206000931924331146767")};
    Point P2 = {BigintOnField("656140758275595462266911470733314801368742930868"), BigintOnField("520016598759908788970158118245083307454255169535")};
    Point PQ = {BigintOnField("793516751077922864463424192516420234268439968247"), BigintOnField("369649928168204083999552433310618595221275787400")};
    ASSERT_EQUAL(el.sum(P, P), P2);
    ASSERT_EQUAL(el.sum(P, Q), PQ);
}

void TestPointMul(){
    auto el = createElipticCurve();
    Point P = {BigintOnField("846930889"), BigintOnField("914895438423939880721154620442986359470249180327")};
    Point P_ = {BigintOnField("744844647859786955424391562333164932941599456494"), BigintOnField("68976613961526457558914601512385766959755651551")};
    ASSERT_EQUAL(el.mul(P, cpp_int("1000")), P_);
}

void TestLine(){
    auto el = createElipticCurve2();
    Point P = el.gen();
    Point Q = el.gen();
    Line l1 = el.lineTroughPoints(P, Q);
    Line l2 = el.lineTroughPoints(P, P);
    Line l3 = el.verticalLine(P);
    
    ASSERT(l1.value(P) == 0);
    ASSERT(l1.value(Q) == 0);
    ASSERT(l2.value(P) == 0);
    ASSERT(l3.value(P) == 0);
}




void TestWeil(){
    LOG_DURATION("LOL_PRIKOL");
    Point basePoint = {cpp_int("0xDB4FF10EC057E9AE26B07D0280B7F4341DA5D1B1EAE06C7D"), cpp_int("0x9B2F2F6D9C5628A7844163D015BE86344082AA88D95E2F9D")};
    cpp_int order = cpp_int("0xFFFFFFFFFFFFFFFFFFFFFFFE26F2FC170F69466A74DEFD8D");
    auto el = createElipticCurve2();
    Point P = el.mul(basePoint, 2343432432243L);
    Point Q = el.mul(basePoint, 4233434343432443243L);
    auto first = weil_pairing(el, P, Q, order);
    auto second = weil_pairing(el, el.sum(P, P), Q, order);
    ASSERT(first * first == second);
}

void TestWeil2(){
    auto el = createElipticCurve3();
    Point P = {36, 60};
    Point Q = {121, 387};
    cpp_int order = 5;
    Point S = {0, 36};
    auto first = weil_pairing(el, P, Q, order);
    auto second = weil_pairing(el, el.sum(P, P), Q, order);
    ASSERT(first * first == second);
}



cpp_int BigintOnField::primeModulo = readFromHexadecimal("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFEE37");
int main() {
    {
        LOG_DURATION("Test time");
        TestRunner tr;
        RUN_TEST(tr, TestPointGeneration);
        RUN_TEST(tr, TestPointAddition);
        RUN_TEST(tr, TestPointMul);
        RUN_TEST(tr, TestLine);
        RUN_TEST(tr, TestWeil);
    }
    return 0;
}
