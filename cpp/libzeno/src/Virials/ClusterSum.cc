// ================================================================
//
/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
// 
// ================================================================

// ================================================================
// 
// Authors:
// Created:
//
// ================================================================

#include "Virials/ClusterSum.h"

/// Constructs the class to compute cluster sum.
///

using namespace zeno;

template <class T>
ClusterSum<T>::
ClusterSum(std::vector<Particle<T> *> * particles, OverlapTester<T> const * overlapTester):
particles(particles),
overlapTester(overlapTester){
}

template <class T>
ClusterSum<T>::
  ~ClusterSum() {
}

/// Constructs a sub class of ClusterSum to compute cluster sum for chains.
///

template <class T>
ClusterSumChain<T>::
ClusterSumChain(std::vector<Particle<T> *> * particles, double diameter,
                double ringFac, double chainFac):
        ClusterSum<T>(particles, NULL), diameter(diameter), ringFac(ringFac), chainFac(chainFac){
}

template <class T>
ClusterSumChain<T>::
~ClusterSumChain() {
}

/// Computes cluster sum for chains.
///
template <class T>
double
ClusterSumChain<T>::
value(){
    const int n = ClusterSum<T>::particles->size();
    const int nf1 = (1 << (n - 1));
    double fValues[n][n];
    for(int iMol1 = 0; iMol1 < n; ++iMol1){
        for(int iMol2 = iMol1 + 1; iMol2 < n; ++iMol2){
            Vector3<T> x = ClusterSum<T>::particles->at(iMol1)->getCenter();
            Vector3<T> y = ClusterSum<T>::particles->at(iMol2)->getCenter();
            Vector3<T> distCenterVec = x - y;
            T distCenterSqr = distCenterVec.getMagnitudeSqr();
            bool overlapped = distCenterSqr < diameter*diameter;
            fValues[iMol1][iMol2] = fValues[iMol2][iMol1] = overlapped ? 1 : 0;
        }
    }
    double nC[n-1][nf1];
    //nC[m][i] is the number of chains beginning at last vertex and ending at m, traversing all points in i
    //Start with all pairwise paths from last vertex to each other vertex
    for(int m = 0; m < (n - 1); ++m){
        nC[m][(1 << m)] = fValues[m][n - 1];
    }
    //All other paths
    for(int i = 1; i < nf1; ++i) {//i excludes the last bit, which is implicit in all partitions
        //the following two loops generate all pairs formed by each bit in i with each bit not in i
        //loop over bits not in i; start with full complement of i (i^(nf1-1)), and in each iteration
        //get lowest bit (im=(iC&-iC)) and strip it from complement (iC^=im) until complement is empty (iC=0)
        for(int iC = i ^ (nf1 - 1), im = (iC & -iC); iC > 0; iC ^= im, im = (iC & -iC)) {
            int m = log2(im);
            int iim = i|im;
            nC[m][iim] = 0;
            //loop over bits in i, in same manner as loop over complement
            for (int it = i, ik = (it & -it); ik > 0; it ^= ik, ik = (it & -it)) {
                int k = log2(ik);
                nC[m][iim] += fValues[m][k] * nC[k][i];
            }
        }//end for(iC)
    }//end for(i)
    double ringValue = 0;
    double chainValue = 0;
    if (ringFac != 0.0) {
        for (int m = 0; m < n - 1; m++) {
            ringValue += nC[m][nf1-1] * fValues[m][n-1];
        }
    }
    if (chainFac != 0.0) {
        //Sum chains in which last (n-1) vertex is not a leaf.
        //Consider all partitions, counting paths beginning in one partition and ending in its complement
        //Use same looping structure as employed above
        for (int iS = 1; iS < nf1; iS += 2) {//keep 1 in iS-partition to prevent double counting
            //loop over bits not in iS
            int iSComp = iS^(nf1-1);
            for (int iC = iSComp, im = (iC & -iC); iC > 0; iC ^= im, im = (iC & -iC)) {
                int m = log2(im);
                //loop over bits in iS
                for (int it = iS, ik = (it & -it); ik > 0; it ^= ik, ik = (it & -it)) {
                    int k = log2(ik);
                    chainValue += nC[m][iSComp] * nC[k][iS];
                }
            }
        }
        //Sum chains where last (n-1) vertex is a leaf
        for (int m = 0; m < n - 1; m++) {
            chainValue += nC[m][nf1-1];
        }
    }//end if(chainFrac)
    return chainFac*chainValue + ringFac*ringValue;
}

/// Constucts a sub class of ClusterSum to compute cluster sum using Wheatley Recursion.
///
template <class T>
ClusterSumWheatleyRecursion<T>::
ClusterSumWheatleyRecursion(std::vector<Particle<T> *> * particles, OverlapTester<T> const * overlapTester):
ClusterSum<T>(particles, overlapTester){
    const int n = ClusterSum<T>::particles->size();
    int factorial = 1;
    for (int m = 2; m <= n; ++m){
        factorial *= m;
    }
    preFac = -(n - 1.0)/factorial;
}

template <class T>
ClusterSumWheatleyRecursion<T>::
~ClusterSumWheatleyRecursion() {
}

/// Computes cluster sum using Wheatley Recursion.
///
template <class T>
double
ClusterSumWheatleyRecursion<T>::
value() {
    const int n = ClusterSum<T>::particles->size();
    const int nf = (1 << n);
    double fQ[nf], fC[nf];
    double fA[nf], fB[nf];
    for(int iMol1 = 0; iMol1 < n; ++iMol1){
        int i = 1 << iMol1;
        fQ[i] = 1.0;
        for(int iMol2 = iMol1 + 1; iMol2 < n; ++iMol2){
            bool overlapped = ClusterSum<T>::overlapTester->isOverlapped(ClusterSum<T>::particles->at(iMol1), ClusterSum<T>::particles->at(iMol2));
            fQ[i|(1<<iMol2)] = overlapped ? 0 : 1;
        }
    }
    //generate all partitions and compute
    for (int i = 3; i < nf; ++i){
        int j = i & -i; //lowest bit in i
        if (i == j) continue; //1-point set
        int k = i & ~j;
        if (k == (k & -k)) {
            // 2-point set
            continue;
        }
        fQ[i] = fQ[k];
        if (fQ[i] == 0) {
            continue;
        }
        for (int l = (j << 1); l < i; l = (l << 1)){
            if ( (l&i) == 0 ) continue;
            fQ[i] *= fQ[l|j];
        }
        if (fQ[i] == 0) {
            continue;
        }
    }
    //Compute the fC's
    for (int i = 1; i < nf; ++i){
        fC[i] = fQ[i];
        int iLowBit = i & -i;
        int inc = iLowBit << 1;
        for (int j = iLowBit; j < i; j += inc){
            int jComp = i & ~j;
            while ((j|jComp) != i && j<i){
                int jHighBits = j ^ iLowBit;
                int jlow = jHighBits & -jHighBits;
                j += jlow;
                jComp = (i & ~j);
            }
            if (j==i) break;
            fC[i] -= fC[j] * fQ[jComp];
        }
    }
    //find fA1
    for (int i = 2; i < nf; i += 2){
        //all even sets don't contain 1
        fB[i] = fC[i];
    }
    fA[1] = 0;
    fB[1] = fC[1];
    for (int i = 3; i < nf; i += 2){
        //every set will contain 1
        fA[i] = 0;
        fB[i] = fC[i];
        int ii = i - 1;//all bits in i but lowest
        int iLow2Bit = (ii & -ii);//next lowest bit
        int jBits = 1 | iLow2Bit;
        if (jBits == i) continue;
        int iii = ii ^ iLow2Bit; //i with 2 lowest bits off
        int jInc = (iii & -iii);//3rd lowest bit, alsso increment for j
        for (int j = jBits; j < i; j += jInc){//sum over partitions of i containing j Bits
            int jComp = (i & ~j);//subset of i complementing j
            while ((j|jComp) != i && j<i){
                int jHighBits = j ^ jBits;
                int jlow = jHighBits & -jHighBits;
                j += jlow;
                jComp = (i & ~j);
            }
            if (j == i) break;
            fA[i] += fB[j] * fC[jComp|1];
        }
        //remove from B graphs that contain articulation point 0.
        fB[i] -= fA[i];
    }
    for (int v = 1; v < n; ++v){
        int vs1 = 1 << v;
        for (int i = vs1 + 1; i < nf; ++i){
            fA[i] = 0;
            if ( (i & vs1) == 0 ) continue;
            int iLowBit = (i & -i);
            if ( iLowBit == i ) continue;
            int jBits;
            int ii = i ^ iLowBit;
            int iLow2Bit = (ii & -ii);
            if ( iLowBit!=vs1 && iLow2Bit!=vs1 ){
                //v is not in the lowest 2 bits
                jBits = iLowBit | vs1;
                //we can only increment by the 2nd lowest
                int jInc = iLow2Bit;
                for (int j = jBits; j < i; j += jInc){
                    if ( (j&jBits) != jBits ){
                        j |= vs1;
                        if (j==i) break;
                    }
                    int jComp = i & ~j;
                    while ((j|jComp) != i && j<i){
                        int jHighBits = j^jBits;
                        int jlow = jHighBits & -jHighBits;
                        j += jlow;
                        j |= vs1;
                        jComp = (i & ~j);
                    }
                    if (j==i) break;
                    fA[i] += fB[j] * (fB[jComp|vs1] + fA[jComp|vs1]);
                }
            }
            else{
                //lowest 2 bits contain v
                jBits = iLowBit | iLow2Bit;
                if (jBits == i) continue; // no bits left jComp
                int iii = ii ^ iLow2Bit;
                int jInc = ( iii & -iii);
                //at this point jBits has (lowest bit + v)
                for (int j = jBits; j < i; j += jInc){//sum over partitions of i
                    int jComp = i & ~j;
                    while ((j|jComp) != i && j<i){
                        int jHighBits = j^jBits;
                        int jlow = jHighBits & -jHighBits;
                        j += jlow;
                        jComp = (i & ~j);
                    }
                    if (j==i) break;
                    fA[i] += fB[j] * (fB[jComp|vs1] + fA[jComp|vs1]);
                }
            }
            fB[i] -= fA[i];
        }
    }
    return preFac*fB[nf-1];
}

