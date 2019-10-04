// ================================================================
//
// This software was developed by employees of the National Institute of
// Standards and Technology (NIST), an agency of the Federal Government.
// Pursuant to title 17 United States Code Section 105, works of NIST employees
// are not subject to copyright protection in the United States and are
// considered to be in the public domain. Permission to freely use, copy,
// modify, and distribute this software and its documentation without fee is
// hereby granted, provided that this notice and disclaimer of warranty appears
// in all copies.
//
// THE SOFTWARE IS PROVIDED 'AS IS' WITHOUT ANY WARRANTY OF ANY KIND, EITHER
// EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT LIMITED TO, ANY
// WARRANTY THAT THE SOFTWARE WILL CONFORM TO SPECIFICATIONS, ANY IMPLIED
// WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, AND FREEDOM
// FROM INFRINGEMENT, AND ANY WARRANTY THAT THE DOCUMENTATION WILL CONFORM TO
// THE SOFTWARE, OR ANY WARRANTY THAT THE SOFTWARE WILL BE ERROR FREE. IN NO
// EVENT SHALL NIST BE LIABLE FOR ANY DAMAGES, INCLUDING, BUT NOT LIMITED TO,
// DIRECT, INDIRECT, SPECIAL OR CONSEQUENTIAL DAMAGES, ARISING OUT OF,
// RESULTING FROM, OR IN ANY WAY CONNECTED WITH THIS SOFTWARE, WHETHER OR NOT
// BASED UPON WARRANTY, CONTRACT, TORT, OR OTHERWISE, WHETHER OR NOT INJURY WAS
// SUSTAINED BY PERSONS OR PROPERTY OR OTHERWISE, AND WHETHER OR NOT LOSS WAS
// SUSTAINED FROM, OR AROSE OUT OF THE RESULTS OF, OR USE OF, THE SOFTWARE OR
// SERVICES PROVIDED HEREUNDER.
//
// Distributions of NIST software should also include copyright and licensing
// statements of any third-party software that are legally bundled with the
// code in compliance with the conditions of those licenses.
// 
// ================================================================

// ================================================================
// 
// Authors: Derek Juba <derek.juba@nist.gov>
// Created: Fri Sep 11 09:43:30 2015 EDT
//
// ================================================================

#ifndef VOXELS_H
#define VOXELS_H

#ifdef USE_MPI
#include <mpi.h>
#endif

#include <iostream>
#include <string>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <zlib.h>

#include "Vector3.h"
#include "Cuboid.h"

// ================================================================

namespace zeno {
  
/// Represents an axis-aligned box of voxels.  Voxels are represented as a data
/// array and as cuboids.
///
template <class T>
class Voxels {
public:
  Voxels();
  ~Voxels();

  bool loadFitsGz(std::string const & inputFileName);

  void generateCuboids(std::vector<Cuboid<T> > * cuboids);

  bool isEmpty() const;

  unsigned short getData(int i, int j, int k) const;

  unsigned short const * getDataArray() const;

  Vector3<int> getNumVoxels() const;
  Vector3<double> getVoxelDims() const;

  void mpiSend() const;
  void mpiReceive();

private:
  unsigned short reverseEndianness(unsigned short itemToReverse);

  bool voxelsLoaded;

  Vector3<int> numVoxels;
  Vector3<double> voxelDims;

  unsigned short * data;
};

/// Construct an empty box of voxels.
///
template <class T>
Voxels<T>::Voxels() 
: voxelsLoaded(false),
  numVoxels(),
  voxelDims(),
  data(NULL) {

}

template <class T>
Voxels<T>::~Voxels() {
  delete [] data;
}

/// Load voxel data from the given .fits.gz file.
///
template <class T>
bool 
Voxels<T>::loadFitsGz(std::string const & inputFileName) {

  gzFile inputFile = gzopen(inputFileName.c_str(), "rb");

  if (inputFile == NULL) {
    std::cerr << "Error opening input fits file " << inputFileName << std::endl;
    return false;
  }

  int bitsPerPixel = 0;

  int numVoxelsArr[3];
  double voxelDimsArr[3];
  char unitsArr[3*81];

  int readError = 0;

  char fitsHeaderLine[81];
  readError = readError || (gzgets(inputFile, fitsHeaderLine, 81) == NULL);
  sscanf(fitsHeaderLine, "SIMPLE  = %*c");
  readError = readError || (gzgets(inputFile, fitsHeaderLine, 81) == NULL);
  sscanf(fitsHeaderLine, "BITPIX  = %d", &bitsPerPixel);
  readError = readError || (gzgets(inputFile, fitsHeaderLine, 81) == NULL);
  sscanf(fitsHeaderLine, "NAXIS   = %*d");
  readError = readError || (gzgets(inputFile, fitsHeaderLine, 81) == NULL);
  sscanf(fitsHeaderLine, "NAXIS1  = %d", numVoxelsArr + 0);
  readError = readError || (gzgets(inputFile, fitsHeaderLine, 81) == NULL);
  sscanf(fitsHeaderLine, "NAXIS2  = %d", numVoxelsArr + 1);
  readError = readError || (gzgets(inputFile, fitsHeaderLine, 81) == NULL);
  sscanf(fitsHeaderLine, "NAXIS3  = %d", numVoxelsArr + 2);
  readError = readError || (gzgets(inputFile, fitsHeaderLine, 81) == NULL);
  sscanf(fitsHeaderLine, "CDELT1  = %lf", voxelDimsArr + 0);
  readError = readError || (gzgets(inputFile, fitsHeaderLine, 81) == NULL);
  sscanf(fitsHeaderLine, "CDELT2  = %lf", voxelDimsArr + 1);
  readError = readError || (gzgets(inputFile, fitsHeaderLine, 81) == NULL);
  sscanf(fitsHeaderLine, "CDELT3  = %lf", voxelDimsArr + 2);
  readError = readError || (gzgets(inputFile, fitsHeaderLine, 81) == NULL);
  sscanf(fitsHeaderLine, "CTYPE1  = %s", unitsArr + 0*81);
  readError = readError || (gzgets(inputFile, fitsHeaderLine, 81) == NULL);
  sscanf(fitsHeaderLine, "CTYPE2  = %s", unitsArr + 1*81);
  readError = readError || (gzgets(inputFile, fitsHeaderLine, 81) == NULL);
  sscanf(fitsHeaderLine, "CTYPE3  = %s", unitsArr + 2*81);
  readError = readError || (gzgets(inputFile, fitsHeaderLine, 81) == NULL);
  sscanf(fitsHeaderLine, "END");

  numVoxels.setXYZ(numVoxelsArr[0],
		   numVoxelsArr[1],
		   numVoxelsArr[2]);

  voxelDims.setXYZ(voxelDimsArr[0],
		   voxelDimsArr[1],
		   voxelDimsArr[2]);

  for (int numHeaderLines = 13; numHeaderLines < 36; numHeaderLines++) {
    readError = readError || (gzgets(inputFile, fitsHeaderLine, 81) == NULL);
  }

  if (readError) {
    std::cerr << "Error reading input fits file " << inputFileName << std::endl;
    return false;
  }

  int numElements = numVoxelsArr[0] * numVoxelsArr[1] * numVoxelsArr[2];

  delete [] data;

  data = new unsigned short[numElements];

  if (data == NULL) {
    std::cerr << "Error allocating memory for output data" << std::endl;
    return false;
  }

  if (bitsPerPixel == 8) {
    unsigned char * inputData = new unsigned char[numElements];

    if (inputData == NULL) {
      std::cerr << "Error allocating memory for input data" << std::endl;
      return false;
    }

    int elementsRead = 
      gzread(inputFile, inputData, sizeof(unsigned char) * numElements);

    if (elementsRead != (int)sizeof(unsigned char) * numElements) {
      std::cerr << "Error reading input file " << inputFileName << std::endl;
      return false;
    }

    for (int i = 0; i < numElements; i++) {
      data[i] = inputData[i];
    }

    delete [] inputData;
  }
  else if (bitsPerPixel == 16) {
    unsigned short * inputData = new unsigned short[numElements];

    if (inputData == NULL) {
      std::cerr << "Error allocating memory for input data" << std::endl;
      return false;
    }

    int elementsRead = 
      gzread(inputFile, inputData, sizeof(unsigned short) * numElements);

    if (elementsRead != (int)sizeof(unsigned short) * numElements) {
      std::cerr << "Error reading input file " << inputFileName << std::endl;
      return false;
    }

    for (int i = 0; i < numElements; i++) {
      data[i] = reverseEndianness(inputData[i]);
    }

    delete [] inputData;
  }
  else {
    std::cerr << "Error: invalid number of bits per pixel " << bitsPerPixel
	      << std::endl;
    return false;
  }

  gzclose(inputFile);

  voxelsLoaded = true;

  return true;
}

template <class T>
unsigned short
Voxels<T>::reverseEndianness(unsigned short itemToReverse) {

  unsigned short reversedItem = 0;

  for (unsigned int byteIndex = 0; 
       byteIndex < sizeof(itemToReverse); 
       byteIndex++) {
    ((char *) &reversedItem)[byteIndex] = 
      ((char *) &itemToReverse)[sizeof(itemToReverse) - 1 - byteIndex];
  }

  return reversedItem;
}

/// Generate cuboids representing the voxels and push them into the provided
/// vector (existing cuboids in the vector are not affected).
///
template <class T>
void 
Voxels<T>::generateCuboids
(std::vector<Cuboid<T> > * cuboids) {
  for (int k = 0; k < numVoxels.getK(); k++) {
  for (int j = 0; j < numVoxels.getJ(); j++) {
  for (int i = 0; i < numVoxels.getI(); i++) {

    if (getData(i, j, k) != 0) {
      Vector3<T> minCoords(i,   j,   k);
      Vector3<T> maxCoords(i+1, j+1, k+1);

      minCoords *= voxelDims;
      maxCoords *= voxelDims;

      cuboids->emplace_back(minCoords, maxCoords);
    }
  }
  }
  }
}

template <class T>
bool
Voxels<T>::isEmpty()
const {
  return !voxelsLoaded;
}

/// Return the value of the voxel at the given coordinates.
///
template <class T>
unsigned short 
Voxels<T>::getData(int i, int j, int k) 
const {

  assert(i >= 0 && i < numVoxels.getI() &&
	 j >= 0 && j < numVoxels.getJ() &&
	 k >= 0 && k < numVoxels.getK());

  int index =
    k * numVoxels.getJ() * numVoxels.getI() +
    j * numVoxels.getI() +
    i;

  return data[index];
}

/// Return the values of the voxels as a raw, linear array.
///
template <class T>
unsigned short const * 
Voxels<T>::getDataArray() 
const {
  return data;
}

/// Return the number of voxels along each dimension.
///
template <class T>
Vector3<int> 
Voxels<T>::getNumVoxels() 
const {
  return numVoxels;
}

/// Return the dimensions of a voxel in length units.
///
template <class T>
Vector3<double> 
Voxels<T>::getVoxelDims() 
const {
  return voxelDims;
}

/// Broadcasts the voxels over MPI.
///
template <class T>
void 
Voxels<T>::mpiSend() const {
#ifdef USE_MPI
  int numDataElements = 0;

  if (voxelsLoaded) {
    numDataElements = 
      numVoxels.getX() * 
      numVoxels.getY() * 
      numVoxels.getZ();
  }
  else {
    numDataElements = -1;
  }

  MPI_Bcast(&numDataElements, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (voxelsLoaded) {
    MPI_Bcast(data, numDataElements, MPI_UNSIGNED_SHORT, 0, MPI_COMM_WORLD);

    int numVoxelsArray[3];
    double voxelDimsArray[3];

    for (int i = 0; i < 3; ++i) {
      numVoxelsArray[i] = numVoxels.get(i);
      voxelDimsArray[i] = voxelDims.get(i);
    }

    MPI_Bcast(numVoxelsArray, 3, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(voxelDimsArray, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }
#endif
}

/// Receives a set of voxels over MPI.  Any existing voxel data is replaced.
///
template <class T>
void 
Voxels<T>::mpiReceive() {
#ifdef USE_MPI
  int numDataElements = 0;

  MPI_Bcast(&numDataElements, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (numDataElements == -1) {
    delete [] data;

    data = NULL;

    voxelsLoaded = false;
  }
  else {
    delete [] data;

    data = new unsigned short[numDataElements];

    MPI_Bcast(data, numDataElements, MPI_UNSIGNED_SHORT, 0, MPI_COMM_WORLD);

    int numVoxelsArray[3];
    double voxelDimsArray[3];

    MPI_Bcast(numVoxelsArray, 3, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(voxelDimsArray, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    for (int i = 0; i < 3; ++i) {
      numVoxels.set(i, numVoxelsArray[i]);
      voxelDims.set(i, voxelDimsArray[i]);
    }

    voxelsLoaded = true;
  }
#endif
}

}

// ================================================================

#endif

