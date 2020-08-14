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
// Created: 2019-05-30
//
// ================================================================

#ifndef PARAMETERS_LOCAL_H_
#define PARAMETERS_LOCAL_H_

// ================================================================

#include <string>

// ================================================================

/// Collects the parameters that are used by the zeno command-line program.
///
/// For some parameters, tracks whether they 
/// have been manually set or are still at their default value.
///
class ParametersLocal
{
public:
  ParametersLocal();
  ~ParametersLocal();

  void setInputFileName(std::string const & inputFileName);
  std::string getInputFileName() const;

  void setXyzInputFileName(std::string const & xyzInputFileName);
  std::string getXyzInputFileName() const;
  bool getXyzInputFileNameWasSet() const;
  
  void setMapInputFileName(std::string const & mapInputFileName);
  std::string getMapInputFileName() const;
  bool getMapInputFileNameWasSet() const;
  
  void setCsvOutputFileName(std::string const & csvOutputFileName);
  std::string getCsvOutputFileName() const;
  bool getCsvOutputFileNameWasSet() const;
  
  int getMpiSize() const;
  void setMpiSize(int mpiSize);

  int getMpiRank() const;
  void setMpiRank(int mpiRank);

  void setSurfacePointsFileName(std::string const & surfacePointsFileName);
  std::string getSurfacePointsFileName() const;

  void setInteriorPointsFileName(std::string const & interiorPointsFileName);
  std::string getInteriorPointsFileName() const;

  void setPrintCounts(bool printCounts);
  bool getPrintCounts() const;

  void setPrintBenchmarks(bool printBenchmarks);
  bool getPrintBenchmarks() const;

  void mpiBroadcast(int root);

private:
  void serializeMpiBroadcast(int root) const;
  void mpiBroadcastDeserialize(int root);
  
  std::string inputFileName;

  std::string xyzInputFileName;
  bool xyzInputFileNameWasSet;

  std::string mapInputFileName;
  bool mapInputFileNameWasSet;
  
  std::string csvOutputFileName;
  bool csvOutputFileNameWasSet;

  int mpiSize;
  int mpiRank;

  std::string surfacePointsFileName;
  std::string interiorPointsFileName;

  bool printCounts;
  bool printBenchmarks;
};

// ================================================================

#endif  // #ifndef PARAMETERS_LOCAL_H_

