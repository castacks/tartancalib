/**
* Copyright (C) 2017-present, Facebook, Inc.
*
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
* This library is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
* Lesser General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public
* License along with this library; if not, write to the Free Software
* Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
#pragma once
#ifndef INCLUDE_CHECKERBOARD_DETECTOR_TAGFAMILY_H
#define INCLUDE_CHECKERBOARD_DETECTOR_TAGFAMILY_H

#include <string>
#include <vector>

namespace AprilTags {

class TagCodes {
public:
  int bits;
  int minHammingDistance;
  std::vector<unsigned long long> codes;
  const char *family;

public:
  // created vector for all entries of codesA
  TagCodes(const char *family, int bits, int minHammingDistance,
           const unsigned long long *codesA, int num)
      : bits(bits), minHammingDistance(minHammingDistance),
        codes(codesA, codesA + num), family(family) {}
};
};

namespace orp {

namespace calibration {

struct TagDetection {

  //! Constructor
  TagDetection()
      : good(false), obsCode(), code(), id(), hammingDistance(), rotation() {}

  //! Constructor for manually creating tags in a world map
  TagDetection(int id)
      : good(false), obsCode(), code(), id(id), hammingDistance(), rotation() {}

  //! which tag family is this tag from
  const char *tagFamily;

  //! Is the detection good enough?
  bool good;

  //! Observed code
  long long obsCode;

  //! Matched code
  long long code;

  //! What was the ID of the detected tag?
  int id;

  //! The hamming distance between the detected code and the true code
  int hammingDistance;

  //! How many 90 degree rotations were required to align the code (internal use
  //! only)
  int rotation;
};

//! Generic class for all tag encoding families
class TagFamily {
public:
  //! The codes array is not copied internally and so must not be modified
  //! externally.
  TagFamily(const AprilTags::TagCodes &tagCodes, const float blackBorder);

  void setErrorRecoveryBits(int b);

  void setErrorRecoveryFraction(float v);

  /* if the bits in w were arranged in a d*d grid and that grid was
   * rotated, what would the new bits in w be?
   * The bits are organized like this (for d = 3):
   *
   *  8 7 6       2 5 8      0 1 2
   *  5 4 3  ==>  1 4 7 ==>  3 4 5    (rotate90 applied twice)
   *  2 1 0       0 3 6      6 7 8
   */
  static unsigned long long rotate90(unsigned long long w, int d);

  // where to take bit i from in the rotated output
  static const int triangularRotationTable[8][49];
  /*
   const int rotate2d[4]  = {3, 1, 0,
   2};
   const int rotate3d[9]  = {8, 6, 5, 1, 0,
   7, 3, 2,
   4};
   const int rotate4d[16] = {15, 13, 12, 8, 7, 1, 0,
   14, 10, 9, 3, 2,
   11, 5, 4,
   6};
   const int rotate5d[25] = {24, 22, 21, 17, 16, 10, 9, 1, 0,
   23, 19, 18, 12, 11, 3, 2,
   20, 14, 13,  5, 4,
   15,  7,  6,
   8};
   const int rotate6d[36] = {35, 33, 32, 28, 27, 21, 20, 12, 11, 1, 0,
   34, 30, 29, 23, 22, 14, 13,  3, 2,
   31, 25, 24, 16, 15,  5,  4,
   26, 18, 17,  7,  6,
   19,  9,  8,
   10};
   const int rotate7d[49] = {48, 46, 45, 41, 40, 34, 33, 25, 24, 14, 13,  1,  0,
   47, 43, 42, 36, 35, 27, 26, 16, 15,  3,  2,
   44, 38, 37, 29, 28, 18, 17,  5,  4,
   39, 31, 30, 20, 19,  7,  6,
   32, 22, 21,  9,  8,
   23, 11, 10,
   12};
   */
  /* if the bits in w were arranged in a d*d triangle and that grid was
   * rotated, what would the new bits in w be?
   * The bits are organized like this (for d = 3):
   *
   * 8                          4              0
   *   6                      3                  1
   * 5   7        ==>       7   2       ==>    2   5
   *   1   3              6   1                  3   6
   * 0   2   4          8   5   0              4   7   8
   */
  //  static unsigned long long rotateTriangle(unsigned long long w, int d);
  static unsigned long long rotateTriangle(unsigned long long w, int d);

  //! Computes the hamming distance between two unsigned long longs.
  static int hammingDistance(unsigned long long a, unsigned long long b);

  //! How many bits are set in the unsigned long long?
  static unsigned char popCountReal(unsigned long long w);

  static int popCount(unsigned long long w);

  //! Given an observed tag with code 'rCode', try to recover the id.
  /*  The corresponding fields of TagDetection will be filled in. */
  void decode(TagDetection &det, unsigned long long rCode) const;

  //! Given an observed tag with code 'rCode', try to recover the id.
  /*  The corresponding fields of TagDetection will be filled in. */
  void decodeTritag(TagDetection &det, unsigned long long rCode) const;

  //! Prints the hamming distances of the tag codes.
  void printHammingDistances() const;

  //! Unique name of the family
  const char *name;

  int tagFamilyOffset;

  //! Numer of pixels wide of the inner black border.
  float blackBorder;

  //! Number of bits in the tag. Must be n^2.
  int bits;

  //! Dimension of tag. e.g. for 16 bits, dimension=4. Must be sqrt(bits).
  int dimension;

  //! Minimum hamming distance between any two codes.
  /*  Accounting for rotational ambiguity? The code can recover
   *  (minHammingDistance-1)/2 bit errors.
   */
  int minimumHammingDistance;

  /* The error recovery value determines our position on the ROC
   * curve. We will report codes that are within errorRecoveryBits
   * of a valid code. Small values mean greater rejection of bogus
   * tags (but false negatives). Large values mean aggressive
   * reporting of bad tags (but with a corresponding increase in
   * false positives).
   */
  int errorRecoveryBits;

  //! The array of the codes. The id for a code is its index.
  std::vector<unsigned long long> codes;

  static const int popCountTableShift = 12;
  static const unsigned int popCountTableSize = 1 << popCountTableShift;
  static unsigned char popCountTable[popCountTableSize];

  //! Initializes the static popCountTable
  static class TableInitializer {
  public:
    TableInitializer() {
      for (unsigned int i = 0; i < TagFamily::popCountTableSize; i++)
        TagFamily::popCountTable[i] = TagFamily::popCountReal(i);
    }
  } initializer;
};

} // namespace calibration
} // namespace orp

#endif
