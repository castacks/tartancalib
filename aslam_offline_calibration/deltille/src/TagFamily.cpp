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

#include <iostream>

#include <algorithm>
#include <cmath>
#include <deltille/TagFamily.h>
#include <limits.h>
/**

 // example of instantiation of tag family:

 #include "apriltags/Tag36h11.h"
 #include "apriltags/TagFamily.h"
 TagFamily *tag36h11 = new TagFamily(tagCodes36h11);

 // available tag families:

 #include "apriltags/Tag16h5.h"
 #include "apriltags/Tag16h5_other.h"
 #include "apriltags/Tag25h7.h"
 #include "apriltags/Tag25h9.h"
 #include "apriltags/Tag36h11.h"
 #include "apriltags/Tag36h11_other.h"
 #include "apriltags/Tag36h9.h"

 */

namespace orp {
namespace calibration {

TagFamily::TagFamily(const AprilTags::TagCodes &tagCodes,
                     const float blackBorder)
    : blackBorder(blackBorder), bits(tagCodes.bits),
      dimension((int)std::sqrt((float)bits)),
      minimumHammingDistance(tagCodes.minHammingDistance), errorRecoveryBits(1),
      codes() {
  if (bits != dimension * dimension)
    std::cerr << "Error: TagFamily constructor called with bits=" << bits
              << "; must be a square number!" << std::endl;
  codes = tagCodes.codes;
  name = tagCodes.family;
}

void TagFamily::setErrorRecoveryBits(int b) { errorRecoveryBits = b; }

void TagFamily::setErrorRecoveryFraction(float v) {
  errorRecoveryBits = (int)(((int)(minimumHammingDistance - 1) / 2) * v);
}

unsigned long long TagFamily::rotate90(unsigned long long w, int d) {
  unsigned long long wr = 0;
  const unsigned long long oneLongLong = 1;

  for (int r = d - 1; r >= 0; r--) {
    for (int c = 0; c < d; c++) {
      int b = r + d * c;
      wr = wr << 1;

      if ((w & (oneLongLong << b)) != 0)
        wr |= 1;
    }
  }
  return wr;
}

const int TagFamily::triangularRotationTable[8][49] = {
    /* 0 */ {0},
    /* 1 */ {0},
    /* 2 */ {3, 1, 0, 2},
    /* 3 */ {8, 6, 5, 1, 0, 7, 3, 2, 4},
    /* 4 */ {15, 13, 12, 8, 7, 1, 0, 14, 10, 9, 3, 2, 11, 5, 4, 6},
    /* 5 */ {24, 22, 21, 17, 16, 10, 9, 1, 0,  23, 19, 18, 12,
             11, 3,  2,  20, 14, 13, 5, 4, 15, 7,  6,  8},
    /* 6 */ {35, 33, 32, 28, 27, 21, 20, 12, 11, 1,  0,  34,
             30, 29, 23, 22, 14, 13, 3,  2,  31, 25, 24, 16,
             15, 5,  4,  26, 18, 17, 7,  6,  19, 9,  8,  10},
    /* 7 */ {48, 46, 45, 41, 40, 34, 33, 25, 24, 14, 13, 1,  0,  47, 43, 42, 36,
             35, 27, 26, 16, 15, 3,  2,  44, 38, 37, 29, 28, 18, 17, 5,  4,  39,
             31, 30, 20, 19, 7,  6,  32, 22, 21, 9,  8,  23, 11, 10, 12}};

unsigned long long TagFamily::rotateTriangle(unsigned long long w, int d) {
  unsigned long long out = 0;
  if (d < 2 || d > 7)
    // we do not support small and rotation of more than 7 x 7 triangles
    return w;
  const int *rotate = triangularRotationTable[d];
  d *= d;
  for (int i = 0; i < d; ++i)
    out |= ((uint64_t(1) << rotate[i]) & w) ? (uint64_t(1) << i) : 0;
  return out;
}

int TagFamily::hammingDistance(unsigned long long a, unsigned long long b) {
  return popCount(a ^ b);
}

unsigned char TagFamily::popCountReal(unsigned long long w) {
  unsigned char cnt = 0;
  while (w != 0) {
    w &= (w - 1);
    ++cnt;
  }
  return cnt;
}

int TagFamily::popCount(unsigned long long w) {
  int count = 0;
  while (w != 0) {
    count += popCountTable[(unsigned int)(w & (popCountTableSize - 1))];
    w >>= popCountTableShift;
  }
  return count;
}

void TagFamily::decode(TagDetection &det, unsigned long long rCode) const {
  int bestId = -1;
  int bestHamming = INT_MAX;
  int bestRotation = 0;
  unsigned long long bestCode = 0;

  unsigned long long rCodes[4];
  rCodes[0] = rCode;
  rCodes[1] = rotate90(rCodes[0], dimension);
  rCodes[2] = rotate90(rCodes[1], dimension);
  rCodes[3] = rotate90(rCodes[2], dimension);

  for (unsigned int id = 0; id < codes.size(); id++) {
    for (unsigned int rot = 0; rot < 4; rot++) {
      int thisHamming = hammingDistance(rCodes[rot], codes[id]);
      if (thisHamming < bestHamming) {
        bestHamming = thisHamming;
        bestRotation = rot;
        bestId = id;
        bestCode = codes[id];
      }
    }
  }
  det.id = bestId;
  det.hammingDistance = bestHamming;
  det.rotation = bestRotation;
  det.good = (det.hammingDistance <= errorRecoveryBits);
  det.obsCode = rCode;
  det.code = bestCode;
  det.tagFamily = name;
}

void TagFamily::decodeTritag(TagDetection &det,
                             unsigned long long rCode) const {
  int bestId = -1;
  int bestHamming = INT_MAX;
  int bestRotation = 0;
  unsigned long long bestCode = 0;

  unsigned long long rCodes[3];
  rCodes[0] = rCode;
  rCodes[1] = rotateTriangle(rCodes[0], dimension);
  rCodes[2] = rotateTriangle(rCodes[1], dimension);

  for (unsigned int id = 0; id < codes.size(); id++) {
    for (unsigned int rot = 0; rot < 3; rot++) {
      int thisHamming = hammingDistance(rCodes[rot], codes[id]);
      if (thisHamming < bestHamming) {
        bestHamming = thisHamming;
        bestRotation = rot;
        bestId = id;
        bestCode = codes[id];
      }
    }
  }
  det.id = bestId;
  det.tagFamily = name;
  det.hammingDistance = bestHamming;
  det.rotation = bestRotation;
  det.good = (det.hammingDistance <= errorRecoveryBits);
  det.obsCode = rCode;
  det.code = bestCode;
}

void TagFamily::printHammingDistances() const {
  std::vector<int> hammings(dimension * dimension + 1);
  for (unsigned i = 0; i < codes.size(); i++) {
    unsigned long long r0 = codes[i];
    unsigned long long r1 = rotate90(r0, dimension);
    unsigned long long r2 = rotate90(r1, dimension);
    unsigned long long r3 = rotate90(r2, dimension);
    for (unsigned int j = i + 1; j < codes.size(); j++) {
      int d = std::min(std::min(hammingDistance(r0, codes[j]),
                                hammingDistance(r1, codes[j])),
                       std::min(hammingDistance(r2, codes[j]),
                                hammingDistance(r3, codes[j])));
      hammings[d]++;
    }
  }

  for (unsigned int i = 0; i < hammings.size(); i++)
    printf("hammings: %u = %d\n", i, hammings[i]);
}

unsigned char TagFamily::popCountTable[TagFamily::popCountTableSize];

TagFamily::TableInitializer TagFamily::initializer;

} // namespace
}
