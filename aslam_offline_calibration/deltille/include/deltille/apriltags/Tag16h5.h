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

/** Tag family with 30 distinct codes.
    bits: 16,  minimum hamming: 5,  minimum complexity: 5

    Max bits corrected       False positive rate
            0                    0.045776 %
            1                    0.778198 %
            2                    6.271362 %

    Generation time: 0.309000 s

    Hamming distance between pairs of codes (accounting for rotation):

       0  0
       1  0
       2  0
       3  0
       4  0
       5  120
       6  172
       7  91
       8  33
       9  13
      10  6
      11  0
      12  0
      13  0
      14  0
      15  0
      16  0
**/

#pragma once

namespace AprilTags {

const unsigned long long t16h5[] = {
    0x231bLL, 0x2ea5LL, 0x346aLL, 0x45b9LL, 0x79a6LL, 0x7f6bLL,
    0xb358LL, 0xe745LL, 0xfe59LL, 0x156dLL, 0x380bLL, 0xf0abLL,
    0x0d84LL, 0x4736LL, 0x8c72LL, 0xaf10LL, 0x093cLL, 0x93b4LL,
    0xa503LL, 0x468fLL, 0xe137LL, 0x5795LL, 0xdf42LL, 0x1c1dLL,
    0xe9dcLL, 0x73adLL, 0xad5fLL, 0xd530LL, 0x07caLL, 0xaf2eLL};

static const TagCodes tagCodes16h5 =
    TagCodes("t16h5", 16, 5, t16h5, sizeof(t16h5) / sizeof(t16h5[0]));
}
