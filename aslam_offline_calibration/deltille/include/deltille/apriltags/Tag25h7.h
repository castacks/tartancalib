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

/** Tag family with 242 distinct codes.
 bits: 25,  minimum hamming: 7,  minimum complexity: 8

 Max bits corrected       False positive rate
 0                    0.000721 %
 1                    0.018752 %
 2                    0.235116 %
 3                    1.893914 %

 Generation time: 72.585000 s

 Hamming distance between pairs of codes (accounting for rotation):

 0  0
 1  0
 2  0
 3  0
 4  0
 5  0
 6  0
 7  2076
 8  4161
 9  5299
 10  6342
 11  5526
 12  3503
 13  1622
 14  499
 15  114
 16  16
 17  3
 18  0
 19  0
 20  0
 21  0
 22  0
 23  0
 24  0
 25  0
 **/

#pragma once

namespace deltille {
namespace AprilTags {

const unsigned long long t25h7[] = {
    0x4b770dLL,  0x11693e6LL, 0x1a599abLL, 0xc3a535LL,  0x152aafaLL,
    0xaccd98LL,  0x1cad922LL, 0x2c2fadLL,  0xbb3572LL,  0x14a3b37LL,
    0x186524bLL, 0xc99d4cLL,  0x23bfeaLL,  0x141cb74LL, 0x1d0d139LL,
    0x1670aebLL, 0x851675LL,  0x150334eLL, 0x6e3ed8LL,  0xfd449dLL,
    0xaa55ecLL,  0x1c86176LL, 0x15e9b28LL, 0x7ca6b2LL,  0x147c38bLL,
    0x1d6c950LL, 0x8b0e8cLL,  0x11a1451LL, 0x1562b65LL, 0x13f53c8LL,
    0xd58d7aLL,  0x829ec9LL,  0xfaccf1LL,  0x136e405LL, 0x7a2f06LL,
    0x10934cbLL, 0x16a8b56LL, 0x1a6a26aLL, 0xf85545LL,  0x195c2e4LL,
    0x24c8a9LL,  0x12bfc96LL, 0x16813aaLL, 0x1a42abeLL, 0x1573424LL,
    0x1044573LL, 0xb156c2LL,  0x5e6811LL,  0x1659bfeLL, 0x1d55a63LL,
    0x5bf065LL,  0xe28667LL,  0x1e9ba54LL, 0x17d7c5aLL, 0x1f5aa82LL,
    0x1a2bbd1LL, 0x1ae9f9LL,  0x1259e51LL, 0x134062bLL, 0xe1177aLL,
    0xed07a8LL,  0x162be24LL, 0x59128bLL,  0x1663e8fLL, 0x1a83cbLL,
    0x45bb59LL,  0x189065aLL, 0x4bb370LL,  0x16fb711LL, 0x122c077LL,
    0xeca17aLL,  0xdbc1f4LL,  0x88d343LL,  0x58ac5dLL,  0xba02e8LL,
    0x1a1d9dLL,  0x1c72eecLL, 0x924bc5LL,  0xdccab3LL,  0x886d15LL,
    0x178c965LL, 0x5bc69aLL,  0x1716261LL, 0x174e2ccLL, 0x1ed10f4LL,
    0x156aa8LL,  0x3e2a8aLL,  0x2752edLL,  0x153c651LL, 0x1741670LL,
    0x765b05LL,  0x119c0bbLL, 0x172a783LL, 0x4faca1LL,  0xf31257LL,
    0x12441fcLL, 0x0d3748LL,  0xc21f15LL,  0xac5037LL,  0x180e592LL,
    0x7d3210LL,  0xa27187LL,  0x2beeafLL,  0x26ff57LL,  0x690e82LL,
    0x77765cLL,  0x1a9e1d7LL, 0x140be1aLL, 0x1aa1e3aLL, 0x1944f5cLL,
    0x19b5032LL, 0x169897LL,  0x1068eb9LL, 0xf30dbcLL,  0x106a151LL,
    0x1d53e95LL, 0x1348ceeLL, 0xcf4fcaLL,  0x1728bb5LL, 0xdc1eecLL,
    0x69e8dbLL,  0x16e1523LL, 0x105fa25LL, 0x18abb0cLL, 0xc4275dLL,
    0x6d8e76LL,  0xe8d6dbLL,  0xe16fd7LL,  0x1ac2682LL, 0x77435bLL,
    0xa359ddLL,  0x3a9c4eLL,  0x123919aLL, 0x1e25817LL, 0x02a836LL,
    0x1545a4LL,  0x1209c8dLL, 0xbb5f69LL,  0x1dc1f02LL, 0x5d5f7eLL,
    0x12d0581LL, 0x13786c2LL, 0xe15409LL,  0x1aa3599LL, 0x139aad8LL,
    0xb09d2aLL,  0x54488fLL,  0x13c351cLL, 0x976079LL,  0xb25b12LL,
    0x1addb34LL, 0x1cb23aeLL, 0x1175738LL, 0x1303bb8LL, 0xd47716LL,
    0x188ceeaLL, 0xbaf967LL,  0x1226d39LL, 0x135e99bLL, 0x34adc5LL,
    0x2e384dLL,  0x90d3faLL,  0x232713LL,  0x17d49b1LL, 0xaa84d6LL,
    0xc2ddf8LL,  0x1665646LL, 0x4f345fLL,  0x2276b1LL,  0x1255dd7LL,
    0x16f4cccLL, 0x4aaffcLL,  0xc46da6LL,  0x85c7b3LL,  0x1311fcbLL,
    0x9c6c4fLL,  0x187d947LL, 0x8578e4LL,  0xe2bf0bLL,  0xa01b4cLL,
    0xa1493bLL,  0x7ad766LL,  0xccfe82LL,  0x1981b5bLL, 0x1cacc85LL,
    0x562cdbLL,  0x15b0e78LL, 0x8f66c5LL,  0x3332bfLL,  0x12ce754LL,
    0x096a76LL,  0x1d5e3baLL, 0x27ea41LL,  0x14412dfLL, 0x67b9b4LL,
    0xdaa51aLL,  0x1dcb17LL,  0x4d4afdLL,  0x6335d5LL,  0xee2334LL,
    0x17d4e55LL, 0x1b8b0f0LL, 0x14999e3LL, 0x1513dfaLL, 0x765cf2LL,
    0x56af90LL,  0x12e16acLL, 0x1d3d86cLL, 0xff279bLL,  0x18822ddLL,
    0x99d478LL,  0x8dc0d2LL,  0x34b666LL,  0xcf9526LL,  0x186443dLL,
    0x7a8e29LL,  0x19c6aa5LL, 0x1f2a27dLL, 0x12b2136LL, 0xd0cd0dLL,
    0x12cb320LL, 0x17ddb0bLL, 0x05353bLL,  0x15b2cafLL, 0x1e5a507LL,
    0x120f1e5LL, 0x114605aLL, 0x14efe4cLL, 0x568134LL,  0x11b9f92LL,
    0x174d2a7LL, 0x692b1dLL,  0x39e4feLL,  0xaaff3dLL,  0x96224cLL,
    0x13c9f77LL, 0x110ee8fLL, 0xf17beaLL,  0x99fb5dLL,  0x337141LL,
    0x02b54dLL,  0x1233a70LL};

static const TagCodes tagCodes25h7 =
    TagCodes("t25h7", 25, 7, t25h7, sizeof(t25h7) / sizeof(t25h7[0]));
}
}