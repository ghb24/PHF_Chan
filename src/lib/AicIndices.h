/* Copyright (c) 2012  Gerald Knizia
 * 
 * This file is part of the bfint program
 * (See http://www.theochem.uni-stuttgart.de/~knizia)
 * 
 * bfint is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3.
 * 
 * bfint is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with bfint (LICENSE). If not, see http://www.gnu.org/licenses/
 */

/* AicIndices.h v20121011 EST [charge, Gerald Knizia] */
// This is: lc -> components required to form lower lc shells in OsrrB.
// (i.e., it encodes the OsrrBX recursion path)
typedef uint32_t
   FCompFlag;
uint OsrrB_iTabRequiredCompC[7] = {
   0, 1, 3, 6, 10, 15, 21
};
FCompFlag OsrrB_TabRequiredCompC[28] = {
   0x0, 0x1, 0x7, 0x1, 0x7, 0x3f, 0x1, 0x7, 0x27, 0x3ff, 0x1, 0x7, 0x7, 0x173, 0x7fff, 0x1,
   0x7, 0x7, 0x133, 0xbbb, 0x1fffff, 0x1, 0x7, 0x7, 0x131, 0x3b9, 0xf07e5, 0xfffffff
};
// This is: [iCartFull(Comp)] -> (Preferred Reduction Direction i, iCartRel(Comp-1i)).
// -1 is emitted when the reduced component does not exist (it is good if it doesn't).
char iCartRD[84][2] = {
   {0,-1}, {0,0}, {1,0}, {2,0}, {0,0}, {1,1},
   {2,2}, {0,1}, {0,2}, {1,2}, {0,0}, {0,1},
   {0,2}, {1,0}, {1,1}, {1,2}, {2,0}, {2,1},
   {2,2}, {0,5}, {0,0}, {0,1}, {2,6}, {1,4},
   {1,5}, {2,8}, {1,0}, {0,4}, {0,5}, {2,0},
   {2,1}, {0,8}, {1,6}, {2,4}, {1,8}, {0,0},
   {0,1}, {2,9}, {0,3}, {0,4}, {0,5}, {1,0},
   {0,7}, {0,8}, {1,3}, {1,4}, {1,5}, {0,9},
   {2,1}, {0,11}, {2,3}, {2,4}, {2,5}, {1,9},
   {2,7}, {2,8}, {0,0}, {1,6}, {0,2}, {1,7},
   {1,8}, {0,5}, {1,9}, {1,10}, {2,16}, {2,17},
   {0,6}, {0,7}, {0,8}, {0,9}, {0,10}, {1,5},
   {2,0}, {1,18}, {2,2}, {1,19}, {0,16}, {0,17},
   {0,18}, {0,19}, {2,8}, {2,9}, {1,16}, {1,17}

};
// This is: [iCartFull(Comp)] -> (nx,ny,nz) for monomial x^nx y^ny z^nz.
unsigned char iCartPow[84][3] = {
   {0,0,0}, {1,0,0}, {0,1,0}, {0,0,1}, {2,0,0}, {0,2,0},
   {0,0,2}, {1,1,0}, {1,0,1}, {0,1,1}, {3,0,0}, {1,2,0},
   {1,0,2}, {2,1,0}, {0,3,0}, {0,1,2}, {2,0,1}, {0,2,1},
   {0,0,3}, {1,1,1}, {4,0,0}, {2,2,0}, {2,0,2}, {0,4,0},
   {0,2,2}, {0,0,4}, {3,1,0}, {1,3,0}, {1,1,2}, {3,0,1},
   {1,2,1}, {1,0,3}, {2,1,1}, {0,3,1}, {0,1,3}, {5,0,0},
   {3,2,0}, {3,0,2}, {1,4,0}, {1,2,2}, {1,0,4}, {4,1,0},
   {2,3,0}, {2,1,2}, {0,5,0}, {0,3,2}, {0,1,4}, {4,0,1},
   {2,2,1}, {2,0,3}, {0,4,1}, {0,2,3}, {0,0,5}, {3,1,1},
   {1,3,1}, {1,1,3}, {6,0,0}, {4,2,0}, {4,0,2}, {2,4,0},
   {2,2,2}, {2,0,4}, {0,6,0}, {0,4,2}, {0,2,4}, {0,0,6},
   {5,1,0}, {3,3,0}, {3,1,2}, {1,5,0}, {1,3,2}, {1,1,4},
   {5,0,1}, {3,2,1}, {3,0,3}, {1,4,1}, {1,2,3}, {1,0,5},
   {4,1,1}, {2,3,1}, {2,1,3}, {0,5,1}, {0,3,3}, {0,1,5}

};
// This is: [iCartFull(Comp)] -> (iCartFull(Comp-1_x),iCartFull(Comp-1_y),iCartFull(Comp-1_z))
short iCartM1[84][3] = {
   {-1,-1,-1}, {0,-1,-1}, {-1,0,-1}, {-1,-1,0}, {1,-1,-1}, {-1,2,-1},
   {-1,-1,3}, {2,1,-1}, {3,-1,1}, {-1,3,2}, {4,-1,-1}, {5,7,-1},
   {6,-1,8}, {7,4,-1}, {-1,5,-1}, {-1,6,9}, {8,-1,4}, {-1,9,5},
   {-1,-1,6}, {9,8,7}, {10,-1,-1}, {11,13,-1}, {12,-1,16}, {-1,14,-1},
   {-1,15,17}, {-1,-1,18}, {13,10,-1}, {14,11,-1}, {15,12,19}, {16,-1,10},
   {17,19,11}, {18,-1,12}, {19,16,13}, {-1,17,14}, {-1,18,15}, {20,-1,-1},
   {21,26,-1}, {22,-1,29}, {23,27,-1}, {24,28,30}, {25,-1,31}, {26,20,-1},
   {27,21,-1}, {28,22,32}, {-1,23,-1}, {-1,24,33}, {-1,25,34}, {29,-1,20},
   {30,32,21}, {31,-1,22}, {-1,33,23}, {-1,34,24}, {-1,-1,25}, {32,29,26},
   {33,30,27}, {34,31,28}, {35,-1,-1}, {36,41,-1}, {37,-1,47}, {38,42,-1},
   {39,43,48}, {40,-1,49}, {-1,44,-1}, {-1,45,50}, {-1,46,51}, {-1,-1,52},
   {41,35,-1}, {42,36,-1}, {43,37,53}, {44,38,-1}, {45,39,54}, {46,40,55},
   {47,-1,35}, {48,53,36}, {49,-1,37}, {50,54,38}, {51,55,39}, {52,-1,40},
   {53,47,41}, {54,48,42}, {55,49,43}, {-1,50,44}, {-1,51,45}, {-1,52,46}

};
