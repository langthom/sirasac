
#include "RandomNumber.h"

/*------------------------------------------------------------------------------------------------------------------*
 * Re-implementation of the Mersenne twister PRNG functionality from the Cython implementation of random.random().  *
 * We decided to do this since the C rand() function is known to produce biased results and since we made the       *
 * decision to create a C library and not to use C++, we do not have access to the random functions therein.        *
 * Please note that that code is literally taken from the Python implementation and was not checked by use.         *
 * Source:  https://github.com/python/cpython/blob/530f506ac91338b55cf2be71b1cdf50cb077512f/Modules/_randommodule.c *
 *------------------------------------------------------------------------------------------------------------------*/

/* ------------------------------------------------------------------
   The code in this module was based on a download from:
      http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MT2002/emt19937ar.html
   It was modified in 2002 by Raymond Hettinger as follows:
    * the principal computational lines untouched.
    * renamed genrand_res53() to random_random() and wrapped
      in python calling/return code.
    * genrand_int32() and the helper functions, init_genrand()
      and init_by_array(), were declared static, wrapped in
      Python calling/return code.  also, their global data
      references were replaced with structure references.
    * unused functions from the original were deleted.
      new, original C python code was added to implement the
      Random() interface.

   The following are the verbatim comments from the original code:
   A C-program for MT19937, with initialization improved 2002/1/26.
   Coded by Takuji Nishimura and Makoto Matsumoto.
   Before using, initialize the state by using init_genrand(seed)
   or init_by_array(init_key, key_length).
   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
    notice, this list of conditions and the following disclaimer.
     2. Redistributions in binary form must reproduce the above copyright
    notice, this list of conditions and the following disclaimer in the
    documentation and/or other materials provided with the distribution.
     3. The names of its contributors may not be used to endorse or promote
    products derived from this software without specific prior written
    permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

   Any feedback is very welcome.
   http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
   email: m-mat @ math.sci.hiroshima-u.ac.jp (remove space)
*/

#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfU
#define UPPER_MASK 0x80000000U
#define LOWER_MASK 0x7fffffffU

typedef unsigned int uint32_t;

typedef struct {
  int index;
  uint32_t state[N];
} RandomObject;


static RandomObject g_randomObject;

static uint32_t genrand_int32() {
  uint32_t y;
  static const uint32_t mag01[2] = { 0x0U, MATRIX_A };
  uint32_t* mt = g_randomObject.state;

  if (g_randomObject.index >= N) {
    int i;

    for (i = 0; i < N - M; ++i) {
      y = (mt[i] & UPPER_MASK) | (mt[i + 1] & LOWER_MASK);
      mt[i] = mt[i + M] ^ (y >> 1) ^ mag01[y & 0x1U];
    }

    for (; i < N - 1; ++i) {
      y = (mt[i] & UPPER_MASK) | (mt[i + 1] & LOWER_MASK);
      mt[i] = mt[i + (M-N)] ^ (y >> 1) ^ mag01[y & 0x1U];
    }

    y = (mt[N-1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
    mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1U];

    g_randomObject.index = 0;
  }

  y = mt[g_randomObject.index++];
  y ^= (y >> 11);
  y ^= (y <<  7) & 0x9d2c5680U;
  y ^= (y << 15) & 0xefc60000U;
  y ^= (y >> 18);
  return y;
}

/* ----------------------------------- The (internal) interface functions. --------------------------------- */

double rand01() {
  uint32_t a = genrand_int32() >> 5, b = genrand_int32() >> 6;
  return (a*67108864.0+b) * (1.0/9007199254740992.0);
}

void setSeed(unsigned int seed) {
  uint32_t* mt = g_randomObject.state;
  mt[0] = seed;
  int mti;

  for (mti = 1; mti < N; ++mti) {
    mt[mti] = (1812433253U * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti);
  }

  g_randomObject.index = mti;
}
