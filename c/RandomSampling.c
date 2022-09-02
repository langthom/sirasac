
/**********************************************************************************
 * MIT License                                                                    *
 *                                                                                *
 * Copyright (c) 2022 Dr. Thomas Lang                                             *
 *                                                                                *
 * Permission is hereby granted, free of charge, to any person obtaining a copy   *
 * of this software and associated documentation files (the "Software"), to deal  *
 * in the Software without restriction, including without limitation the rights   *
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell      *
 * copies of the Software, and to permit persons to whom the Software is          *
 * furnished to do so, subject to the following conditions:                       *
 *                                                                                *
 * The above copyright notice and this permission notice shall be included in all *
 * copies or substantial portions of the Software.                                *
 *                                                                                *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR     *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,       *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE    *
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER         *
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,  *
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE  *
 * SOFTWARE.                                                                      *
 **********************************************************************************/

#include "RandomSampling.h"
#include "RandomNumber.h"

#include <assert.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdio.h>

/* Size of the chunks read and processed at once. This is rather small, which seems good for smaller and larger files on average sytems. */
#define CHUNK_SIZE 64ull << 20

/* Byte sizes for an element of the specified type. */
size_t typeSizes[] = {
  [TYPE_ID_WRAPPER(INT8)]  = sizeof(char),          [TYPE_ID_WRAPPER(INT16)]  = sizeof(short),          [TYPE_ID_WRAPPER(INT32)]  = sizeof(int),          [TYPE_ID_WRAPPER(INT64)]  = sizeof(long),
  [TYPE_ID_WRAPPER(UINT8)] = sizeof(unsigned char), [TYPE_ID_WRAPPER(UINT16)] = sizeof(unsigned short), [TYPE_ID_WRAPPER(UINT32)] = sizeof(unsigned int), [TYPE_ID_WRAPPER(UINT64)] = sizeof(unsigned long),
  [TYPE_ID_WRAPPER(FLOAT32)] = sizeof(float), [TYPE_ID_WRAPPER(FLOAT64)] = sizeof(double),
};

/* The span of representable values of each type. To avoid numerical issues, we employed
 * that for signed types the domain is [ -2^(p-1) .. 2^(p-1)-1 ], i.e., the data type span
 * is 2^(p-1)-1+2^(p-1) == 2^p - 1, which is the span of the unsigned types.                */
double __type_spans[] = {
  [TYPE_ID_WRAPPER(INT8)]  = UCHAR_MAX, [TYPE_ID_WRAPPER(INT16)]  = USHRT_MAX, [TYPE_ID_WRAPPER(INT32)]  = UINT_MAX, [TYPE_ID_WRAPPER(INT64)]  = ULONG_MAX,
  [TYPE_ID_WRAPPER(UINT8)] = UCHAR_MAX, [TYPE_ID_WRAPPER(UINT16)] = USHRT_MAX, [TYPE_ID_WRAPPER(UINT32)] = UINT_MAX, [TYPE_ID_WRAPPER(UINT64)] = ULONG_MAX,
  [TYPE_ID_WRAPPER(FLOAT32)] = FLT_MAX, [TYPE_ID_WRAPPER(FLOAT64)] = DBL_MAX,
};

/* The minimum representable values for each supported type. */
double __type_mins[] = {
  [TYPE_ID_WRAPPER(INT8)] = SCHAR_MIN, [TYPE_ID_WRAPPER(INT16)] = SHRT_MIN, [TYPE_ID_WRAPPER(INT32)] = INT_MIN, [TYPE_ID_WRAPPER(INT64)] = LONG_MIN,
  [TYPE_ID_WRAPPER(UINT8)] = 0, [TYPE_ID_WRAPPER(UINT16)] = 0, [TYPE_ID_WRAPPER(UINT32)] = 0, [TYPE_ID_WRAPPER(INT64)] = 0,
  [TYPE_ID_WRAPPER(FLOAT32)] = FLT_MIN, [TYPE_ID_WRAPPER(FLOAT64)] = DBL_MIN,
};

/* --------------------------------------------------------------------- Algorithm L ----------------------------------------------------------------------- */

/* The regular Algorithm L. */
enum RandomSamplingError AlgorithmL(void* sample, IOConfiguration const* ioConfiguration, SamplingConfiguration const* configuration)
{
  size_t sampleSize = configuration->sampleSize;
  size_t elementSize = typeSizes[ioConfiguration->typeTag];

  FILE* blobFile;
  if(ioConfiguration->inputIsFile) {
    blobFile = fopen(ioConfiguration->filename, "rb");
    if(!blobFile) return SAMPLING__SOURCE_FILE_CANNOT_BE_OPENED;

    if(fseek(blobFile, (long)ioConfiguration->headerOffset, SEEK_SET) != 0) {
      fclose(blobFile);
      return SAMPLING__COULD_NOT_JUMP_PAST_HEADER;
    }
  }

  // First, fill up the reservoir.
  if(ioConfiguration->inputIsFile) {
    size_t elementsRead = fread(sample, elementSize, sampleSize, blobFile);
    if(elementsRead < sampleSize) {
      fclose(blobFile);
      return SAMPLING__COULD_NOT_FILL_RESERVOIR;
    }
  } else {
    memcpy(sample, ioConfiguration->blob, sampleSize * elementSize);
  }

  double W = exp(log(rand01())/sampleSize);
  size_t i = sampleSize - 1;
  size_t numElements = ioConfiguration->numElements;
  char* samplePtr = (char*)sample;

  /* Now, perform the sampling with skipping as defined in Algorithm L. */
  while(i < numElements) {
    i += (size_t)(floor(log(rand01())/log(1.0-W))) + 1;

    if(i < numElements) {
      size_t r = (size_t)( rand01() * (sampleSize - 1) );

      if(ioConfiguration->inputIsFile) {
        if(fseek(blobFile, (long)(i * elementSize), SEEK_SET) != 0) {
          fclose(blobFile);
          return SAMPLING__COULD_NOT_JUMP_TO_NEXT_ELEMENT;
        }
        if(fread(samplePtr + r * elementSize, elementSize, 1, blobFile) != 1) {
          fclose(blobFile);
          return SAMPLING__COULD_NOT_READ_ELEMENT;
        }
      } else {
        memcpy(samplePtr + r * elementSize, (char*)ioConfiguration->blob + i * elementSize, elementSize);
      }

      W *= exp(log(rand01())/sampleSize);
    }
  }

  return SAMPLING__NO_ERROR;
}

/* ----------------------------------------------------------- Stratified sampling utilities --------------------------------------------------------------- */
inline void fillStratification__linear(size_t* stratification, int numberOfStrati, enum TypeTag typeTag) {
  for(int i = 0; i < numberOfStrati; ++i) {
    double stratVal = i * (1. - 1./__type_spans[typeTag])/((double)numberOfStrati - 1) + 1./__type_spans[typeTag];
    stratification[i] = (size_t)(stratVal * __type_spans[typeTag] + __type_mins[typeTag]);
  }
}

inline void fillStratification__exponential(size_t* stratification, int numberOfStrati, enum TypeTag typeTag) {
  double cumsum = exp2(1-numberOfStrati);
  stratification[0] = (size_t)(cumsum * __type_spans[typeTag] + __type_mins[typeTag]);

  for(int i = 1; i < numberOfStrati; ++i) {
    cumsum += exp2(i-numberOfStrati);
    stratification[i] = (size_t)(cumsum * __type_spans[typeTag] + __type_mins[typeTag]);
  }
}

inline void fillStratification__mixed(size_t* stratification, int numberOfStrati, enum TypeTag typeTag, size_t threshold) {
  int i = 0, N2 = numberOfStrati / 2;
  double powStep  = (log10((double)threshold)-1)/N2;
  double thresh01 = ((double)threshold)/__type_spans[typeTag];
  double linStep  = (1.0 - thresh01)/((double)numberOfStrati-N2-1);

  for(; i < N2; ++i) {
    stratification[i] = (size_t)(pow(10, 1 + i*powStep));
  }
  for(int j = 0; j < numberOfStrati - N2; ++j) {
    double stratVal = thresh01 + j * linStep;
    stratification[i++] = (size_t)(stratVal * __type_spans[typeTag] + __type_mins[typeTag]);
  }
}

size_t* get_stratification(SamplingConfiguration const* configuration, enum TypeTag typeTag) {
  int numberOfStrati = configuration->numberOfStrati;
  size_t* stratification = (size_t*)calloc(numberOfStrati, sizeof(size_t));
  assert(stratification);

  switch(configuration->stratificationType) {
  case LINEAR:      fillStratification__linear(     stratification, numberOfStrati, typeTag); break;
  case EXPONENTIAL: fillStratification__exponential(stratification, numberOfStrati, typeTag); break;
  case MIXED:       fillStratification__mixed(      stratification, numberOfStrati, typeTag, configuration->stratification.mixedThreshold); break;
  case PRECOMPUTED: 
    for(int i = 0; i < numberOfStrati; ++i) {
      stratification[i] = (size_t)(configuration->stratification.strati[i] * __type_spans[typeTag] + __type_mins[typeTag]);
    }
    break;
  default: /* do nothing */;
  }

  return stratification;
}

/* ------------------------------------------------------------------------------ */

/* Here we create one function that does the counting for each voxel type. This enables us to perform the switch
   over the voxel types before looping, which should be a little bit more performant.                            */

#define STR(x) #x
#define STRINGIFY(x) STR(x)
#define CONCATENATE(X) X

#define SAMPLE_SIZES_FUNCTION(typenam,type)                                                                                          \
inline void getSampleSizes_ ## typenam(size_t* sampleSizes, void const* _chunk, size_t chunkSize, size_t* stratification, int numStrati) {  \
  type* const chunk = (type* const) _chunk;                                                                                          \
                                                                                                                                     \
  _Pragma(STRINGIFY(CONCATENATE(omp parallel if(chunkSize >= 64))))                                                                  \
  {                                                                                                                                  \
    size_t* threadLocalCounts = (size_t*)calloc(numStrati, sizeof(size_t));                                                          \
    assert(threadLocalCounts);                                                                                                       \
    int i, stratumIndex;                                                                                                             \
                                                                                                                                     \
    _Pragma(STRINGIFY(CONCATENATE(omp for)))                                                                                         \
    for (i = 0; i < chunkSize; ++i) {                                                                                                \
      double value = (double)chunk[i];                                                                                               \
                                                                                                                                     \
      for (size_t stratumIdx = 0; stratumIdx < numStrati; ++stratumIdx) {                                                            \
        if (value <= stratification[stratumIdx]) {                                                                                   \
          threadLocalCounts[stratumIdx]++;                                                                                           \
          break;                                                                                                                     \
        }                                                                                                                            \
      }                                                                                                                              \
    }                                                                                                                                \
                                                                                                                                     \
    _Pragma(STRINGIFY(CONCATENATE(omp critical)))                                                                                    \
    for (stratumIndex = 0; stratumIndex < numStrati; ++stratumIndex) {                                                               \
      sampleSizes[stratumIndex] += threadLocalCounts[stratumIndex];                                                                  \
    }                                                                                                                                \
                                                                                                                                     \
    free(threadLocalCounts);                                                                                                         \
  }                                                                                                                                  \
}

SAMPLE_SIZES_FUNCTION(int8,    char)
SAMPLE_SIZES_FUNCTION(int16,   short)
SAMPLE_SIZES_FUNCTION(int32,   int)
SAMPLE_SIZES_FUNCTION(int64,   long)
SAMPLE_SIZES_FUNCTION(uint8,   unsigned char)
SAMPLE_SIZES_FUNCTION(uint16,  unsigned short)
SAMPLE_SIZES_FUNCTION(uint32,  unsigned int)
SAMPLE_SIZES_FUNCTION(uint64,  unsigned long)
SAMPLE_SIZES_FUNCTION(float32, float)
SAMPLE_SIZES_FUNCTION(float64, double)


inline void getSampleSizesForChunk(size_t* sampleSizes, void const* chunk, size_t chunkSize, size_t* stratification, int stratificationSize, enum TypeTag typeTag) {
  switch(typeTag) {
  case TYPE_ID_WRAPPER(INT8):    getSampleSizes_int8(   sampleSizes, chunk, chunkSize, stratification, stratificationSize); break;
  case TYPE_ID_WRAPPER(INT16):   getSampleSizes_int16(  sampleSizes, chunk, chunkSize, stratification, stratificationSize); break;
  case TYPE_ID_WRAPPER(INT32):   getSampleSizes_int32(  sampleSizes, chunk, chunkSize, stratification, stratificationSize); break;
  case TYPE_ID_WRAPPER(INT64):   getSampleSizes_int64(  sampleSizes, chunk, chunkSize, stratification, stratificationSize); break;
  case TYPE_ID_WRAPPER(UINT8):   getSampleSizes_uint8(  sampleSizes, chunk, chunkSize, stratification, stratificationSize); break;
  case TYPE_ID_WRAPPER(UINT16):  getSampleSizes_uint16( sampleSizes, chunk, chunkSize, stratification, stratificationSize); break;
  case TYPE_ID_WRAPPER(UINT32):  getSampleSizes_uint32( sampleSizes, chunk, chunkSize, stratification, stratificationSize); break;
  case TYPE_ID_WRAPPER(UINT64):  getSampleSizes_uint64( sampleSizes, chunk, chunkSize, stratification, stratificationSize); break;
  case TYPE_ID_WRAPPER(FLOAT32): getSampleSizes_float32(sampleSizes, chunk, chunkSize, stratification, stratificationSize); break;
  case TYPE_ID_WRAPPER(FLOAT64): getSampleSizes_float64(sampleSizes, chunk, chunkSize, stratification, stratificationSize); break;
  }
}

/* ------------------------------------------------------------------------------ */

int getSampleSizes(size_t* sampleSizes, void* source, IOConfiguration const* ioConfiguration, size_t* stratification, size_t sampleSize, int numStrati) {
  size_t chunkSize = CHUNK_SIZE;
  size_t numChunks = ioConfiguration->numElements / chunkSize;
  size_t elementSize = typeSizes[ioConfiguration->typeTag];

  /* First, get the number of elements contained in each stratum. */
  if(ioConfiguration->inputIsFile) {
    void* chunk = malloc(chunkSize * elementSize);
    for(size_t chunkIndex = 0; chunkIndex < numChunks; ++chunkIndex) {
      size_t read = fread(chunk, elementSize, chunkSize, (FILE*)source);
      assert(read == chunkSize); /* Assumes that file I/O does not fail. See argumentation below. */
      getSampleSizesForChunk(sampleSizes, chunk, chunkSize, stratification, numStrati, ioConfiguration->typeTag);
    }

    size_t rest = ioConfiguration->numElements - (numChunks * chunkSize);
    size_t read = fread(chunk, elementSize, rest, (FILE*)source);
    if(read < rest) {
      return 1;
    }
    getSampleSizesForChunk(sampleSizes, chunk, rest, stratification, numStrati, ioConfiguration->typeTag);

    free(chunk); // that's dirty!
  } else {
    char* chunk = source;

    for(size_t chunkIndex = 0; chunkIndex < numChunks; ++chunkIndex) {
      getSampleSizesForChunk(sampleSizes, chunk, chunkSize, stratification, numStrati, ioConfiguration->typeTag);
      chunk += chunkSize * elementSize;
    }

    size_t rest = ioConfiguration->numElements - (numChunks * chunkSize);
    getSampleSizesForChunk(sampleSizes, chunk, rest, stratification, numStrati, ioConfiguration->typeTag);
  }

  /* Secondly, compute the frequencies of the elements in each stratum and scale it by the sample size. 
   * The result are samples sizes for each individual stratum, which might be zero.                      */

  size_t sampleSizeSum = 0;
  for(int stratum = 0; stratum < numStrati; ++stratum) {
    sampleSizeSum += sampleSizes[stratum];
  }

  for(int stratum = 0; stratum < numStrati; ++stratum) {
    sampleSizes[stratum] = (size_t)( sampleSize * (double)sampleSizes[stratum] / sampleSizeSum );
  }
  return 0;
}

/* ---------------------------------------------------------------------------------------------------------------------------------------------------------- */

#define STRATUM_INDEX_FUNCTION(typenam,type)                                                                                                   \
inline void getStratumIndices_ ## typenam(size_t* stratumIdcs, void const* _chunk, size_t chunkSize, size_t* stratification, int numStrati) {  \
  type* const chunk = (type* const) _chunk;                                                                                                    \
  int chunkIndex;                                                                                                                              \
                                                                                                                                               \
  _Pragma(STRINGIFY(CONCATENATE(omp parallel for if(chunkSize >= 64))))                                                                        \
  for (chunkIndex = 0; chunkIndex < chunkSize; ++chunkIndex) {                                                                                 \
    double value = (double)chunk[chunkIndex];                                                                                                  \
                                                                                                                                               \
    for (size_t stratumIdx = 0; stratumIdx < numStrati; ++stratumIdx) {                                                                        \
      if (value <= stratification[stratumIdx]) {                                                                                               \
        stratumIdcs[chunkIndex] = stratumIdx;                                                                                                  \
        break;                                                                                                                                 \
      }                                                                                                                                        \
    }                                                                                                                                          \
  }                                                                                                                                            \
}

STRATUM_INDEX_FUNCTION(int8,    char)
STRATUM_INDEX_FUNCTION(int16,   short)
STRATUM_INDEX_FUNCTION(int32,   int)
STRATUM_INDEX_FUNCTION(int64,   long)
STRATUM_INDEX_FUNCTION(uint8,   unsigned char)
STRATUM_INDEX_FUNCTION(uint16,  unsigned short)
STRATUM_INDEX_FUNCTION(uint32,  unsigned int)
STRATUM_INDEX_FUNCTION(uint64,  unsigned long)
STRATUM_INDEX_FUNCTION(float32, float)
STRATUM_INDEX_FUNCTION(float64, double)


/* Stratified sampling on a single chunk */
void stratifiedSamplingOnChunk(char* sample, char* chunk, size_t chunkSize, size_t* stratification, int numStrati, size_t* sampleSizes, size_t* js, size_t* jsAfterSkip, unsigned char* skipElements, double* Ws, enum TypeTag typeTag, size_t* samplesPerStratum) {
  size_t elementSize = typeSizes[typeTag];
  char* c = chunk;

  size_t* stratumIndices = (size_t*)calloc(chunkSize, sizeof(size_t));
  switch(typeTag) {
  case TYPE_ID_WRAPPER(INT8):    getStratumIndices_int8(   stratumIndices, chunk, chunkSize, stratification, numStrati); break;
  case TYPE_ID_WRAPPER(INT16):   getStratumIndices_int16(  stratumIndices, chunk, chunkSize, stratification, numStrati); break;
  case TYPE_ID_WRAPPER(INT32):   getStratumIndices_int32(  stratumIndices, chunk, chunkSize, stratification, numStrati); break;
  case TYPE_ID_WRAPPER(INT64):   getStratumIndices_int64(  stratumIndices, chunk, chunkSize, stratification, numStrati); break;
  case TYPE_ID_WRAPPER(UINT8):   getStratumIndices_uint8(  stratumIndices, chunk, chunkSize, stratification, numStrati); break;
  case TYPE_ID_WRAPPER(UINT16):  getStratumIndices_uint16( stratumIndices, chunk, chunkSize, stratification, numStrati); break;
  case TYPE_ID_WRAPPER(UINT32):  getStratumIndices_uint32( stratumIndices, chunk, chunkSize, stratification, numStrati); break;
  case TYPE_ID_WRAPPER(UINT64):  getStratumIndices_uint64( stratumIndices, chunk, chunkSize, stratification, numStrati); break;
  case TYPE_ID_WRAPPER(FLOAT32): getStratumIndices_float32(stratumIndices, chunk, chunkSize, stratification, numStrati); break;
  case TYPE_ID_WRAPPER(FLOAT64): getStratumIndices_float64(stratumIndices, chunk, chunkSize, stratification, numStrati); break;
  }

  for(size_t i = 0; i < chunkSize; ++i, c += elementSize) {
    size_t stratumIndex = stratumIndices[i];

    size_t __js = js[stratumIndex];
    char* sampleForThisStratum = sample + samplesPerStratum[stratumIndex] * elementSize;

    if(__js < sampleSizes[stratumIndex]) {
      /* If the running indices are smaller than the individual sample sizes, we still need to fill up the reservoir. */
      memcpy(sampleForThisStratum + __js * elementSize, c, elementSize);
    } else {
      /* Reservoir sampling logic. */
      if (skipElements[stratumIndex]) {
        /* If skipElements is true, we need to compute the elements to skip. */
        jsAfterSkip[stratumIndex] = js[stratumIndex] + (size_t)(floor(log(rand01()) / log(1 - Ws[stratumIndex])) + 1);
        skipElements[stratumIndex] = 0;
      } else {
        /* Else, only do the sampling step if we have skipped the necessary elements.
         * Additionally, if the sample size is empty for a stratum, do nothing here. This case appears if
         * the sample size is rather small and the percentage of actual elements in the stratum of the
         * whole volume is small as well, such that the integer scaled sample size for this stratum is zero. */
        if (__js == jsAfterSkip[stratumIndex] && sampleSizes[stratumIndex] > 0) {
          size_t r = (size_t)(rand01() * (sampleSizes[stratumIndex] - 1));
          memcpy(sampleForThisStratum + r * elementSize, c, elementSize);
          Ws[stratumIndex] *= exp(log(rand01()) / sampleSizes[stratumIndex]);
          skipElements[stratumIndex] = 1;
        }
      }
    }

    js[stratumIndex]++;
  }

  free(stratumIndices);
}

/* ---------------------------------------------------------------------------------------------------------------------------------------------------------- */

/* ========================================================================================================================================================== */

enum RandomSamplingError extractSample(void* sample, IOConfiguration const* ioConfiguration, SamplingConfiguration const* configuration)
{
  assert(sample);
  assert(ioConfiguration);
  assert(configuration);

  size_t numSourceElements = ioConfiguration->numElements;

  if(numSourceElements <= configuration->sampleSize) {
    /* The dataset size is smaller than the sample size, thus we just retrieve the entire dataset. */

    if(ioConfiguration->inputIsFile) {
      FILE* blobFile = fopen(ioConfiguration->filename, "rb");
      if(!blobFile) return SAMPLING__SOURCE_FILE_CANNOT_BE_OPENED;

      if(fseek(blobFile, (long)ioConfiguration->headerOffset, SEEK_SET) != 0) {
        fclose(blobFile);
        return SAMPLING__COULD_NOT_JUMP_PAST_HEADER;
      }

      fread(sample, typeSizes[ioConfiguration->typeTag], ioConfiguration->numElements, blobFile);
    } else {
      memcpy(sample, ioConfiguration->blob, numSourceElements * typeSizes[ioConfiguration->typeTag]);
    }
    return SAMPLING__NO_ERROR;
  }

  setSeed(configuration->randomSeed ? *configuration->randomSeed : time(NULL));

  /* No stratification, i.e., we use the regular Algorithm L. */
  if(configuration->stratificationType == NONE) {
    return AlgorithmL(sample, ioConfiguration, configuration);
  }

  void* source;
  if(ioConfiguration->inputIsFile) {
    source = fopen(ioConfiguration->filename, "rb");
    if(!source) return SAMPLING__SOURCE_FILE_CANNOT_BE_OPENED;

    if(fseek(source, (long)ioConfiguration->headerOffset, SEEK_SET) != 0) {
      fclose(source);
      return SAMPLING__COULD_NOT_JUMP_PAST_HEADER;
    }
  } else {
    source = ioConfiguration->blob;
  }

  /* Of how much elements does a single chunk consist? */
  size_t chunkSize = CHUNK_SIZE;

  /* Else, get the stratification */
  size_t* stratification      = get_stratification(configuration, ioConfiguration->typeTag);

  /* The sizes of the samples for each stratum. */
  size_t* sampleSizes         = (size_t*)calloc(configuration->numberOfStrati, sizeof(size_t));
  /* Running indices for each stratum . */
  size_t* js                  = (size_t*)calloc(configuration->numberOfStrati, sizeof(size_t));
  /* The indicies after skipping due to Algorithm L for each stratum. */
  size_t* jsAfterSkip         = (size_t*)calloc(configuration->numberOfStrati, sizeof(size_t));
  /* Coefficients involved in computing how far to skip. */
  double* Ws                  = (double*)calloc(configuration->numberOfStrati, sizeof(double));
  /* Indices pointing to the beginning of the "samples per stratum". Expressed in byte offsets. */
  size_t* samplesPerStratum   = (size_t*)calloc(configuration->numberOfStrati, sizeof(size_t));
  /* Indicators if elements of the current stratum should be skipped */
  unsigned char* skipElements = (unsigned char*)calloc(configuration->numberOfStrati, sizeof(unsigned char));
  memset(skipElements, 1, configuration->numberOfStrati * sizeof(unsigned char)); /* internal note: works since type is unsigned char */

  int error = SAMPLING__NO_ERROR;

  /* First pass over the volume, get the counts for each stratum. */
  if(getSampleSizes(sampleSizes, source, ioConfiguration, stratification, configuration->sampleSize, configuration->numberOfStrati) != 0) {
    error = SAMPLING__COULD_NOT_READ_ELEMENT;
    goto CLEANUP;
  }

  /* Offsets to the samples per stratum. All samples per stratum are already inside the global sample 
     for memory efficiency and avoiding repeated copies.   */
  samplesPerStratum[0] = 0;
  for(size_t i = 1; i < configuration->numberOfStrati; ++i) {
    samplesPerStratum[i] = samplesPerStratum[i-1] + sampleSizes[i-1];
  }

  // Initialize stuff based on the sample sizes.
  for(int i = 0; i < configuration->numberOfStrati; ++i) {
    jsAfterSkip[i] = sampleSizes[i];
    Ws[i] = exp(log(rand01())/(sampleSizes[i] == 0 ? 1 : sampleSizes[i]));
  }

  /* Second pass over the volume, does the actual sampling.
     For efficiency reasons, we pursue a chunked approach.  */

  size_t numChunks = ioConfiguration->numElements / chunkSize;
  size_t elementSize = typeSizes[ioConfiguration->typeTag];

  char* chunk;
  if(ioConfiguration->inputIsFile) {
    // If we have a file, get back to the start of the binary blob.
    if(fseek((FILE*)source, (long)ioConfiguration->headerOffset, SEEK_SET) != 0) {
      fclose(source);
      error = SAMPLING__COULD_NOT_JUMP_PAST_HEADER;
      goto CLEANUP;
    }

    chunk = (char*)malloc(chunkSize * elementSize);
  } else {
    chunk = (char*)ioConfiguration->blob;
  }

  int const isNotFile = !ioConfiguration->inputIsFile;

  for(size_t chunkIndex = 0; chunkIndex < numChunks; ++chunkIndex) {
    if(ioConfiguration->inputIsFile) {
      size_t read = fread(chunk, elementSize, chunkSize, (FILE*)source);
      /* Normally, we would do some error handling here. However, the profiler detects this as quite slow.
       * Since we computed the chunk iterations such that each chunk is fully filled in this loop, the
       * only possibility of failure here is if the file I/O fails. Here, we assume that this does not happen.
       * If, however, this occurs, include the following error handling code block:
       *
       *   if(read < chunkSize) {
       *     error = SAMPLING__COULD_NOT_READ_ELEMENT;
       *     free(chunk);
       *     goto CLEANUP;
       *   }
       */
      assert(read == chunkSize);
    }

    stratifiedSamplingOnChunk(sample, chunk, chunkSize, stratification, configuration->numberOfStrati, sampleSizes, js, jsAfterSkip, skipElements, Ws, ioConfiguration->typeTag, samplesPerStratum);

    if(isNotFile) {
      chunk += chunkSize * elementSize;
    }
  }

  /* Do the stratified sampling over the remaining elements. */

  size_t rest = ioConfiguration->numElements - (numChunks * chunkSize);

  if(ioConfiguration->inputIsFile) {
    size_t read = fread(chunk, elementSize, rest, (FILE*)source);
    if(read < rest) {
      error = SAMPLING__COULD_NOT_READ_ELEMENT;
      free(chunk);
      goto CLEANUP;
    }
  }

  stratifiedSamplingOnChunk(sample, chunk, rest, stratification, configuration->numberOfStrati, sampleSizes, js, jsAfterSkip, skipElements, Ws, ioConfiguration->typeTag, samplesPerStratum);

  if(ioConfiguration->inputIsFile) free(chunk);

CLEANUP:
  if(ioConfiguration->inputIsFile) fclose(source);
  free(sampleSizes);
  free(js);
  free(jsAfterSkip);
  free(Ws);
  free(stratification);
  free(samplesPerStratum);
  free(skipElements);
  return error;
}
