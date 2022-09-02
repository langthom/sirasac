
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "../sirasac.h"

void sampleFromFile(const char* dataFile, enum StratificationType type, int numberOfStrati, int sampleSize, int numberOfVoxels, int dtype, size_t headerOffset, double* stratificationValues) {
  enum TypeTag __data_type_codes[] = {
    [7] = TYPE_ID_WRAPPER(INT8),[15] = TYPE_ID_WRAPPER(INT16),[31] = TYPE_ID_WRAPPER(INT32),[63] = TYPE_ID_WRAPPER(INT64),
    [8] = TYPE_ID_WRAPPER(UINT8),[16] = TYPE_ID_WRAPPER(UINT16),[32] = TYPE_ID_WRAPPER(UINT32),[64] = TYPE_ID_WRAPPER(UINT64),
    [24] = TYPE_ID_WRAPPER(FLOAT32),[53] = TYPE_ID_WRAPPER(FLOAT64),
  };

  size_t __data_type_sizes[] = {
    [7] = sizeof(char),[15] = sizeof(short),[31] = sizeof(int),[63] = sizeof(long),
    [8] = sizeof(unsigned char),[16] = sizeof(unsigned short),[32] = sizeof(unsigned int),[64] = sizeof(unsigned long),
    [24] = sizeof(float),[53] = sizeof(float),
  };

  void* dest = malloc(sampleSize * __data_type_sizes[dtype]);

  unsigned int seed = 42;
  enum TypeTag typeTag = __data_type_codes[dtype];
  size_t voxelSizeInBytes = __data_type_sizes[dtype];

  SamplingConfiguration configuration;
  configuration.sampleSize = sampleSize;
  configuration.randomSeed = &seed;
  configuration.stratificationType = type;
  configuration.numberOfStrati = numberOfStrati;
  if (stratificationValues != NULL) {
    configuration.stratification.strati = stratificationValues;
  } else {
    configuration.stratification.mixedThreshold = (size_t)(0.076 * (1ull << (1ull << (3 * voxelSizeInBytes))));
  }

  IOConfiguration ioConf;
  ioConf.inputIsFile  = 1;
  ioConf.numElements  = numberOfVoxels;
  ioConf.typeTag      = typeTag;
  ioConf.headerOffset = headerOffset;
  ioConf.filename     = dataFile;
  
  clock_t start = clock();
  int err = sirasac_sample(dest, &ioConf, &configuration);
  clock_t end = clock();

  double duration = 1000.0 * (end - start) / CLOCKS_PER_SEC;
  printf("Sampling took %d [ms].\n", (int)duration);

  if(err != SAMPLING__NO_ERROR) {
    printf(">>> Sampling failed, error code is %d\n", err);
  } else if(configuration.stratificationType != NONE) {
    size_t* stratification = get_stratification(&configuration, typeTag);

    printf("Stratification: ");
    for(int i = 0; i < configuration.numberOfStrati; ++i) printf("%zd ", stratification[i]);
    printf("\n");

    double* values_per_stratum = (double*)calloc(configuration.numberOfStrati, sizeof(double));

    unsigned char* s = (unsigned char*)dest;

    for(int i = 0; i < sampleSize; ++i, s += voxelSizeInBytes) {
      double valueNum = 0;
      switch(typeTag) {
      case TYPE_ID_WRAPPER(INT8):    valueNum = (double)(*(char*)          s); break;
      case TYPE_ID_WRAPPER(INT16):   valueNum = (double)(*(short*)         s); break;
      case TYPE_ID_WRAPPER(INT32):   valueNum = (double)(*(int*)           s); break;
      case TYPE_ID_WRAPPER(INT64):   valueNum = (double)(*(long*)          s); break;
      case TYPE_ID_WRAPPER(UINT8):   valueNum = (double)(*(unsigned char*) s); break;
      case TYPE_ID_WRAPPER(UINT16):  valueNum = (double)(*(unsigned short*)s); break;
      case TYPE_ID_WRAPPER(UINT32):  valueNum = (double)(*(unsigned int*)  s); break;
      case TYPE_ID_WRAPPER(UINT64):  valueNum = (double)(*(unsigned long*) s); break;
      case TYPE_ID_WRAPPER(FLOAT32): valueNum = (double)(*(float*)         s); break;
      case TYPE_ID_WRAPPER(FLOAT64): valueNum = (double)(*(double*)        s); break;
      }

      for(int j = 0; j < configuration.numberOfStrati; ++j) {
        if(valueNum <= stratification[j]) {
          values_per_stratum[j] += valueNum;
          values_per_stratum[j+5]++;
          break;
        }
      }
    }

    free(values_per_stratum);
    free(stratification);

    printf("Means of stratum samples: ");
    for(int i = 0; i < 5; ++i) {
      printf("%f ", values_per_stratum[i] == 0 ? 0 : values_per_stratum[i] / values_per_stratum[i+5]);
    }
    printf("\n");
  }

  free(dest);
}


void sampleFromBuffer(enum StratificationType type, int numberOfStrati, int sampleSize, double* stratificationValues) {
  int N = 1ull << 16;
  unsigned short* data = (unsigned short*)malloc(N * sizeof(unsigned short));
  srand(42);
  for (int i = 0; i < N; ++i) {
    data[i] = (unsigned short)(rand() * 65535 / RAND_MAX);
  }

  void* dest = malloc(sampleSize * sizeof(unsigned short));

  unsigned int seed = 42;
  enum TypeTag typeTag = TYPE_ID_WRAPPER(UINT16);
  size_t voxelSizeInBytes = sizeof(unsigned short);

  SamplingConfiguration configuration;
  configuration.sampleSize = sampleSize;
  configuration.randomSeed = &seed;
  configuration.stratificationType = type;
  configuration.numberOfStrati = numberOfStrati;
  if (stratificationValues != NULL) {
    configuration.stratification.strati = stratificationValues;
  } else {
    configuration.stratification.mixedThreshold = (size_t)(0.076 * (1ull << (1ull << (3 * voxelSizeInBytes))));
  }

  IOConfiguration ioConf;
  ioConf.inputIsFile  = 0;
  ioConf.numElements  = N;
  ioConf.typeTag      = typeTag;
  ioConf.headerOffset = 0;
  ioConf.blob         = data;
  
  clock_t start = clock();
  int err = sirasac_sample(dest, &ioConf, &configuration);
  clock_t end = clock();

  double duration = 1000.0 * (end - start) / CLOCKS_PER_SEC;
  printf("Sampling took %d [ms].\n", (int)duration);

  if(err != SAMPLING__NO_ERROR) {
    printf(">>> Sampling failed, error code is %d\n", err);
  } else if(configuration.stratificationType != NONE) {
    size_t* stratification = get_stratification(&configuration, typeTag);

    printf("Stratification: ");
    for(int i = 0; i < configuration.numberOfStrati; ++i) printf("%zd ", stratification[i]);
    printf("\n");

    double* values_per_stratum = (double*)calloc(configuration.numberOfStrati, sizeof(double));

    unsigned char* s = (unsigned char*)dest;

    for(int i = 0; i < sampleSize; ++i, s += voxelSizeInBytes) {
      double valueNum = 0;
      switch(typeTag) {
      case TYPE_ID_WRAPPER(INT8):    valueNum = (double)(*(char*)          s); break;
      case TYPE_ID_WRAPPER(INT16):   valueNum = (double)(*(short*)         s); break;
      case TYPE_ID_WRAPPER(INT32):   valueNum = (double)(*(int*)           s); break;
      case TYPE_ID_WRAPPER(INT64):   valueNum = (double)(*(long*)          s); break;
      case TYPE_ID_WRAPPER(UINT8):   valueNum = (double)(*(unsigned char*) s); break;
      case TYPE_ID_WRAPPER(UINT16):  valueNum = (double)(*(unsigned short*)s); break;
      case TYPE_ID_WRAPPER(UINT32):  valueNum = (double)(*(unsigned int*)  s); break;
      case TYPE_ID_WRAPPER(UINT64):  valueNum = (double)(*(unsigned long*) s); break;
      case TYPE_ID_WRAPPER(FLOAT32): valueNum = (double)(*(float*)         s); break;
      case TYPE_ID_WRAPPER(FLOAT64): valueNum = (double)(*(double*)        s); break;
      }

      for(int j = 0; j < configuration.numberOfStrati; ++j) {
        if(valueNum <= stratification[j]) {
          values_per_stratum[j] += valueNum;
          values_per_stratum[j+5]++;
          break;
        }
      }
    }

    free(values_per_stratum);
    free(stratification);

    printf("Means of stratum samples: ");
    for(int i = 0; i < 5; ++i) {
      printf("%f ", values_per_stratum[i] == 0 ? 0 : values_per_stratum[i] / values_per_stratum[i+5]);
    }
    printf("\n");
  }

  free(dest);
  free(data);
}


int main(int argc, char** argv)
{
  if (argc < 6) {
    printf("Usage: %s <file> <stratification-type> <number-of-strati> <sample-size> <number-of-voxels> [<header-offset = 2048> | <stratification-values>...]\n\n", argv[0]);
    printf("Options:\n");
    printf("  <file>                    Full path to the data file containing the binary voxel data.\n");
    printf("  <stratification-type>     Which stratification to apply.\n");
    printf("                            Must be an integer between 0 and 4 inclusive, where ...\n");
    printf("                              * ... 0 corresponds to _no_ stratification, regular random sampling is applied. \n");
    printf("                              * ... 1 corresponds to stratification using a linear subdivision of the grayscale value range.\n");
    printf("                              * ... 2 corresponds to stratification using an exponential subdivision of the grayscale value range.\n");
    printf("                              * ... 3 corresponds to stratification using a mixed subdivision of the grayscale value range (exponential followed by linear).\n");
    printf("                              * ... 4 corresponds to stratification using a precomputed subdivision of the grayscale value range.\n");
    printf("                                    In this case, all numbers following the last command line paramter will be interpreted as stratification interval limits between 0 and 1 (floating point).\n");
    printf("  <number-of-strati>        Number of strati to apply. Typically a value of 5 or 7 suffices.\n");
    printf("  <sample-size>             Number of voxels to sample, typically 4096 elements or similar.\n");
    printf("  <number-of-voxels>        Number of voxels contained in the volume.\n");
    printf("  <dtype>                   Data type of each individual voxel. Must be one of [7,8,15,16,31,32,63,64,24,53].\n");
    printf("  <header-offset>           Offset in bytes to skip over initially in the passed data file. Defaults to 2048\n");
    printf("  <stratification-values>   Optional list of stratification upper bin limits. Typically, the last value is 1.\n");
    return EXIT_FAILURE;
  }

  int stratificationType = atoi(argv[2]);
  int headerOffset = 2048;
  double* stratificationValues = NULL;

  if (stratificationType != 4) {
    if (argc == 8) {
      headerOffset = atoi(argv[7]);
    }
  } else {
    int numberOfStratificationValues = argc - 7;
    stratificationValues = (double*)malloc(numberOfStratificationValues * sizeof(double));
    for (int i = 0; i < numberOfStratificationValues; ++i) {
      sscanf(argv[i+7], "%lf", &stratificationValues[i]);
    }
  }

  if (0) {
    sampleFromBuffer(stratificationType, atoi(argv[3]), atoi(argv[4]), stratificationValues);
  } else {
    sampleFromFile(argv[1], stratificationType, atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), atoi(argv[6]), headerOffset, stratificationValues);
  }


  if (stratificationValues) {
    free(stratificationValues);
  }
  return EXIT_SUCCESS;
}