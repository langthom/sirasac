# MIT License
# 
# Copyright (c) 2022 Dr. Thomas Lang
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import numpy as np
cimport numpy as np
cimport c_random_sampling
cimport cython

from libc.stdint cimport uint8_t

# First, define the API interface in Cython. Basically, defines the necessary structs and declares the sampling function.

cdef extern from "sirasac.h":
  cdef enum TypeTag:
    pass
  
  ctypedef struct IOConfiguration:
    int inputIsFile
    Py_ssize_t numElements
    TypeTag typeTag
    Py_ssize_t headerOffset
    const char* filename
    void* blob
  
  ctypedef union stratification_conf:
    Py_ssize_t mixedThreshold
    double* strati
  
  ctypedef struct SamplingConfiguration:
    Py_ssize_t sampleSize
    int numberOfStrati
    int stratificationType
    unsigned int* randomSeed
    stratification_conf stratification
  
  int sirasac_sample(void* dest, IOConfiguration* ioConfiguration, SamplingConfiguration* samplingConfiguration)


# Extracts a random sample from the given configuration by calling the fast C implementation.
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef extract_random_sample(IOConfiguration* ioConf, py_samplingConf):
  dtypes = "bhiqBHIQfd"
  dtype = dtypes[<int>ioConf.typeTag]
  
  cdef SamplingConfiguration samplingConf
  samplingConf.sampleSize         = py_samplingConf.get('sampleSize',           4096)
  samplingConf.stratificationType = py_samplingConf.get('stratification_type',     0)
  samplingConf.numberOfStrati     = py_samplingConf.get('number_of_strati',        5)
  
  cdef unsigned int seed = 0
  if 'seed' in py_samplingConf:
    seed = py_samplingConf['seed']
    samplingConf.randomSeed = &seed
  else:
    samplingConf.randomSeed = NULL
  
  cdef double[::1] strati
  
  if samplingConf.stratificationType == 3: # MIXED
    samplingConf.stratification.mixedThreshold = py_samplingConf.get('mixed_threshold', 5_000)
  elif samplingConf.stratificationType == 4: # PRECOMPUTED
    _precomputed_stratification = py_samplingConf['strati']
    strati = _precomputed_stratification
    samplingConf.stratification.strati = &strati[0]
  
  cdef np.ndarray[uint8_t, ndim=1, mode="c"] sample = np.zeros(samplingConf.sampleSize * np.dtype(dtype).itemsize, dtype='B', order='C')
  
  error_code = sirasac_sample(<void*>&sample[0], ioConf, &samplingConf)
  
  return error_code, np.empty(0,dtype=dtype) if error_code != 0 else np.frombuffer(sample, dtype=dtype)

################################################################
### The actual interface implementations
################################################################

#----------------------------------------------
# Extract sample from file.

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef extract_random_sample_from_file(str filename, extract_metadata_func, samplingConf):
  """
  Extracts a (potentially stratified) random sample of a fixed size from a file containing the binary data.
  
  Arguments:
    filename : str
      Full path to the file containing the binary data from which the sample shall be extracted.
      The file contents may begin with some header information, followed by binary data only.
    
    extract_metadata_func : function
      Custom function responsible for extracting metadata about the data in the passed file.
      Specifically, this function must return a type tag, the number of elements (not bytes!) to read
      and the header offset (in bytes) to skip initially when reading the file.
      
      The type tag needs to be one of the following numeric values corresponding to the datatype of each read element.
        * 0 .. int8
        * 1 .. int16
        * 2 .. int32
        * 3 .. int64
        * 4 .. uint8
        * 5 .. uint16
        * 6 .. uint32
        * 7 .. uint64
        * 8 .. float32
        * 9 .. float64
    
    samplingConf : dict
      A dictionary containing information how to sample from the data. The following keys must be present:
        * sampleSize          .. The number of elements to sample from the input data.
                                 If the input data size is smaller than this number, everything will be sampled.
                                 Defaults to 4096.
        * stratification_type .. Type of stratification to use. A stratification subdivides the available grayscale value
                                 range (which is dependent on the element data type) into bins, and the stratification type
                                 specifies how this subdivison takes place.
                                 On a general note, since the grayscale value range is well-defined for integral types only,
                                 one should use the PRECOMPUTED stratification in case of floating point values.
                                 
                                 Must be one of the following values:
                                   0 .. No stratification. A random sample according to Algorithm L will be extracted.
                                   1 .. LINEAR stratification. The value range is divided evenly into bins of equal width.
                                   2 .. EXPONENTIAL stratification. The value range is divided first in half, then the smaller
                                        value range half is divided again in half, and so on until the number of bins is reached.
                                        Intended to reflect some kind of exponential value decay.
                                   3 .. MIXED stratification. Uses an even subdivision in log10 space until a certain
                                        value (see argument 'mixed_threshold'), followed by a linear subdivision on the remaining range.
                                   4 .. PRECOMPUTED stratification, uses user-defined values. See argument 'strati'.
                                 This parameter defaults to 0.
        * number_of_strati    .. The number of strati by which the grayscale value shall be divided into. Defaults to 5.
        * seed                .. Seed for random number generation. Optional, if missing no seed will be set.
        * mixed_threshold     .. Threshold at which grayscale value the transition from logarithmic to linear subdivision takes place.
                                 This value is only read if the stratification type is "MIXED". Defaults to 5000 then.
        * strati              .. Fixed list of stratification / bin upper values. Read only if the stratification type is "PRECOMPUTED".
                                 If used, the last value typically should be either the largest representable value of the 
                                 element data type, or alternatively the maximum voxel value. Any value bigger than the last value will be ignored.
  Returns:
    Returns the extracted sample as a numpy array
  """
  _fname_encoded = filename.encode('utf-8')
  cdef const char* c_fname = _fname_encoded
  typeTag, numElements, headerOffset = extract_metadata_func(filename)
  
  cdef IOConfiguration ioConf
  ioConf.inputIsFile  = 1
  ioConf.numElements  = numElements
  ioConf.headerOffset = headerOffset
  ioConf.typeTag      = typeTag
  ioConf.filename     = c_fname
  
  return extract_random_sample(&ioConf, samplingConf)


# --------------------------------------------
# Extract sample from buffer.

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef extract_random_sample_from_buffer(np.ndarray buffer, extract_metadata_func, samplingConf):
  """
  Extracts a (potentially stratified) random sample of a fixed size from a file containing the binary data.
  This function performs exactly like 'extract_random_sample_from_file' except there is not file involved.
  For a detailed description of its arguments, see the documentation of the mentioned function.
  """
  
  typeTag, numElements, headerOffset = extract_metadata_func(buffer)
  
  # Interpret the Numpy array as a memory view of unsigned char type, i.e., its raw bytes. 
  # This is necessary to successfully cast it as  void*  below.
  cdef np.ndarray[uint8_t, ndim=1, mode='c'] _buffer = buffer.view('B')

  cdef IOConfiguration ioConf
  ioConf.inputIsFile  = 0
  ioConf.numElements  = numElements
  ioConf.headerOffset = headerOffset
  ioConf.typeTag      = typeTag
  ioConf.blob         = <void*> &_buffer[0]
  
  return extract_random_sample(&ioConf, samplingConf)


