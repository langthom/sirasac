#ifndef SIRASAC__H
#define SIRASAC__H

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

#include "sirasac_export.h"

/* Wrap the names with a prefix to avoid double definitions in Cython context.
 * Specifically, the Windows basic signed types defines INT8 and many others.
 * These result in double definitions when compiling the Cython wrappers.     */
#define TYPE_ID_WRAPPER(type) type ## _T

/* \enum TypeTag
 * Enumeration of all element types.
 */
enum TypeTag {
  TYPE_ID_WRAPPER(INT8),  TYPE_ID_WRAPPER(INT16),  TYPE_ID_WRAPPER(INT32),  TYPE_ID_WRAPPER(INT64),
  TYPE_ID_WRAPPER(UINT8), TYPE_ID_WRAPPER(UINT16), TYPE_ID_WRAPPER(UINT32), TYPE_ID_WRAPPER(UINT64),
  TYPE_ID_WRAPPER(FLOAT32), TYPE_ID_WRAPPER(FLOAT64),
};


/* \enum RandomSamplingError */
enum RandomSamplingError {
  SAMPLING__NO_ERROR = 0,                    /**< No error. All is well.                                                      */
  SAMPLING__SOURCE_FILE_CANNOT_BE_OPENED,    /**< I/O: The binary source file cannot be opened.                               */
  SAMPLING__COULD_NOT_JUMP_PAST_HEADER,      /**< I/O: Cannot skip past the binary header in the file.                        */
  SAMPLING__COULD_NOT_FILL_RESERVOIR,        /**< I/O: Could not read enough elements to fill up the reservoir (Algorithm L). */
  SAMPLING__COULD_NOT_JUMP_TO_NEXT_ELEMENT,  /**< I/O: Could not jump to the next element to be read.                         */
  SAMPLING__COULD_NOT_READ_ELEMENT,          /**< I/O: Could not read the next element from the input data.                   */
};


/* \struct IOConfiguration
 * Describes how the input data can be read. Basically, the input can be either a
 * raw buffer, or a binary file.
 */
typedef struct {
  int inputIsFile;            /**< Indicates if the input is a file (if nonzero) or a raw buffer. */
  size_t numElements;         /**< Number of elements to sample from.                             */
  enum TypeTag typeTag;       /**< Type tag representing the input element type.                  */
  union {
    struct {
      size_t headerOffset;    /**< Offset inside the binary file to skip initially.               */
      const char* filename;   /**< Path to the binary file to read from.                          */
    };
    void* blob;               /**< Raw buffer to read from.                                       */
  };
} IOConfiguration;


/* \struct SamplingConfiguration
 * Describes how the random sampling is performed. Basically, the configuration chooses
 * between accelerated random sampling (stratification type \a NONE), and stratified random
 * sampling. In the latter case, a few different stratifications are available, where the
 * first three types (\a LINEAR, \a EXPONENTIAL, and \a MIXED) automatically compute a 
 * stratification, while the last available type (\a PRECOMPUTED) uses a priori provided
 * stratification values.
 */
typedef struct {
  size_t sampleSize;          /**< The number of sample elements to retrieve.         */
  int numberOfStrati;         /**< The number of strati the stratification considers. */
  unsigned int* randomSeed;   /**< Random seed. If NULL, no seed is set.              */

  /* \enum StratificationType
   * The available stratificationt types:
   * <ol>
   *   <li>NONE denotes no stratification, accelerated random sampling by means of Algorithm L.</li>
   *   <li>LINEAR stratification subdivides the grayscale value range evenly.</li>
   *   <li>EXPONENTIAL stratification subdivides the graysvale value range evenly in log space.</li>
   *   <li>MIXED stratification combines an exponential stratification of the interval up to some
           user-passed threshold, followed by a linear stratification of the remaining interval.</li>
   *   <li>PRECOMPUTED uses user-passed stratification values directly.</li>
   * </ol>
   */
  enum StratificationType {
    NONE, LINEAR, EXPONENTIAL, MIXED, PRECOMPUTED,
  } stratificationType;
  
  union {
    size_t mixedThreshold;  /**< Threshold that divides the MIXED stratification into exponential and linear parts. */
    double* strati;         /**< Application-dependent stratification values (incl. start and end values), used in PRECOMPUTED. */
  } stratification;
} SamplingConfiguration;


/* \brief Extracts a random sample from its input given some sampling configuration.
 * \param[out] sample           Target buffer to write the sample to. Must be properly allocated to store at least
 *                              SamplingConfiguration::sampleSize elements where each element requires the same 
 *                              number of bytes as in the read data.
 * \param[in]  ioConfiguration  I/O configuration specifying how the input data can be read.
 * \param[in]  configuration    Sampling configuration specifying the sampling process.
 * \return Returns an error code from the \a RandomSamplingError enumeration.
 * \sa RandomSamplingError
 * \sa IOConfiguration
 * \sa SamplingConfiguration
 */
sirasac_API int sirasac_sample(void* dest, IOConfiguration const* ioConfiguration, SamplingConfiguration const* samplingConfiguration);

/* \brief Computes the stratification used, if any.
 * This function allocates the stratification as a buffer of element type \a size_t on the heap. It is the
 * caller's responsibility to free the allocated memory via C's  \a free function.
 * \param[in]  configuration    Sampling configuration specifying the sampling process.
 * \param[in]  typeTag          Type tag representing the element type of a data item read.
 * \return Returns a pointer to an allocation containing the stratification values.
 */
sirasac_API size_t* get_stratification(SamplingConfiguration const* configuration, enum TypeTag typeTag);

#endif // SIRASAC__H