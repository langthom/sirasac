#ifndef RANDOM_SAMPLING__H
#define RANDOM_SAMPLING__H

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

#include "sirasac.h"

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
enum RandomSamplingError extractSample(void* sample, IOConfiguration const* ioConfiguration, SamplingConfiguration const* configuration);

#endif // RANDOM_SAMPLING__H