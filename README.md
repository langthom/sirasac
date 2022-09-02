# sirasac
This repository forms the implementation of the *Simple Random Sampling Clustering (SiRaSaC)* library.
The goal of that library is to extract a (stratified) random sample from either an existing buffer or from a binary file.

The procedure is described in our paper [Clustering large 3D volumes: A sampling-based approach](TODO).

## Building sirasac
Building the C library is done via CMake, please check the options therein.
It also supports building Python bindings to the library, requiring Python 3.7 installed and the *wheel* and *numpy* packages.
The CMake target `INSTALL` will install the C library and the potentially generated Python wheel package to the specified `CMAKE_INSTALL_PREFIX`.

## Using sirasac
The usage of the C library directly is straight-forward. Below, we discuss the usage of the Python bindings.

The installation of the module is executed via *pip*:

```python
  pip install ../path/to/install/dir/sirasac-0.0.1-cp310-cp310-win_amd64.whl
```

where the file name depends on the Python version run (here: 3.10) and the platform (here: Windows, AMD64).


Within Python, one can import the sampling functions from the `sirasac` module.
In order to apply the sampling, the definition of a sampling configuration and a function extracing meta information is necessary.
The meta information is a triple describing the element type, the overall number of elements, and an offset to the beginning of the data in case of file-based sampling.
Additional information can be found in the Cython file.

The following example shows the setup for a `numpy.ndarray`:

```python
from sirasac import extract_random_sample_from_buffer

def meta_from_buf(buf):
  # returns index describing element type, the number of elements, and the offset
  vt = 'bhiqBHIQfd'
  return vt.index(buf.dtype.char), buf.size, 0


STRATIFICATION_TYPES = [ 'NONE', 'LINEAR', 'EXPONENTIAL', 'MIXED', 'PRECOMPUTED' ]

samplingConf = {
  'sampleSize':          1<<13,                                 # The number of elements to sample.
  'stratification_type': STRATIFICATION_TYPES.index('MIXED'),   # The stratification type to use.
  'number_of_strati':    7,                                     # The number of strati to use.
  'mixed_threshold':     5_000                                  # Threshold for the 'MIXED' stratification scheme.
}
```

Given that setup, the sampling from a numpy array is executed via

```python
buffer = np.arange(65535, dtype='H')
error_code, sample = extract_random_sample_from_buffer(buffer, meta_from_buf, samplingConf)
```

The Python bindings file contains the interpretations of the possible error codes.
The obtained sample could be used in any way, like to train some AI-based method.


## Citing sirasac
Please cite *SiRaSaC* as follows:

<p align="center">T. Lang, <em>Clustering large 3D volumes: A sampling-based approach</em>, 59th Annual British Conference on Non-Destructive Testing, 2022</p>


