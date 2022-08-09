:py:mod:`ALLCools._open`
========================

.. py:module:: ALLCools._open

.. autoapi-nested-parse::

   I/O Classes for ALLC and BAM files

   - read and write allc file
   - parallel writing
   - read region from allc
   - separate process gives better performance


   This file is modified from xopen 0.3.4
   Here is the licence

   Copyright (c) 2010-2018 Marcel Martin <mail@marcelm.net>

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in
   all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
   THE SOFTWARE.



Module Contents
---------------

.. py:data:: PIGZ
   :annotation: = True

   

.. py:data:: BGZIP
   :annotation: = True

   

.. py:class:: Closing

   Closeable class

   Inherit from this class and implement a close() method to offer context
   manager functionality.

   .. py:method:: close()
      :abstractmethod:


   .. py:method:: __enter__()


   .. py:method:: __exit__(*exc_info)


   .. py:method:: __del__()



.. py:class:: PipedGzipWriter(path, mode='wt', compresslevel=6, threads=1)

   Bases: :py:obj:`Closing`

   Piped Gzip Writer.

   Write gzip-compressed files by running an external gzip or pigz process and
   piping into it. On Python 2, this is faster than using gzip.open(). On
   Python 3, it allows to run the compression in a separate process and can
   therefore also be faster.

   .. py:method:: write(arg)


   .. py:method:: close()



.. py:class:: PipedGzipReader(path, region=None, mode='r')

   Bases: :py:obj:`Closing`

   Closeable class

   Inherit from this class and implement a close() method to offer context
   manager functionality.

   .. py:method:: close()


   .. py:method:: __iter__()


   .. py:method:: readline()


   .. py:method:: _raise_if_error()

      Raise an exception if the gzip process has exited with an error.

      Raise IOError if process is not running anymore and the
      exit code is nonzero.


   .. py:method:: read(*args)



.. py:class:: PipedBamReader(path, region=None, mode='r', include_header=True, samtools_parms_str=None)

   Bases: :py:obj:`Closing`

   Closeable class

   Inherit from this class and implement a close() method to offer context
   manager functionality.

   .. py:method:: close()


   .. py:method:: __iter__()


   .. py:method:: readline()


   .. py:method:: _raise_if_error()

      Raise IOError if process is not running anymore and the exit code is nonzero.


   .. py:method:: read(*args)



.. py:class:: PipedBamWriter(path, mode='wt', threads=1)

   Bases: :py:obj:`Closing`

   Closeable class

   Inherit from this class and implement a close() method to offer context
   manager functionality.

   .. py:method:: write(arg)


   .. py:method:: close()



.. py:function:: open_bam(file_path, mode='r', region=None, include_header=True, samtools_parms_str=None, threads=1)


.. py:function:: open_gz(file_path, mode='r', compresslevel=3, threads=1, region=None)


.. py:function:: open_allc(file_path, mode='r', compresslevel=3, threads=1, region=None)

   Open a .allc file.

   A replacement for the "open" function that can also open files that have
   been compressed with gzip, bzip2 or xz. If the file_path is '-', standard
   output (mode 'w') or input (mode 'r') is returned.

   The file type is determined based on the file_path: .gz is gzip, .bz2 is bzip2 and .xz is
   xz/lzma.

   When writing a gzip-compressed file, the following methods are tried in order to get the
   best speed 1) using a pigz (parallel gzip) subprocess; 2) using a gzip subprocess;
   3) gzip.open. A single gzip subprocess can be faster than gzip.open because it runs in a
   separate process.

   Uncompressed files are opened with the regular open().

   mode can be: 'rt', 'rb', 'at', 'ab', 'wt', or 'wb'. Also, the 't' can be omitted,
   so instead of 'rt', 'wt' and 'at', the abbreviations 'r', 'w' and 'a' can be used.

   threads is the number of threads for pigz. If None, then the pigz default is used.
   multi-thread only apply to writer, reader (decompression) can't be paralleled


.. py:function:: has_tabix(filename)


.. py:function:: has_bai(filename)


