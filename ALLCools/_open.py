"""
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
"""
import codecs
import gzip
import os
import pathlib
import sys
import time
from subprocess import Popen, PIPE, run

try:
    run(['pigz', '--version'],
        stdout=PIPE,
        stderr=PIPE)
    PIGZ = True
except OSError:
    PIGZ = False

try:
    run(['bgzip', '--version'],
        stdout=PIPE,
        stderr=PIPE)
    BGZIP = True
except OSError:
    BGZIP = False


class Closing(object):
    """
    Inherit from this class and implement a close() method to offer context
    manager functionality.
    """

    def close(self):
        raise NotImplementedError

    def __enter__(self):
        return self

    def __exit__(self, *exc_info):
        self.close()

    def __del__(self):
        try:
            self.close()
        except OSError:
            pass


class PipedGzipWriter(Closing):
    """
    Write gzip-compressed files by running an external gzip or pigz process and
    piping into it. On Python 2, this is faster than using gzip.open(). On
    Python 3, it allows to run the compression in a separate process and can
    therefore also be faster.
    """

    def __init__(self, path, mode='wt', compresslevel=6, threads=1):
        """
        mode -- one of 'w', 'wt', 'wb', 'a', 'at', 'ab'
        compresslevel -- gzip compression level
        threads (int) -- number of pigz threads (None means to let pigz decide)
        """
        if mode not in ('w', 'wt', 'wb', 'a', 'at', 'ab'):
            raise ValueError("Mode is '{0}', but it must be 'w', 'wt', 'wb', 'a', 'at' or 'ab'".format(mode))

        self.outfile = open(path, mode)
        self.devnull = open(os.devnull, mode)
        self.closed = False
        self.name = path

        kwargs = dict(stdin=PIPE, stdout=self.outfile, stderr=self.devnull)
        # Setting close_fds to True in the Popen arguments is necessary due to
        # <http://bugs.python.org/issue12786>.
        # However, close_fds is not supported on Windows. See
        # <https://github.com/marcelm/cutadapt/issues/315>.
        if sys.platform != 'win32':
            kwargs['close_fds'] = True

        # only apply to gzip and pigz, bgzip use -l
        if 'w' in mode and compresslevel != 6:
            extra_args = ['-' + str(compresslevel)]
        else:
            extra_args = []

        if 'b' not in mode:
            encoding = None
        else:
            encoding = 'utf8'
        kwargs['encoding'] = encoding

        try:
            if BGZIP:
                bgzip_args = ['bgzip']
                if threads is not None and threads > 0:
                    bgzip_args += ['-@', str(threads), '-l', str(compresslevel)]
                self.process = Popen(bgzip_args, **kwargs)
                self.program = 'bgzip'
            elif PIGZ:
                pigz_args = ['pigz']
                if threads is not None and threads > 0:
                    pigz_args += ['-p', str(threads)]
                self.process = Popen(pigz_args + extra_args, **kwargs)
                self.program = 'pigz'
            else:
                # pigz not found, try regular gzip
                self.process = Popen(['gzip'] + extra_args, **kwargs)
                self.program = 'gzip'
        except OSError:
            self.outfile.close()
            self.devnull.close()
            raise
        self._file = codecs.getwriter('utf-8')(self.process.stdin)

    def write(self, arg):
        self._file.write(arg)

    def close(self):
        self.closed = True
        self._file.close()
        return_code = self.process.wait()
        self.outfile.close()
        self.devnull.close()
        if return_code != 0:
            raise IOError(f"Output {self.program} process terminated with exit code {return_code}")


class PipedGzipReader(Closing):
    # decompression can't be parallel even in pigz, so there is not thread/cpu parameter
    def __init__(self, path, region=None, mode='r'):
        if mode not in ('r', 'rt', 'rb'):
            raise ValueError(f"Mode is {mode}, but it must be 'r', 'rt' or 'rb'")
        if 'b' not in mode:
            encoding = 'utf8'
        else:
            encoding = None

        if region is None:
            self.process = Popen(['gzip', '-cd', path],
                                 stdout=PIPE,
                                 stderr=PIPE,
                                 encoding=encoding)
        else:
            self.process = Popen(['tabix', path] + region.split(' '),
                                 stdout=PIPE,
                                 stderr=PIPE,
                                 encoding=encoding)

        self.name = path
        self._file = self.process.stdout
        self._stderr = self.process.stderr
        self.closed = False
        # Give gzip a little bit of time to report any errors
        # (such as a non-existing file)
        time.sleep(0.01)
        self._raise_if_error()

    def close(self):
        self.closed = True
        return_code = self.process.poll()
        if return_code is None:
            # still running
            self.process.terminate()
        self._raise_if_error()

    def __iter__(self):
        for line in self._file:
            yield line
        self.process.wait()
        self._raise_if_error()

    def readline(self):
        return self._file.readline()

    def _raise_if_error(self):
        """
        Raise IOError if process is not running anymore and the
        exit code is nonzero.
        """
        return_code = self.process.poll()
        if return_code is not None and return_code != 0:
            message = self._stderr.read().strip()
            raise IOError(message)

    def read(self, *args):
        data = self._file.read(*args)
        if len(args) == 0 or args[0] <= 0:
            # wait for process to terminate until we check the exit code
            self.process.wait()
        self._raise_if_error()
        return data


class PipedBamReader(Closing):
    # decompression can't be parallel even in pigz, so there is not thread/cpu parameter
    def __init__(self, path, region=None, mode='r', include_header=True, samtools_parms_str=None):
        if mode not in ('r', 'rt', 'rb'):
            raise ValueError(f"Mode is {mode}, but it must be 'r', 'rt' or 'rb'")
        if 'b' not in mode:
            encoding = 'utf8'
        else:
            encoding = None

        command_list = ['samtools', 'view']
        if include_header:
            command_list.append('-h')
        if samtools_parms_str is not None:
            command_list.extend(samtools_parms_str.split(' '))

        if region is None:
            self.process = Popen(command_list + [path],
                                 stdout=PIPE,
                                 stderr=PIPE,
                                 encoding=encoding)
        else:
            self.process = Popen(command_list + [path] + region.split(' '),
                                 stdout=PIPE,
                                 stderr=PIPE,
                                 encoding=encoding)

        self.name = path
        self.file = self.process.stdout
        self._stderr = self.process.stderr
        self.closed = False
        # Give gzip a little bit of time to report any errors
        # (such as a non-existing file)
        time.sleep(0.01)
        self._raise_if_error()

    def close(self):
        self.closed = True
        return_code = self.process.poll()
        if return_code is None:
            # still running
            self.process.terminate()
        self._raise_if_error()

    def __iter__(self):
        for line in self.file:
            yield line
        self.process.wait()
        self._raise_if_error()

    def readline(self):
        return self.file.readline()

    def _raise_if_error(self):
        """
        Raise IOError if process is not running anymore and the
        exit code is nonzero.
        """
        return_code = self.process.poll()
        if return_code is not None and return_code != 0:
            message = self._stderr.read().strip()
            raise IOError(message)

    def read(self, *args):
        data = self.file.read(*args)
        if len(args) == 0 or args[0] <= 0:
            # wait for process to terminate until we check the exit code
            self.process.wait()
        self._raise_if_error()
        return data


class PipedBamWriter(Closing):
    def __init__(self, path, mode='wt', threads=1):
        """
            mode -- one of 'w', 'wt', 'wb', 'a', 'at', 'ab'
            threads (int) -- number of samtools threads
        """
        if mode not in ('w', 'wt', 'wb', 'a', 'at', 'ab'):
            raise ValueError("Mode is '{0}', but it must be 'w', 'wt', 'wb', 'a', 'at' or 'ab'".format(mode))

        self.outfile = open(path, mode)
        self.devnull = open(os.devnull, mode)
        self.closed = False
        self.name = path

        kwargs = dict(stdin=PIPE, stdout=self.outfile, stderr=self.devnull)
        # Setting close_fds to True in the Popen arguments is necessary due to
        # <http://bugs.python.org/issue12786>.
        # However, close_fds is not supported on Windows. See
        # <https://github.com/marcelm/cutadapt/issues/315>.
        if sys.platform != 'win32':
            kwargs['close_fds'] = True

        try:
            samtools_args = ['samtools', 'view', '-b', '-@', str(threads)]
            self.process = Popen(samtools_args, **kwargs)
            self.program = 'samtools'
        except OSError:
            self.outfile.close()
            self.devnull.close()
            raise
        self._file = codecs.getwriter('utf-8')(self.process.stdin)
        return

    def write(self, arg):
        self._file.write(arg)

    def close(self):
        self.closed = True
        self._file.close()
        return_code = self.process.wait()
        self.outfile.close()
        self.devnull.close()
        if return_code != 0:
            raise IOError(f"Output {self.program} process terminated with exit code {return_code}")


def open_bam(file_path, mode='r', region=None, include_header=True, samtools_parms_str=None, threads=1):
    if 'r' in mode:
        if region is not None:
            if not has_bai(file_path):
                raise FileNotFoundError(f'.bai index file not found for {file_path}. Index before region query.')
        try:
            return PipedBamReader(file_path, region=region, mode=mode,
                                  include_header=include_header,
                                  samtools_parms_str=samtools_parms_str)
        except OSError as e:
            raise e
    else:
        return PipedBamWriter(file_path, mode, threads=threads)


def open_gz(file_path, mode='r', compresslevel=3, threads=1, region=None):
    if region is not None:
        if not isinstance(region, str):
            raise TypeError('region parameter need to be string.')

    if 'r' in mode:
        try:
            return PipedGzipReader(file_path, region=region, mode=mode)
        except OSError:
            # gzip not installed
            return gzip.open(file_path, mode)
    else:
        try:
            return PipedGzipWriter(file_path, mode, compresslevel, threads=threads)
        except OSError:
            return gzip.open(file_path, mode, compresslevel=compresslevel)


def open_allc(file_path, mode='r', compresslevel=3, threads=1,
              region=None):
    """
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
    """
    file_path = pathlib.Path(file_path).resolve()
    file_exist = file_path.exists()
    if 'w' not in mode:
        if not file_exist:
            raise FileNotFoundError(f'{file_path} does not exist.')
    else:
        if not file_path.parent.exists():
            raise OSError(f'File directory {file_path.parent} does not exist.')
    file_path = str(file_path)

    if mode in ('r', 'w', 'a'):
        mode += 't'
    if mode not in ('rt', 'rb', 'wt', 'wb', 'at', 'ab'):
        raise ValueError("mode '{0}' not supported".format(mode))
    if compresslevel not in range(1, 10):
        raise ValueError("compresslevel must be between 1 and 9")
    if region is not None:
        # unzipped file
        if not file_path.endswith('gz'):
            raise ValueError(f'File must be compressed by bgzip to use region query. File path {file_path}')
        # normal gzipped file
        if not has_tabix(file_path):
            raise ValueError(f'Tried inspect {file_path}, '
                             'File is compressed by normal gzip, but region query only apply to bgzip')

        if not os.path.exists(file_path + '.tbi'):
            raise FileNotFoundError('region query provided, but .tbi index not found')

    if file_path.endswith('gz'):
        return open_gz(file_path, mode, compresslevel, threads, region=region)
    else:
        return open(file_path, mode)


def has_tabix(filename):
    tabix_path = filename + '.tbi'
    if os.path.exists(tabix_path):
        return True
    else:
        return False


def has_bai(filename):
    index_path = filename + '.bai'
    if os.path.exists(index_path):
        return True
    else:
        return False
