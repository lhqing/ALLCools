:py:mod:`ALLCools.__main__`
===========================

.. py:module:: ALLCools.__main__

.. autoapi-nested-parse::

   CLI defined here

   When adding new function:
   1. add a func_register_subparser function to register the subparser
   2. add a condition in main func about this new func name, import the real func as func in main



Module Contents
---------------

.. py:data:: log
   

   

.. py:data:: DESCRIPTION
   :annotation: = Multiline-String

    .. raw:: html

        <details><summary>Show Value</summary>

    .. code-block:: text
        :linenos:

        
        The ALLCools command line toolkit contains multiple functions to manipulate the ALLC format, 
        a core file format that stores single base level methylation information.
        Throughout this toolkit, we use bgzip/tabix to compress and index the ALLC file to allow 
        flexible data query from the ALLC file.

        Current Tool List in ALLCools:

        [Generate ALLC]
        bam-to-allc          - Generate 1 ALLC file from 1 position sorted BAM file via 
                               samtools mpileup.

        [Manipulate ALLC]
        standardize-allc     - Validate 1 ALLC file format, standardize the chromosome names, 
                               compression format (bgzip) and index (tabix).
        tabix-allc           - A simple wrapper of tabix command to index 1 ALLC file.
        profile-allc         - Generate some summary statistics of 1 ALLC
        merge-allc           - Merge N ALLC files into 1 ALLC file
        extract-allc         - Extract information (strand, context) from 1 ALLC file

        [Get Region Level]
        allc-to-bigwig       - Generate coverage (cov) and ratio (mc/cov) bigwig track files 
                               from 1 ALLC file
        allc-to-region-count - Count region level mc, cov by genome bins or provided BED files.
        generate-mcds        - Generate methylation dataset (MCDS) for a group of ALLC file and 
                               different region sets. This is a convenient wrapper function for 
                               a bunch of allc-to-region-count and xarray integration codes. 
                               MCDS is inherit from xarray.DataSet
        generate-mcad        - Generate mCG hypo-methylation score AnnData dataset (MCAD) for 
                               a group of ALLC file and one region set.


    .. raw:: html

        </details>

   

.. py:data:: EPILOG
   :annotation: = Multiline-String

    .. raw:: html

        <details><summary>Show Value</summary>

    .. code-block:: text
        :linenos:

        
        Author: Hanqing Liu

        See ALLCools documentation here: https://lhqing.github.io/ALLCools/intro.html


    .. raw:: html

        </details>

   

.. py:class:: NiceFormatter(fmt=None, datefmt=None, style='%', validate=True)

   Bases: :py:obj:`logging.Formatter`

   From Cutadapt https://github.com/marcelm/cutadapt
   Do not prefix "INFO:" to info-level log messages (but do it for all other
   levels).
   Based on http://stackoverflow.com/a/9218261/715090 .

   .. py:method:: format(self, record)

      Format the specified record as text.

      The record's attribute dictionary is used as the operand to a
      string formatting operation which yields the returned string.
      Before formatting the dictionary, a couple of preparatory steps
      are carried out. The message attribute of the record is computed
      using LogRecord.getMessage(). If the formatting string uses the
      time (as determined by a call to usesTime(), formatTime() is
      called to format the event time. If there is exception information,
      it is formatted using formatException() and appended to the message.



.. py:function:: validate_environment()


.. py:function:: setup_logging(stdout=False, quiet=False, debug=False)

   From Cutadapt https://github.com/marcelm/cutadapt
   Attach handler to the global logger object


.. py:function:: _str_to_bool(v: str) -> bool


.. py:function:: bam_to_allc_register_subparser(subparser)


.. py:function:: standardize_allc_register_subparser(subparser)


.. py:function:: tabix_allc_register_subparser(subparser)


.. py:function:: profile_allc_register_subparser(subparser)


.. py:function:: merge_allc_register_subparser(subparser)


.. py:function:: extract_context_allc_register_subparser(subparser)


.. py:function:: allc_to_region_count_register_subparser(subparser)


.. py:function:: allc_to_bigwig_register_subparser(subparser)


.. py:function:: generate_mcds_register_subparser(subparser)


.. py:function:: generate_mcad_register_subparser(subparser)


.. py:function:: convert_mcds_to_zarr_register_subparser(subparser)


.. py:function:: main()


