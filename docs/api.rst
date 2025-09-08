API Reference
=============

This page contains the comprehensive API documentation for the set of modules, each 
collected into their own section below. 

Batch Processing
-----------------------
The batch modules handle the submission of batch jobs to the cluster, as well as the
conversion of the ROOT and AmpTools fit results into csv files for later analysis.

.. autosummary::
   :toctree: generated
   :recursive:

   neutralb1.batch

Analysis
---------------
These modules provide statistical, preprocessing, plotting, and other analysis tools
for analyzing the csv data files produced by the batch processing modules.

.. autosummary::
   :toctree: generated
   :recursive:

   neutralb1.analysis


Utilities
---------
This singular module contains a variety of functions applicable across each stage of the
analysis.

.. autosummary::
   :toctree: generated
   :recursive:

   neutralb1.utils

Selection Module
----------------

.. autosummary::
   :toctree: generated
   :recursive:

   neutralb1.selection