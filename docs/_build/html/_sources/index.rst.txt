.. currentmodule:: MAPT

Microbiome Amplification Preference Tool
==========================================
Welcome! The goal of MAPT is to allow researchers greater flexibility when designing oligonucleotides to block unwanted sequences. We focus on PNAs, peptide nucleic acids, that are frequently designed to block host DNA but can be applied widely, as researchers see fit. With our tool scientists can add the primers they plan on using, and the sequence that they want to block as well as find out what sequences are at risk of being blocked unnecessarily. 

Check how we recapitulated the design of different PNAs!
---------------------------------------------------------
.. toctree::
   :maxdepth: 1

   Pages/gPNA
   Pages/mPNA
   Pages/pPNA


Design a PNA
----------------------------------------
PNA_Designer is the primary class for designing
your own custom PNA.

.. autoclass:: MAPT.PNA_Designer
   :members:

Additional Resources in this module
----------------------------------------
To increase the tractability offered to researchers, we
have made the k_mapper datatype available for users.

.. toctree::
   :maxdepth: 2

   Pages/extra



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

