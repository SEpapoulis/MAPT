.. pna_designer documentation master file, created by
   sphinx-quickstart on Mon Aug 12 13:35:40 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. currentmodule:: pna_designer

Welcome to pna_designer's documentation!
========================================
Something general about the PNA


Check how we recapitulated the design of different PNAs!
---------------------------------------------------------
.. toctree::
   :maxdepth: 2

   Pages/mPNA
   Pages/pPNA


Design a PNA
----------------------------------------
PNA_Designer is the
.. autoclass:: pna_designer.PNA_Designer
   :members:

Mapping k-mers
----------------------------------------
To increase the tractibility offered to researchers, we
have made the k_mapper datatype avalible for users.

.. autoclass:: pna_designer.k_mapper
    :members:

Managing Sequence Data from Silva
----------------------------------------
.. autoclass:: pna_designer.silva_manager
    :members:

.. toctree::
   :maxdepth: 2
   :caption: Contents:



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
