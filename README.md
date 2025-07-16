# RDynLib

Creating and expanding an annotated spectral database generated from in house LC-MS/MS and LC-MSn experiments.

In this repository we restructure the existing dynlib database that was on csv format and we convert it to sql subdatabases in the following files:

-   "dynlib-to-rforms-ftms-neg.qmd" : in this file we convert the cvs format of the dynlib ftms negative database to a sql format

<!-- -->

-   "dynlib-to-rforms-ftms-pos.qmd" : in this file we convert the cvs format of the dynlib ftms positive database to a sql format

-   "dynlib-to-rforms-qtof-neg.qmd" : in this file we convert the cvs format of the dynlib qtof negative database to a sql format

-   "dynlib-to-rforms-qtof-pos.qmd" : in this file we convert the cvs format of the dynlib qtof positive database to a sql format

-   "DyL_SQL_to_MongoEmbadded2.qmd" : In this file, we created four MongoDB databases from the previous four SQL DynLib databases, along with a benchmark comparing the performance of both database types.

-   "Similarity" : this folder contains two sub folders in which we calculate the similarity between flax data and dynlib database, in the first subfolder "FTMS" we calculate the similarity between flax ftms negative data and dynlib ftms negative data, in addition we created a script for adding ftmsneg flax to the dynlib ftmsneg sql database. The second folder "QTOF" we calculate the similarity between flax qtof negative data and dynlib qtof negative data, in addition we created a script for adding qtofneg flax to the dynlib qtofneg sql database.
