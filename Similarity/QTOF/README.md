------------------------------------------------------------------------

------------------------------------------------------------------------

**The QTOF folder contains the following files :**

-   "sim_functions_qtofneg.R" : contains all the similarity functions.

-   "QTOF_neg_call.qmd" : in this file we call the similarity functions from the "sim_functions_qtofneg.R" to calculate the similarity between flax qtof negative and Dynlib qtof negative data.

-   "test_real_spectra.qmd" : here we test the similarity functions from "sim_functions_qtofneg.R" on some spectra examples and we compare the results with other similarity measures such as compareSpectra() from the spectra package, and NeutralLossesCosine() from the Python's matchms library.

-   "unit_test_qtofneg.R" : We used the "testthat" library to do a unit test for the functions in the "sim_functions_qtofneg.R" file.

-   "flax_qtof_to_Dynlib.qmd" : In this file we added the flax qtof negative data to the sql database of the Dynlib qtof negative.

**The execution order:**

1.  "sim_functions_qtofneg.R"

2.  "QTOF_neg_call.qmd"

3.  "unit_test_qtofneg.R"

4.  "test_real_spectra.qmd"

5.  "flax_qtof_to_Dynlib.qmd"
