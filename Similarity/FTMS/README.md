---

---

**The FTMS folder contains the following files :**

-   "sim_functions_ftmsneg.R" : contains all the similarity functions.

-   "FTMS_neg_call.qmd" : in this file we call the similarity functions from the "sim_functions_ftmsneg.R" to calculate the similarity between flax ftms negative and Dynlib ftms negative data.

-   "test_real_spectra.qmd" : here we test the similarity functions from "sim_functions_ftmsneg.R" on some spectra examples and we compare the results with other similarity measures such as compareSpectra() from the spectra package, and NeutralLossesCosine() from the Python's matchms library.

-   "unit_test_ftmsneg" : We used the "testthat" library to do a unit test for the functions in the "sim_functions_ftmsneg.R" file.

-   "flax_ftms_to_Dynlib.qmd" : In this file we added the flax ftms negative data to the sql database of the Dynlib ftms negative.
