library(testthat)
source("C:/Users/amentag/Desktop/these/These/clones/RDynlib/RDynLib/Similarity/sim_functions.R")

a <- cbind(mz = c(13.3, 15.4, 15.6, 23.4, 26.2),
           intensity = 1:5)
b <- cbind(mz = c(15.3, 33.4),
           intensity = 1:2)

test_that("remove_duplicates works", {
  res <- remove_duplicates(a[, 1], a[, 2])
  expect_true(is.list(res))
  expect_equal(names(res), c("mz", "intensity"))
  ## Here I'm not sure. Would maybe have preferred that 15.4 and 15.6 get
  ## merged into the same? i.e. result would be rather:
  ## c(13, 15, 23, 26), c(1, 3, 4, 5)?
  expect_equal(res$mz, c(13, 15, 16, 23, 26))
  expect_equal(res$intensity, c(1, 2, 3, 4, 5))
})

test_that("dynlib_map works", {
  
  ## Empty input. Expect empty result, but in a specific format
  res <- dynlib_map(matrix(numeric(), nrow = 0), matrix(numeric(), nrow = 0))
  expect_true(is.list(res))
  expect_equal(names(res), c("x", "y"))
  expect_true(nrow(res$x) == 0)
  expect_true(nrow(res$y) == 0)
  
  ## Only direct matches, no neutral loss matches
  res <- dynlib_map(a, b, 12, 12)
  expect_true(is.list(res))
  expect_equal(names(res), c("x", "y"))
  x <- res$x
  expect_true(attributes(x)$wintensity_sum > 0)
  attributes(x)$wintensity_sum <- NULL
  y <- res$y
  expect_true(attributes(y)$wintensity_sum > 0)
  attributes(y)$wintensity_sum <- NULL
  expect_equal(x, cbind(mz = 15, intensity = 2))
  expect_equal(y, cbind(mz = 15, intensity = 1))
  
  ## Direct and neutral loss matches
  res <- dynlib_map(a, b, 0, 10)
  expect_true(is.list(res))
  expect_equal(names(res), c("x", "y"))
  x <- res$x
  expect_true(attributes(x)$wintensity_sum > 0)
  attributes(x)$wintensity_sum <- NULL
  expect_equal(x, cbind(mz = c(15, 23), intensity = c(2, 4)))
  y <- res$y
  expect_true(attributes(y)$wintensity_sum > 0)
  attributes(y)$wintensity_sum <- NULL
  expect_equal(y, cbind(mz = c(15, 33), intensity = c(1, 2)))
  
  ## No matches
  d <- cbind(mz = c(33.3, 46.1), intensity = c(1, 2))
  res <- dynlib_map(a, d, 0, 0)
  expect_true(is.list(res))
  expect_equal(names(res), c("x", "y"))
  x <- res$x
  expect_true(attributes(x)$wintensity_sum > 0)
  attributes(x)$wintensity_sum <- NULL
  y <- res$y
  expect_true(attributes(y)$wintensity_sum > 0)
  attributes(y)$wintensity_sum <- NULL
  expect_equal(x, cbind(mz = numeric(), intensity = numeric()))
  expect_equal(y, cbind(mz = numeric(), intensity = numeric()))
  
  ## Only neutral loss matches
  res <- dynlib_map(a, d, 0, 20)
  expect_true(is.list(res))
  expect_equal(names(res), c("x", "y"))
  x <- res$x
  expect_true(attributes(x)$wintensity_sum > 0)
  attributes(x)$wintensity_sum <- NULL
  y <- res$y
  expect_true(attributes(y)$wintensity_sum > 0)
  attributes(y)$wintensity_sum <- NULL
  expect_equal(x, cbind(mz = c(13, 26), intensity = c(1, 5)))
  expect_equal(y, cbind(mz = c(33, 46), intensity = c(1, 2)))
})




test_that("dynlib_symmetric_dotproduct works", {
  x <- cbind(mz = c(15, 17),
             intensity = c(2, 19))
  y <- cbind(mz = c(15, 17),
             intensity = c(2, 19))
  m <- dynlib_map(x, y, 0, 0)
  
  res <- dynlib_symmetric_dotproduct(m$x, m$y)
  expect_equal(res, 1)
  
  x <- cbind(mz = c(10, 17, 5, 9),
             intensity = c(2, 19, 6, 5))
  y <- cbind(mz = c(10, 17, 5, 9),
             intensity = c(2, 19, 7, 8))
  m <- dynlib_map(x, y, 0, 0)
  res <- dynlib_symmetric_dotproduct(m$x, m$y)
  expect_equal(res, 0.9984423)
  
  x <- cbind(mz = c(10, 17, 5, 9),
             intensity = c(8, 9, 6, 5))
  y <- cbind(mz = c(10, 17, 5, 9),
             intensity = c(2, 19, 7, 8))
  m <- dynlib_map(x, y, 0, 0)
  res <- dynlib_symmetric_dotproduct(m$x, m$y)
  expect_equal(res, 0.7850193)
})



#First test of the sim_final_final() function

x <- DataFrame(
               precursorMz = 0)
x$mz <- list(c(15,17))
x$intensity = list(c(2, 19))
spx <- Spectra(x)

y <- DataFrame(
  precursorMz = 0)
y$mz <- list(c(15,17))
y$intensity = list(c(2, 19))
spy <- Spectra(y)

res <- sim_final_final(spx, spy)
res


#second test

spd <- DataFrame(
  msLevel = c(2L, 2L, 2L),
  polarity = c(0L, 0L, 0L),
  precursorMz = 0,
  id = c("HMDB0000001", "HMDB0000001", "HMDB0001847"),
  name = c("1-Methylhistidine", "1-Methylhistidine", "Caffeine"))

## Assign m/z and intensity values.
spd$mz <- list(
  c(109.2, 124.2, 124.5, 170.16, 170.52),
  c(83.1, 96.12, 97.14, 109.14, 124.08, 125.1, 170.16),
  c(56.0494, 69.0447, 83.0603, 109.0395, 110.0712,
    111.0551, 123.0429, 138.0662, 195.0876)
)
spd$intensity <- list(
  c(3.407, 47.494, 3.094, 100.0, 13.240),
  c(6.685, 4.381, 3.022, 16.708, 100.0, 4.565, 40.643),
  c(0.459, 2.585, 2.446, 0.508, 8.968, 0.524, 0.974, 100.0, 40.994)
 )

spd <- Spectra(spd)

spc <- DataFrame(
  msLevel = c(2L, 2L, 2L),
  polarity = c(0L, 0L, 0L),
  precursorMz = 0,
  id = c("HMDB0000001", "HMDB0000001", "HMDB0001847"),
  name = c("1-Methylhistidine", "1-Methylhistidine", "Caffeine"))

## Assign m/z and intensity values.
spc$mz <- list(
  c(10.2, 1233.2, 124.5, 170.16, 170.52),
  c(83.1, 96.12, 97.14, 109.14, 124.08, 125.1, 170.16),
  c(56.0494, 69.0447, 83.0603, 109.0395, 110.0712,
    111.0551, 123.0429, 138.0662, 195.0876)
  )
spc$intensity <- list(
  c(3.407, 47.494, 88.094, 100.0, 13.240),
  c(6.685, 4.381, 55.022, 16.708, 100.0, 4.565, 40.643),
  c(0.459, 2.585, 2.446, 0.508, 8.968, 0.524, 0.974, 100.0, 40.994)
 )

spc <- Spectra(spc)

res2 <- sim_final_final(spd, spc)









