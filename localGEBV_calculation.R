####load dependencies####

library(purrr)
library(furrr)
library(future)
library(parallel)
library(dplyr)
library(progressr)

load(file = "Example_Files/gapit_marker_effects.R")
load(file = "Example_Files/gapit_marker_pecov.R")
