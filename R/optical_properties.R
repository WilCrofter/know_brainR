# Utilities for data from
# Steven L Jacques, Optical properties of biological tissues: a review.
# 2013 Phys. Med. Biol. 58 R37
# http://iopscience.iop.org/0031-9155/58/11/R37

readcsv <- function(filename)read.table(filename, header=TRUE, sep=",", comment.char="#")


