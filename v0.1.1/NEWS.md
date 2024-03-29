# sfaR 0.1.1
Changes in 'sfaR' version 0.1.1 (2022-05-05)

## NEW FEATURES

* The package is enriched with a `NEWS` file

## BUG FIXES

* The `fitted()` function now works properly for objects of class `sfacross`

* When `extraPar = TRUE`, the `coef()` function now works properly on objects of class `lcmcross` for models with a number of classes (argument `lcmClasses`) equal to 4 or more.

## DEPRECATED & DEFUNCT
None

## OTHER USER-VISIBLE CHANGES

* The `efficiencies()` documentation is enriched with an `lcmcross` example.

* Updated `DESCRIPTION` file.

* Updated startup message.

***
# sfaR 0.1.0
Changes in 'sfaR' version 0.1.0 (released on CRAN 2021-05-04)

## NEW FEATURES

* When the function `efficiencies()` is applied to an object of class `lcmcross` an additional value is returned. This new value, `ineff_c#`, returns the conditional inefficiency for observations classified under class# only. Therefore a value is assigned to `ineff_c#` only if the observations is classified under class# (i.e. `Group_c = #`).

## BUG FIXES

* In `efficiencies()` function applied to an object of class `lcmcross`, the values for `u_c`, `teJLMS_c`, `Group_c`, and `PosteriorProb_c` are now returned correctly. In the previous version (v 0.0.91) all observations were applied the same values when the argument `lcmClasses` was set to 3 and more. This bug is now fixed.

## DEPRECATED & DEFUNCT
None

## OTHER USER-VISIBLE CHANGES

* Updated `DESCRIPTION` file.

* In the `efficiencies()` function applied to an object of class `lcmcross`, the returned elements have been reordered in order to better reflect the logic of the `lcmcross` model.

* Updated `efficiencies()` documentation to incorporate the new `ineff_c#` elements returned by the function, as well as the reordering of the returned values.

* Updated startup message.

***
# sfaR 0.0.91
First public release (released on CRAN 2021-03-05)