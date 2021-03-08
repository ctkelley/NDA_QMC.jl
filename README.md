# NDA_QMC.jl
# Open the notebook __NDA_Notebook.ipynb__ to see what to do.

This is the repository for the SN + NDA + QMC + GMRES stuff

This repo is part of the [CEMeNT](https://cement-psaap.github.io) project

If you're not a CEMeNTer, go away.

## Contents

- [Instructions](#What-to-do)
- [Dependencies](#Dependencies) get these packages first
- [Problems?](#Problems)

## What to do.

- Clone this repo
- Install the dependencies
- Fire up __IJulia__ and cd to the directory where you put the repo
- Open the notebook __NDA_Notebook.ipynb__
- Do a ```run all``` this will run __nda_note.jl__ which gets the dependencies organized and puts the right things in your load path
- You are now ready to mess with the sovlers. The codes in the /src directory
  are pretty well documented and ready for you to play with.
  
## The transport solves

I'm duplicting the results from

R.D.M. Garcia and C.E. Siewert.
 Radiative transfer in finite inhomogeneous plane-parallel
  atmospheres.
 <em>J. Quant. Spectrosc. Radiat. Transfer</em>, 27:141-148, 1982.

I run two classes of sovlers so far

- Plain-vanilla source iteration and GMRES accelerated source iteration
- NDA with Picard iteration and NDA with Newton-GMRES
- QMC sweep functions

QMC solution iteration is next on the list.

## Dependencies:
- [SIAMFANLEquations.jl](https://github.com/ctkelley/SIAMFANLEquations.jl)
  - You need v 0.3.0 or higher. Get it from the pkg REPL with 
  
  ```add SIAMFANLEquations```

- FastGaussQuadrature

- LaTeXStrings

- PyPlot

- SuiteSparse

- SpecialFunctions

- LinearAlgebra

- Sobol

- JLD2

- Random


## Problems
If you're CEMeNTer, you know how to find me. Otherwise ... go away.
