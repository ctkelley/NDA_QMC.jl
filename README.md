# NDA_QMC.jl
# Open the notebook __NDA.ipynb__ to see what to do.

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
- Open the notebook __NDA.ipynb__
- Do a ```run all``` this will run __nda_note.jl__ which gets the dependencies organized and puts the right things in your load path
- You are now ready to mess with the sovlers. Right now all that I have in the notebook is __compare.jl__. The other codes in the /src directory
  are pretty well documented and ready for you to play with.

## Dependencies:
- [SIAMFANLEquations.jl](https://github.com/ctkelley/SIAMFANLEquations.jl)
  - You need v 0.3.0 or higher. As of today, you have to get it from the pkg REPL with 
  
  ```add https://github.com/ctkelley/SIAMFANLEquations.jl#master```
  
     because I have not pushed 0.3.0 yet.

- FastGaussQuadrature

- LaTeXStrings

- PyPlot

## Problems
If you're Ryan or one of this students, you know how to find me. Otherwise ... go away.
