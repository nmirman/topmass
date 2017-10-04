# TopMass

This project determines the top quark mass from data collected by the [CMS experiment](https://cms.cern) at CERN.  The resulting paper can be viewed [here](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.96.032002).

### How to run:

```
./DoFit <options>
```

Fit results and diagnostics plots are output to `results\`.  Post-analysis scripts are located in `scripts\`.  Batch jobs can be submitted to condor or crab using the `condorJob` and `crab.cfg` configuration files, respectively.

### Features:

- Maximum likelihood fitting framework, with numerical minimization conducted by [MINUIT2](http://seal.web.cern.ch/seal/snapshot/work-packages/mathlibs/minuit/).
- Calculation of physics observables MT2, Mbl, and MAOS for mass estimation.
- Gaussian process regression framework.
- Code interfaced with the [ROOT](https://root.cern.ch) physics and visualization libraries.
- Support for batch execution on local clusters and the [LHC Computing GRID](http://wlcg.web.cern.ch).
