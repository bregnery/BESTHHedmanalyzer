# ED Analyzer for BEST using HH Monte Carlo

This ED Analyzer applieds the BEST standalone function to privately produced HH Monte Carlo samples.
BEST (Boosted Event Shape Tagger) outputs probabilities for the heavy particle corresponding to 
each AK8 jet in an event. For more information see [BESTAnalysis](https://github.com/justinrpilot/BESTAnalysis/tree/master).

## Installation 

This ED Analyzer has a couple of dependencies: CMSSW_9_X, BEST, and lwtnn.

### Dependencies

First, install CMSSW_9_X on CERN's lxplus or Fermilab's lpc.

```bash
cmsrel CMSSW_9_4_8
cd CMSSW_9_4_8/src/
scram b -j8
```

Then, install and compile [lwtnn](https://github.com/demarley/lwtnn/tree/CMSSW_8_0_X-compatible#cmssw-compatibility)

```bash
cd CMSSW_9_4_8/src/
mkdir lwtnn
cd lwtnn
git clone https://github.com/demarley/lwtnn.git
cd lwtnn
git checkout CMSSW_8_0_X-compatible
scram b -j8
```

Now, install, clean, and compile [BEST](https://github.com/justinrpilot/BESTAnalysis/tree/master)

```bash
cd CMSSW_9_4_8/src/
git clone https://github.com/justinrpilot/BESTAnalysis.git 
cd BESTAnalysis
git checkout master
rm -r -f BESTAnalyzer
rm -r -f BESTProducer
cd ..
scram b -j8
```

Finally, update the dnnFile path in ``BESTAnalysis/BoostedEventShapeTagger/data/config.txt`` to the full
working directory and compile.

```bash
cd CMSSW_9_4_8/src/
scram b -j8
```

### Install ED Analyzer

The ED Analyzer simply needs to be cloned and compiled.

```bash
cd CMSSSW_9_4_8/src/
mkdir BESTHHedmanalyzer
cd BESTHHedmanalyzer
git clone https://github.com/bregnery/BESTHHedmanalyzer.git
cd ..
scram b -j8
```

Now, the ED Analyzer should be ready to use!

## Run Instructions

The actual analyzer is located in the ``BESTHHedmanalyzer/BESTHHedmanalyzer/plugins/BESTHHedmanlyzer.cc`` and
the run instructions are located in ``BESTHHedmanalyzer/BESTHHedmanalyzer/test/run_HH.py``. To run, use the
cms environment to run the ``run_HH.py`` file. 

```bash
cd BESTHHedmanalyzer/BESTHHedmanalyzer/test/
cmsenv
cmsRun run_HH.py
```

Be sure to update any file locations in ``run_HH.py``!!

