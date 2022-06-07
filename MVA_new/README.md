# WWG Analysis
Based on NanoAOD Tools: <https://github.com/cms-nanoAOD/nanoAOD-tools>

Dedicated for WWG analysis on lxplus or cms-connect environment

Introduction slides:
<br>
<https://docs.google.com/presentation/d/1Td-nFotttG4cUfyNEeJJGHtRFbVfc7ZF6uo8iMnB8z0/edit?usp=sharing>

--------------
## content

- [Download and setup](#Download-and-setup)
- [Generate Signal Sample](#Generate-Signal-Sample)
- [Local mode](#Baseline-selection)
- [Crab mode](#Crab-mode)
- [Condor mode](#Condor-mode)
- [Recent updates](#Recent-updates)


--------------
<br>

## <span id="Download-and-setup"> Download and setup </span> 

```bash
cmsrel CMSSW_10_6_20
cd CMSSW_10_6_20/src
cmsenv
git clone https://github.com/cms-nanoAOD/nanoAOD-tools.git PhysicsTools/NanoAODTools
cd PhysicsTools/NanoAODTools
git clone https://github.com/AnYpku/WWG.git
scram b
```
<br>
<br>

## <span id="Generate-Signal-Sample"> Generate Signal Sample </span> 

This will generate signal samples for WWA. you can generate 3 different schemes signal.you can    
```bash
voms-proxy-init -voms cms -valid 192:00

python submit_job.py
or
condor_submit submit_*.jdl
```
<br>

## <span id="Baseline-selection"> Local mode </span>
In WWG_seletor, `WWG_Module.py` is designed for selections. Use `WWG_postproc.py` to test.

```bash
python WWG_postproc.py -f *.root -y 2018 -d (for data)
python WWG_postproc.py -f *.root -y 2018 (for MC)
```

<br>


## <span id="Crab-mode"> Crab mode </span>
for Crab job, you can use creat_cfg.py to generate the crab scripts,and you can choose your datasets by changing dataset_2018_mc_nano_v7.py or dataset_2018_data_nano_v7.py

```bash
voms-proxy-init -voms cms -valid 192:00
python3 create_cfg.py -y 2018 -u 1 -m crab -k MC
crab submit -c cfg2018_mc/DY_cfg.py
```

<br>


## <span id="Condor-mode"> Condor mode </span>
In condor folder, `condor_for_postproc.py` is designed for preparing codes and submitting them to HTcondor. The purpose is to run over samples on DAS in parallel. Similarly, it has `-f` and `-n` arguments.

First you need to setup grid certification
```bash
voms-proxy-init -voms cms -valid 192:00
```
Modify `Proxy_path` in `condor_for_post.py` according to you own settings

```bash
python condor_for_postproc.py -y 2018 -n DY -k MC
```
## <span id="Recent-updates"> Recent Updates </span>
### HLT details

In this analysis, we use mixed HLT.  

you can change the HLT and datasets in  `WWG_Module.py` and `WWG_postproc.py`.

If you use crab mode, you have to modify `cfreate_cfg.py`.

### TMVA for some BKG

Some test versions for BKG like ZGJ are on going.

```bash
root -l  TMVAClassification.C
```




