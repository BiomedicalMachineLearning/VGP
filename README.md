# Scalable PRSUP framework to run the transfer PRS from large population to under-representative population

<p align="center">
  <img src="https://i.ibb.co/brvQX90/run-pipeline-1.png"
    alt="deepreg_logo" title="DeepReg" width="600"/>
</p>


### Prepare dataset

You need to download from this ggdrive link:
- data.zip: [https://drive.google.com/file/d/1EVVIpVa_9sr_S8Raq_54kJ9Y67dxLdde/view?usp=share_link](https://drive.google.com/file/d/1EVVIpVa_9sr_S8Raq_54kJ9Y67dxLdde/view?usp=share_link)

Unzip the data.zip file and put into the `tutorials` folder 

### Install the requirements

I recommend to use conda to setup the environment

```
git clone https://github.com/BiomedicalMachineLearning/VGP.git
cd VGP
conda env create -f environment.yml
```

Besides, it's required to install R packages for LDpred2 and PolyFun/PolyPred. In R console:

```
install.packages(c("bigsnpr","susieR"))
```

### Tutorials

- Multiple populations for transfer PRS: [Link](https://github.com/BiomedicalMachineLearning/VGP/blob/main/Multiple_population_tutorial.ipynb)
