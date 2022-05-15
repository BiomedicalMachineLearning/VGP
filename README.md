# Scalable TransPRS framework to run the trasfer PRS from large population to under-representative population

<p align="center">
  <img src="https://i.ibb.co/brvQX90/run-pipeline-1.png"
    alt="deepreg_logo" title="DeepReg" width="600"/>
</p>


### Prepare dataset

You need to download from this ggdrive link:
- data.zip: [https://drive.google.com/file/d/1ZsbfabvwnssMiji6ECDieg1C7mJNhvYm/view?usp=sharing](https://drive.google.com/file/d/1ZsbfabvwnssMiji6ECDieg1C7mJNhvYm/view?usp=sharing)

Unzip the data.zip file and put into the `tutorials` folder 

### Install the requirements

I recommend to use conda to setup the environment

```
conda env create -f environment.yml
```

Besides, it's required to install R packages for LDpred2 and PolyFun/PolyPred. In R console:

```
install.packages(c("bigsnpr","susieR"))
```

### Tutorials

- Multiple populations for transfer PRS: [Link](https://github.com/BiomedicalMachineLearning/VGP/blob/main/Multiple_population_tutorial.ipynb)
