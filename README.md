# Basic tutorial to run the framework

### Prepare code and dataset
- Clone the github repository:

 ```
 git clone https://github.com/BiomedicalMachineLearning/VGP.git
 ```

- You need to download from this ggdrive link

data.zip: [https://drive.google.com/file/d/1HGFyQDnCKoGZEjXudEShVub3iS_DqpYd/view?usp=sharing](https://drive.google.com/file/d/1HGFyQDnCKoGZEjXudEShVub3iS_DqpYd/view?usp=sharing)

Unzip the data.zip file and put the folder `data` into the `tutorials` folder 

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
