## Uncertainty Quantification Predictive GAN (UQ-PredGAN) - Forecasting spatial variation of COVID-19 infection using GAN

- 1.ReadCompress.ipynb -> Read and compress the training snapshots (time steps) using PCA. Apply the PCA compression to the test datasets. 

- 2.GAN-training.ipynb -> Train a GAN and save the model

- 3.GAN-Prediction.ipynb -> Predict with the GAN 

- 4.GAN-DataAssimilation.ipynb -> Assimilate observed data with the GAN 

- 5.GAN-UncertaintyQuantification.ipynb -> Quantify uncertainty with the GAN (with regularization) 

### First execute in order:
 
1.ReadCompress.ipynb 

2.GAN-training.ipynb (optional - the trained model and scaler are already on the folder *.h5 / *.pkl)

To execute: 

From inside the notebooks 
```
Cell->Run All 
```

From the command line:
```
$ jupyter nbconvert --to notebook --execute <notebookname>.ipynb
```

### Run an example of prediction 

3.GAN-Prediction.ipynb (optional)

### Run an example of data assimilation 

4.GAN-DataAssimilation.ipynb (optional)

### Run the uncertainty quantification 

5.GAN-UncertaintyQuantification.ipynb

## Requirements

To install requirements:

```setup
 $ conda env create -f environment.yml 
 $ conda activate py3ml
 $ python -m ipykernel install --user --name=python3 (optional)
```

Finally, start Jupyter:

```start
 $ jupyter notebook
```


