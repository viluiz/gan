## (outdated see UQ-PredGAN folder) Data Assimilation Predictive GAN (DA-PredGAN) - Forecasting spatial variation of COVID-19 infection using GAN

- 1.Compress_train.ipynb -> Compress the training snapshots (time steps) using PCA. 

- 2.Compress_test.ipynb -> Apply the PCA Compression to the test dataset. 

- 3.GAN-training.ipynb -> Train a GAN and save the model

- 4.GAN-Prediction.ipynb -> Predict with the GAN 

- 5.GAN-DataAssimilation.ipynb -> Assimilate observed data with the GAN 

### First execute in order:

1.Compress_train.ipynb 

2.Compress_test.ipynb

3.GAN-training.ipynb (optional - the trained model and scaler are already on the folder *.h5 / *.pkl)

To execute: 

From inside the notebooks 
```
Cell->Run All 
```

From the command line
```
$ jupyter nbconvert --to notebook --execute <notebookname>.ipynb
```

### Run the prediction 

4.GAN-Prediction.ipynb 

### Run the data assimilation 

5.GAN-DataAssimilation.ipynb

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


