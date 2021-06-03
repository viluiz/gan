## Data Assimilation Predictive GAN (DA-PredGAN) - Forecasting spatial variation of COVID-19 infection using GAN

- 1.Compress_train.ipynb -> Compress the training snapshots (time steps) using PCA. 

- 2.Compress_test.ipynb -> Apply the PCA Compression to the first test dataset. 

- 3.Compress_test-uq.ipynb -> Apply the PCA Compression to the second test dataset. 

- 4.GAN-training.ipynb -> Train a GAN and save the model

- 5.GAN-Prediction.ipynb -> Predict with the GAN 

- 6.GAN-DataAssimilation.ipynb -> Assimilate observed data with the GAN 

- 7.GAN-UncertaintyQuantification.ipynb -> Quantify uncertainty with the GAN 

### Requirements

First execute in order:
 
1.Compress_train.ipynb 

2.Compress_test.ipynb

3.Compress_test-uq.ipynb

4.GAN-training.ipynb (optional - the trained model and scaler are already on the folder *.h5 / *.pkl)

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

5.GAN-Prediction.ipynb (optional)

### Run an example of data assimilation 

6.GAN-DataAssimilation.ipynb (optional)

### Run the uncertainty quantification 

7.GAN-UncertaintyQuantification.ipynb


