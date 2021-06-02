# Data Assimilation Predictive GAN (DA-PredGAN) - Forecasting spatial variation of COVID-19 infection using GAN

- 1.Compress_train.ipynb -> Compress the 800 variables per time step to 15 variables per time step using PCA. 

- 2.Compress_test.ipynb -> Apply the PCA Compression fitted on the training dataset to the test dataset. 

- 3.GAN-training.ipynb -> Train a GAN and save the model 

- 4.GAN-Prediction.ipynb -> Predict with the GAN 

- 5.GAN-DataAssimilation.ipynb -> Assimilate observed data with the GAN 
