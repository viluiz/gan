## Predictive GAN (PredGAN) - Forecasting spatial variation of COVID-19 infection using GAN

- 1.Compression.ipynb -> Compress the 800 variables per time step to 15 variables per time step (The compression I am using is from PCA. However, there are in this file other types of compression using autoencoders). 

- 2.GAN-training.ipynb -> Train a GAN and save the model 

- 3.GAN-Prediction.ipynb -> Predict with the GAN 

To execute: 

From inside the notebooks 
```
Cell->Run All 
```

From the command line
```
$ jupyter nbconvert --to notebook --execute <notebookname>.ipynb
