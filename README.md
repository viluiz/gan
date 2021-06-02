# GAN for time series prediction, data assimilation and uncertainty quantification

This repository is the official implementation of: 

[Digital twins based on bidirectional LSTM and GAN for modelling the COVID-19 pandemic](https://arxiv.org/abs/2102.02664). 

[Data Assimilation Predictive GAN (DA-PredGAN): applied to determine the spread of COVID-19](https://arxiv.org/abs/2105.07729). 

[GAN for time series prediction, data assimilation and uncertainty quantification](https://arxiv.org/abs/2105.13859). 

## Directories:

- **PredGAN**: Prediction using GAN - applied to the spatio-temporal spread of COVID-19 in an idealized town.
- **DA-PredGAN**: Prediction and data assimilation using GAN - applied to the spatio-temporal spread of COVID-19 in an idealized town.
- **UQ-PredGAN**: Prediction, data assimilation and uncertainty quantification using GAN - applied to the spatio-temporal spread of COVID-19 in an idealized town.

## Requirements

To install requirements:

```setup
 $ conda env create -f environment.yml 
 $ conda activate py3ml
 $ python -m ipykernel install --user --name=python3
```

Finally, start Jupyter:

```start
 $ jupyter notebook
```
