# A GAN-based reduced order model for prediction, data Assimilation and uncertainty quantification

This repository is the official implementation of: 

[Digital twins based on bidirectional LSTM and GAN for modelling the COVID-19 pandemic](https://www.sciencedirect.com/science/article/pii/S0925231221015290) (for the Predictive GAN).

[Data Assimilation Predictive GAN (DA-PredGAN): applied to determine the spread of COVID-19](https://arxiv.org/abs/2105.07729). 

[A GAN-based Reduced Order Model for Prediction, Data Assimilation and Uncertainty Quantification](https://arxiv.org/abs/2105.13859). 

## Directories:

- **PredGAN**: Prediction using GAN - applied to the spatio-temporal spread of COVID-19 in an idealized town.
- **DA-PredGAN**: Data assimilation using GAN - applied to the spatio-temporal spread of COVID-19 in an idealized town.
- **UQ-PredGAN**: Uncertainty quantification using GAN - applied to the spatio-temporal spread of COVID-19 in an idealized town.
- **datasets**: Datasets of the spatio-temporal spread of COVID-19 in an idealized town. 

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

