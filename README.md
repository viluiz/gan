# A GAN-based reduced order model for prediction, data Assimilation and uncertainty quantification

This repository is the official implementation of: 

[Digital twins based on bidirectional LSTM and GAN for modelling the COVID-19 pandemic](https://www.sciencedirect.com/science/article/pii/S0925231221015290) (for the Predictive GAN).

[Data Assimilation Predictive GAN (DA-PredGAN) Applied to a Spatio-Temporal Compartmental Model in Epidemiology](https://link.springer.com/article/10.1007/s10915-022-02078-1). 

[A GAN-based Reduced Order Model for Prediction, Data Assimilation and Uncertainty Quantification](https://arxiv.org/abs/2105.13859). 

## Directories:

- **PredGAN**(outdated see UQ-PredGAN folder): Prediction using GAN - applied to the spatio-temporal spread of COVID-19 in an idealized town.
- **DA-PredGAN**: Data assimilation using GAN - applied to the spatio-temporal spread of COVID-19 in an idealized town.
- **UQ-PredGAN**: Uncertainty quantification using GAN - applied to the spatio-temporal spread of COVID-19 in an idealized town.
- **datasets**: Datasets of the spatio-temporal spread of COVID-19 in an idealized town. 
- **GAN_evaluation**: New way of evaluating the GAN training. 
- **Regularization**: Regularization to improve the GAN-based Reduced Order Model. 
- **MCMC**: Comparison between the UQ-PredGAN and the Markov chain Monte Carlo (MCMC) methods.

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

