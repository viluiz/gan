print('---------------------------------------------------------------')
print('Program developed by Vinicius L S Silva (vs1819) - 14.04.2021')
print('---------------------------------------------------------------')

# ==============================================================================
# Import modules
# ==============================================================================

# Python ≥3.5 is required
import sys
assert sys.version_info >= (3, 5)

# Scikit-Learn ≥0.20 is required
import sklearn
assert sklearn.__version__ >= "0.20"

# TensorFlow ≥2.0 is required
import tensorflow as tf
from tensorflow.keras.models import load_model
assert tf.__version__ >= "2.0"

# Common imports
import numpy as np
import os
import time
from collections import deque
import joblib

# ==============================================================================
# Input parameters
# ==============================================================================

# -----------------------------------------
# Load data
# -----------------------------------------

## Load PCA and Scaler
#
pca_compress = joblib.load("../../pca_compress_15.pkl") 
scaler = joblib.load("../../scaler.pkl")
# Number of POD coeffients (Global variable)
pca_size = pca_compress.n_components_
# PCA and scaler variables for obs_loss (Global variables)
data_range_ = tf.constant(scaler.data_range_, dtype=tf.float32)
data_min_ = tf.constant(scaler.data_min_, dtype=tf.float32)
components_ = tf.constant(pca_compress.components_, dtype=tf.float32)
pca_mean_ = tf.constant(pca_compress.mean_, dtype=tf.float32)

# Load GAN model (Global variables)
#
generator, discriminator = load_model("../../gan-tfex-DCGAN-5kernel.h5", compile=False).layers
# Size of the latent space (Global variable)
_, latent_size = generator.input_shape 
# Number of consecutive time steps and the number of POD + model parameters (Global variables)
_, ntimes, codings_size, _  = generator.output_shape 

# Load observed data and initial condition 
#
initial_coding, X_obs, obs_points, forward_steps = joblib.load("input_data.pkl")    
# (groups, lines, columns) Size of a time snapshot in the high dimension space (Global variables)
_, ng, nl, nc = X_obs.shape 

# -----------------------------------------
# Calculate the weights
# -----------------------------------------

# Delta variables 
#
delta_obs = 1200
delta_pod = 2
delta_r0 = 2

# Weights (Global variables)
#
# POD coeficients
weight_pod = pca_compress.singular_values_/pca_compress.singular_values_.sum()
weight_pod = tf.constant(weight_pod, dtype=tf.float32)
print('weight_pod: ', weight_pod, end='\n\n')
# R0s
weight_R0_p = 0.01*((delta_pod/delta_r0)**2)/2
weight_R0_p = tf.constant([weight_R0_p]*2, dtype=tf.float32)
print('weight_R0_p: ', weight_R0_p, end='\n\n')
weight_R0_da = 0.0001*((delta_pod/delta_r0)**2)/2
weight_R0_da = tf.constant([weight_R0_da]*2, dtype=tf.float32)
print('weight_R0_da: ', weight_R0_da, end='\n\n')
weight_R0 = tf.Variable(weight_R0_p, dtype=tf.float32)
# Observed data
weight_obs = 10.0*((delta_pod/delta_obs)**2)*ntimes
weight_obs = tf.constant(weight_obs, dtype=tf.float32)
print('weight_obs', weight_obs, end='\n\n')

# -----------------------------------------
# Set optimization tolerances
# -----------------------------------------

tol = 1e-6 # relative loss function difference 
tol_cv = 1e-4 # relative latent vector norm difference

# ==============================================================================
# Data assimilation functions 
# ==============================================================================

# POD loss
def pod_sse_loss(real_coding, gen_output):
    # -- POD coefficients --
    #
    # apply weights
    inp = tf.reshape(real_coding, [-1, codings_size])[:,:pca_size]
    out = tf.reshape(gen_output, [-1, codings_size])[:,:pca_size]
    
    pod_loss = tf.reduce_sum(tf.math.squared_difference(inp, out)*weight_pod)
    return pod_loss

# R0s loss
def r0_sse_loss(real_coding, gen_output):
    # -- R0s --
    #
    # apply weights
    inp = tf.reshape(real_coding, [-1, codings_size])[:,pca_size:]
    out = tf.reshape(gen_output, [-1, codings_size])[:,pca_size:]
    
    r0_loss = tf.reduce_sum(tf.math.squared_difference(inp, out)*weight_R0)
    return r0_loss

# Observed data loss
def obs_sse_loss(gen_output, obs_data, obs_points_inrange): 
    # -- Observed data --
    # 
    if obs_data.shape[0]>0:
        # scaler.inverse_transform
        gen_output_scaled = ((tf.reshape(gen_output, [-1, codings_size])+1)/2)*data_range_ + data_min_  
        # pca_compress.inverse_transform
        X_generated = tf.reshape(gen_output_scaled[:,:pca_size]@components_+pca_mean_, (len(gen_output_scaled), ng, nl, nc))
        # get the data points
        sim_data = tf.gather_nd(X_generated, obs_points_inrange)
        
        obs_loss = tf.reduce_sum(tf.math.squared_difference(obs_data, sim_data))*weight_obs/tf.size(obs_data, out_type=tf.float32) 
    else:
        obs_loss = 0
        
    return obs_loss

# One time iteration forward without observed data
def optimize_coding_forward_onlyPOD(latent_values, real_coding, obs_data, obs_points_inrange, epochs=1000):
    
    optimizer = tf.keras.optimizers.Adam(1e-2)
    
    @tf.function
    def opt_step_forward(optimizer, latent_values, real_coding, obs_data, obs_points_inrange):
        with tf.GradientTape() as tape:
            tape.watch(latent_values)
            gen_output = generator(latent_values, training=False) 

            pod_loss = pod_sse_loss(real_coding, gen_output[:,:(ntimes - 1),:,:])   
            r0_loss = r0_sse_loss(real_coding, gen_output[:,:(ntimes - 1),:,:])   
            obs_loss = obs_sse_loss(gen_output[:,:(ntimes - 1),:,:], obs_data, obs_points_inrange)
            loss = pod_loss + r0_loss 

        gradient = tape.gradient(loss, latent_values)  
        optimizer.apply_gradients(zip([gradient], [latent_values]))  

        return loss, pod_loss, r0_loss, obs_loss  
    
    loss = []
    loss.append(opt_step_forward(optimizer, latent_values, real_coding, obs_data, obs_points_inrange))
    latent_values_old = latent_values.numpy()
        
    for epoch in range(epochs):
        loss.append(opt_step_forward(optimizer, latent_values, real_coding, obs_data, obs_points_inrange))
        
        err = abs(loss[-1][0]-loss[-2][0])/loss[-2][0]
        err_cv = np.linalg.norm(latent_values.numpy() - latent_values_old)/np.linalg.norm(latent_values_old)
        if err < tol and err_cv < tol_cv:
            break
        latent_values_old = latent_values.numpy()          
        
    return latent_values, loss[-1]  #returns the optimized input that generates the desired output

# One time iteration forward
def optimize_coding_forward(latent_values, real_coding, obs_data, obs_points_inrange, epochs=1000):
    
    optimizer = tf.keras.optimizers.Adam(1e-2)
    
    @tf.function
    def opt_step_forward(optimizer, latent_values, real_coding, obs_data, obs_points_inrange):
        with tf.GradientTape() as tape:
            tape.watch(latent_values)
            gen_output = generator(latent_values, training=False) 

            pod_loss = pod_sse_loss(real_coding, gen_output[:,:(ntimes - 1),:,:])   
            r0_loss = r0_sse_loss(real_coding, gen_output[:,:(ntimes - 1),:,:])   
            obs_loss = obs_sse_loss(gen_output[:,:(ntimes - 1),:,:], obs_data, obs_points_inrange)
            loss = pod_loss + r0_loss + obs_loss 

        gradient = tape.gradient(loss, latent_values)  
        optimizer.apply_gradients(zip([gradient], [latent_values]))  

        return loss, pod_loss, r0_loss, obs_loss
    
    loss = []
    loss.append(opt_step_forward(optimizer, latent_values, real_coding, obs_data, obs_points_inrange))
    latent_values_old = latent_values.numpy()
        
    for epoch in range(epochs):
        loss.append(opt_step_forward(optimizer, latent_values, real_coding, obs_data, obs_points_inrange))
        
        err = abs(loss[-1][0]-loss[-2][0])/loss[-2][0]
        err_cv = np.linalg.norm(latent_values.numpy() - latent_values_old)/np.linalg.norm(latent_values_old)
        if err < tol and err_cv < tol_cv:
            break
        latent_values_old = latent_values.numpy()          
       
    return latent_values, loss[-1]  #returns the optimized input that generates the desired output

# One time iteration backward
def optimize_coding_backward(latent_values, real_coding, obs_data, obs_points_inrange, epochs=1000):
    
    optimizer = tf.keras.optimizers.Adam(1e-2)
    
    @tf.function
    def opt_step_backward(optimizer, latent_values, real_coding, obs_data, obs_points_inrange):
        with tf.GradientTape() as tape:
            tape.watch(latent_values)
            gen_output = generator(latent_values, training=False) 

            pod_loss = pod_sse_loss(real_coding, gen_output[:,1:,:,:])   
            r0_loss = r0_sse_loss(real_coding, gen_output[:,1:,:,:])   
            obs_loss = obs_sse_loss(gen_output[:,1:,:,:], obs_data, obs_points_inrange)
            loss = pod_loss + r0_loss + obs_loss 

        gradient = tape.gradient(loss, latent_values)  
        optimizer.apply_gradients(zip([gradient], [latent_values]))  

        return loss, pod_loss, r0_loss, obs_loss 
    
    loss = []
    loss.append(opt_step_backward(optimizer, latent_values, real_coding, obs_data, obs_points_inrange))
    latent_values_old = latent_values.numpy()
        
    for epoch in range(epochs):
        loss.append(opt_step_backward(optimizer, latent_values, real_coding, obs_data, obs_points_inrange))
        
        err = abs(loss[-1][0]-loss[-2][0])/loss[-2][0]
        err_cv = np.linalg.norm(latent_values.numpy() - latent_values_old)/np.linalg.norm(latent_values_old)
        if err < tol and err_cv < tol_cv:
            break
        latent_values_old = latent_values.numpy()          
        
    return latent_values, loss[-1]  #returns the optimized input that generates the desired output

# Process observed data
def process_obs_data(march_range, X_obs, obs_points):
    # Process the observed data to be in the range of the forward/backward march
    
    X_obs_inrange = X_obs[march_range]
    obs_points_inrange = np.array([x-np.array([march_range[0], 0, 0, 0]) for x in obs_points if x[0] in march_range])
    if obs_points_inrange.size>0:
        obs_data = tf.gather_nd(X_obs_inrange, obs_points_inrange)
    else:
        obs_data = []

    obs_points_inrange = tf.constant(obs_points_inrange)
    obs_points_inrange = tf.cast(obs_points_inrange, dtype=tf.int32)
    obs_data = tf.constant(obs_data)
    obs_data = tf.cast(obs_data, dtype=tf.float32)
    
    return obs_data, obs_points_inrange
       
# Forward and backward marches    
def forward_march(initial_coding, X_obs, obs_points, forward_steps):
    start_pred = time.time()
     
    #-----------------------------------------
    # Forward march
    #-----------------------------------------
    print('-- Time Stepping Forward March --')

    # Initialize variables 
    weight_R0.assign(weight_R0_p) 
    march_range = np.arange(0,(ntimes-1))
    obs_data, obs_points_inrange = process_obs_data(march_range, X_obs, obs_points)
    loss_iteration = []
    forward_loss = []
 
    real_coding = initial_coding
    real_coding = real_coding[:-1,:]
    R0s_run = real_coding[0,-2:]
    real_coding = tf.constant(real_coding)
    real_coding = tf.cast(real_coding, dtype=tf.float32) 

    # Find first latent values
    list_latent_values = []
    list_loss = []
    for k in range(30):
        latent_values = tf.random.normal([1, latent_size], mean=0.0, stddev=0.01*10**(k//10))  
        latent_values = tf.Variable(latent_values) 
        latent_values, loss = optimize_coding_forward_onlyPOD(latent_values, real_coding, obs_data, obs_points_inrange, epochs=1000)
        list_latent_values.append(latent_values)
        list_loss.append(loss)
    latent_values_onlyPOD = list_latent_values[np.argmin(list_loss,axis=0)[0]]
    loss_iteration.append(list(map(float,list_loss[np.argmin(list_loss,axis=0)[0]])))
    print('Initial point losses: ', end='\n')
    for k, l in enumerate(list_loss):
        if k == np.argmin(list_loss,axis=0)[0]:
            print('{:g}: *'.format(k), ['{0:1.2e}'.format(float(x)) for x in l])
        else:
            print('{:g}: '.format(k), ['{0:1.2e}'.format(float(x)) for x in l])
    print()
    print('Loss iteration 0: '+str(['{0:1.2e}'.format(float(x)) for x in list_loss[np.argmin(list_loss,axis=0)[0]]]), end='\n')      

    X_predict_forward_last = list(generator(latent_values_onlyPOD).numpy().reshape(-1,codings_size))
    gen_predict = X_predict_forward_last[-1]
    gen_predict[-2:] = R0s_run

    real_coding = np.concatenate((real_coding, gen_predict.reshape(1,-1)), axis=0)[1:,:]
    real_coding = tf.constant(real_coding)
    real_coding = tf.cast(real_coding, dtype=tf.float32)

    # Time stepping 
    for i in range(1, forward_steps): 
        march_range = np.arange(i,(ntimes-1)+i)
        obs_data, obs_points_inrange = process_obs_data(march_range, X_obs, obs_points)

        start = time.time()
        latent_values_onlyPOD, loss = optimize_coding_forward_onlyPOD(latent_values_onlyPOD, real_coding, obs_data, obs_points_inrange, epochs=1000)
        loss_iteration.append(list(map(float,loss)))
        print('Loss iteration '+str(i)+': '+str(['{0:1.2e}'.format(float(x)) for x in loss]), end=' - ')

        gen_predict = generator(latent_values_onlyPOD).numpy().reshape(-1,codings_size)[-1]
        X_predict_forward_last.append(gen_predict.copy())
        gen_predict[-2:] = R0s_run

        real_coding = np.concatenate((real_coding, gen_predict.reshape(1,-1)), axis=0)[1:,:]
        real_coding = tf.constant(real_coding)
        real_coding = tf.cast(real_coding, dtype=tf.float32)
        print ('{:.0f}s'.format( time.time()-start))
    X_predict_forward_last = np.array(X_predict_forward_last)
    print('Loss iteration mean:', ['{0:1.2e}'.format(float(x)) for x in np.array(loss_iteration).mean(axis=0)])
    
    print ('\nRuntime: {:.0f}s'.format(time.time()-start_pred))
    print ('Iterations: {:g}'.format(0))
    print ('Final obs loss: {:.2e}\n'.format(np.array(loss_iteration).mean(axis=0)[-1]))
    
    forward_loss.append(np.array(loss_iteration).mean(axis=0))
    
    return X_predict_forward_last, forward_loss
# ==============================================================================
# Run case   
# ==============================================================================
np.random.seed(0)
tf.random.set_seed(0)

output_data = forward_march(initial_coding, X_obs, obs_points, forward_steps)
joblib.dump(output_data, "output_data.pkl") 

# ==============================================================================
# End  
# ==============================================================================

















