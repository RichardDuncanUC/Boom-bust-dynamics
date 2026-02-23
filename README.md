# Boom-bust-dynamics
Data and code for analyis of fox and hare bounties as they spread across NSW.

The R script "Simulate travelling boom bust wave.R" simulates a consumer invading a landscape with resources at carrying capacity and draws Fig. 1.

The R scripts "CR and Gompertz models.R" and "CR model with predation.R" contain the Nimble code for fitting the consumer-resource and Gompertz models, including the hare consumer-resource model with fox predation.

The R scripts starting with "Fit...." read in the data and Nimble code, fit the models and write the output to the Model output directory. You need to run "Fit CR model.R" first for both hares and foxes. This fits the consumer-resource models - some of that output is used in fitting the other models.

Having fit all the models and stored the results in the model output directory, the script "Analyse fox and hare results.R" produces all of the figures and tables in the text, except for Figure 1. The model output is already stored in the Model output directory, so you can jsut run this script to reproduce the figures and table.

