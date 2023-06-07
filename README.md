
# Data and code required to reproduce all analyses and figures in
# “Quantifying cellular dynamics in mice using a novel fluorescent division reporter system”
# Eva Lukas, Thea Hogan,Cayman Williams, Benedict Seddon, Andrew J. Yates


## Sequence of events:

- Run data/Data Processing.Rmd
Generates the data frames for all fitting:
"CD4_CD8","demog", "final_data_simple_model","final_data_flow_model","final_data_flow_model_v2"

- Run fitting/Fitting A - Thymic Functions Y and Z.Rmd to get the empirical functions for Y(t) and Z(t)
	The fitting for Z(t) can take a few manual tries to converge.
- Make Figure S.4 with this (.Rmd files for making figures are within the scripts_figures folder); plots Y(t) and Z(t), with data 


- Run fitting/Fitting B - Simple model no Ki67.Rmd in the fitting folder  
- Make Figure 4 (simple model fits, on counts and fractions). Includes mouse age as point size. 


- Run Fitting C to make the function K+(t) 
	… this makes parameters/p_Ki67_high_YFP.RData that are used together with the thymic function Y(t) from above to fit the Ki67 model (Eqn. 5)
- Make Figure S5, showing K+(t) 


- Run Fitting D1, D2, D3 - Ki67 (“flow”) model with fixed thymic inputs Y(t) and K+(t) (eqn 5).
	D1 - gets point estimates for the models with p free, and p set to zero
	D2 - bootstraps for p free
	D3 - bootstraps for p zero

- Make Figure 7 (flow model with fixed inputs). Also includes mouse age as point size.

- Fitting E - Simple model, estimating thymus simultaneously 

- Fitting F - Ki67 (flow) model, estimating thymus simultaneously
	F1 - gets point estimates for the models with p free, and p set to zero
	F2 - bootstraps for p free
	F3 - bootstraps for p zero

- Fitting G - F tests, comparing improvement in fit by adding p in both Ki67 models.

- Make Figure 8

- Make Figure 5 (summary of residence time estimates from all models)

- Make Figure 9 and Fig 9 alt (= Fig S6) - bulk Ki67 predictions.
