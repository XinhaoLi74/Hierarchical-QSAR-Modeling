# Hierarchical-QSAR-Modeling
Code for reproducing the work in paper: **Hierarchical H-QSAR Modeling Approach for Integrating Binary/Multi Classification and Regression Models of Acute Oral Systemic Toxicity**

Some results and trained models are too large to be able to github. I will figure out a way to share the results and trained models.

## Data:
The curated training and test data are avaiable in [train_test_sets](https://github.com/XinhaoLi74/Hierarchical-QSAR-Modeling/tree/master/data/train_test_sets) folder.

## Code:
1. Prepare `labels` for modeling: [labels.ipynb](https://github.com/XinhaoLi74/Hierarchical-QSAR-Modeling/blob/master/notebooks/labels.ipynb).
### Base models
2. Compute chemical descriptors/fingerprints for base models: [descriptors.ipynb](https://github.com/XinhaoLi74/Hierarchical-QSAR-Modeling/blob/master/notebooks/descriptors.ipynb).
3. Descriptor Selection: [descriptors_selections.ipynb](https://github.com/XinhaoLi74/Hierarchical-QSAR-Modeling/blob/master/notebooks/descriptors_selections.ipynb).
4. Hyperparameter tuning of base models: [Base_models_selection.ipynb](https://github.com/XinhaoLi74/Hierarchical-QSAR-Modeling/blob/master/notebooks/Base_models_selection.ipynb).
5. Building base models with optimal hyperparameters: [Base_models.ipynb](https://github.com/XinhaoLi74/Hierarchical-QSAR-Modeling/blob/master/notebooks/Base_models.ipynb).
### Hierarchial Models
6. Meta features: [Hierarchical_features.ipynb](https://github.com/XinhaoLi74/Hierarchical-QSAR-Modeling/blob/master/notebooks/Hierarchical_features.ipynb).
7. Hyperparameter tuning of hierarchical models: [Hierarchical_models_selection.ipynb](https://github.com/XinhaoLi74/Hierarchical-QSAR-Modeling/blob/master/notebooks/Hierarchical_models_selection.ipynb).
8. Build hierarchial models with optimal hyperparameters: [Hierarchical_models.ipynb](https://github.com/XinhaoLi74/Hierarchical-QSAR-Modeling/blob/master/notebooks/Hierarchical_models.ipynb).
### Model Evaluation
9. Evaluate cross-validation and test set performance: [Model_evaluation.ipynb](https://github.com/XinhaoLi74/Hierarchical-QSAR-Modeling/blob/master/notebooks/Model_evaluation.ipynb).

## GUI
![](images/GUI.gif)
### Installation 

`
conda env create -f hqsar_env.yml
pip install streamlit
`

Run the GUI:

`
streamlit run GUI.py
`