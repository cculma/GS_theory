# random forest with less memory
# taken form https://mljar.com/blog/random-forest-memory/

import os
import joblib
import pandas as pd
import numpy as np
# from sklearn.ensemble.forest import RandomForestClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import train_test_split
from sklearn.metrics import log_loss
from matplotlib import pyplot as plt
from sklearn.model_selection import cross_val_score


df = pd.read_csv("https://raw.githubusercontent.com/pplonski/datasets-for-start/master/adult/data.csv", 
                 skipinitialspace=True)
# df.info()

df = df.fillna(df.mode().iloc[0])
for col in df.columns:
    if df[col].dtype == "object":
        encode = LabelEncoder()
        df[col] = encode.fit_transform(df[col])

X = df[df.columns[:-1]]
# Z = df.drop('age', 1) # junst in case to delete column one 'age'
# print(Z)
y = df["income"]

## Let’s use 25% of the data for testing and the rest for training.
# X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.25, random_state=56)

## Let’s use 10% of the data for testing and the rest for training.
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.1, random_state=56)

rf = RandomForestClassifier()
rf.fit(X_train, y_train)
# print(rf.estimators_[0].tree_.max_depth)

depths = [tree.tree_.max_depth for tree in rf.estimators_]
# print(f"Mean tree depth in the Random Forest: {np.round(np.mean(depths))}")

joblib.dump(rf.estimators_[0], "first_tree_from_RF.joblib") 
print(f"Single tree size: {np.round(os.path.getsize('first_tree_from_RF.joblib') / 1024 / 1024, 2) } MB")

joblib.dump(rf, "RandomForest_100_trees.joblib") 
print(f"Random Forest size: {np.round(os.path.getsize('RandomForest_100_trees.joblib') / 1024 / 1024, 2) } MB")

print(f"Cross validation: {np.mean(cross_val_score(estimator, X))}")
###########################
# 10-Fold Cross validation
print np.mean(cross_val_score(rf, X_train, y_train, cv=10))

###########################

# y_predicted = rf.predict_proba(X_test)
# rf_loss = log_loss(y_test, y_predicted)
# print(rf_loss)

# shallow_rf = RandomForestClassifier(max_depth=6)
# shallow_rf.fit(X_train, y_train)

# joblib.dump(shallow_rf.estimators_[0], "first_tree_from_shallow_RF.joblib") 
# # print(f"Single tree size from shallow RF: {np.round(os.path.getsize('first_tree_from_shallow_RF.joblib') / 1024 / 1024, 2) } MB") # this line is the same as 38

# joblib.dump(shallow_rf, "Shallow_RandomForest_100_trees.joblib") 
# print(f"Shallow Random Forest size: {np.round(os.path.getsize('Shallow_RandomForest_100_trees.joblib') / 1024 / 1024, 2) } MB")

# y_predicted = shallow_rf.predict_proba(X_test)
# shallow_rf_loss = log_loss(y_test, y_predicted)
# print(shallow_rf_loss)

# joblib.dump(rf, "RF_compressed.joblib", compress=3)
# # print(f"Compressed Random Forest: {np.round(os.path.getsize('RF_compressed.joblib') / 1024 / 1024, 2) } MB")
