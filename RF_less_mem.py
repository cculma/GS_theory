# random forest with less memory
# taken form https://mljar.com/blog/random-forest-memory/

import os
import joblib
import pandas as pd
import numpy as np
from sklearn.ensemble.forest import RandomForestClassifier
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import train_test_split
from sklearn.metrics import log_loss
from matplotlib import pyplot as plt
