import pandas as pd
import seaborn as sb

import pingouin as pg

# Python equivalent of `ggplot2`
from plotnine import *

# Statistical models, conducting tests and statistical data exploration
import statsmodels.api as sm

# Convenience interface for specifying models using formula strings and DataFrames
import statsmodels.formula.api as smf

USArrests_py = pd.read_csv("/home/participant/Course_Materials/data/CS3-usarrests.csv")

# create scatterplot of the data
myplot = (ggplot(USArrests_py,
         aes(x = "assault",
             y = "murder")) +
     geom_point())

# create a linear model
model = smf.ols(formula= "murder ~ assault", data = USArrests_py)
# and get the fitted parameters of the model
lm_USArrests_py = model.fit()

#print(USArrests_py)

plt = (ggplot(USArrests_py,
        aes(x = "assault", y = "murder")) +
     geom_point() +
     geom_smooth(method = "lm",
                 se = False,
                 colour = "blue"))
#print(plt)

print(USArrests_py.head())
USArrests_cor_py = USArrests_py.corr(numeric_only = True, method = "spearman")

USArrests_cor_py = USArrests_cor_py.rename_axis("var1").reset_index()

print(USArrests_cor_py.head())

USArrests_pear_py = pd.melt(USArrests_cor_py,
        id_vars=['var1'],
        value_vars=['murder', 'assault', 'urban_pop', 'robbery'],
        var_name='var2',
        value_name='cor').round(3)

print(USArrests_pear_py.head())

plt2 = (ggplot(USArrests_pear_py,
        aes(x = "var1", y = "var2", fill = "cor")) +
     geom_tile() +
     geom_text(aes(label = "cor"),
               colour = "white",
               size = 10))

print(plt2)