import pandas as pd
import plotnine
# load the data
cortisol_py = pd.read_csv('~/Documents/data/CS1-twopaired.csv')

print(cortisol_py.head())

(ggplot(cortisol_py,
        aes(x = "time",
            y = "cortisol")) +
     geom_boxplot() +
     geom_jitter(width = 0.05) +
     ylab("Cortisol level (nmol/l)"))

