
# coding: utf-8

# ## Observations and Insights 

# In[106]:


# Dependencies and Setup
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as st
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import linregress
# Study data files
mouse_metadata_path = "data/Mouse_metadata.csv"
study_results_path = "data/Study_results.csv"

# Read the mouse data and the study results
mouse_metadata = pd.read_csv(mouse_metadata_path)
study_results = pd.read_csv(study_results_path)
mouse_metadata = mouse_metadata.dropna()
study_results = study_results.dropna()
# Combine the data into a single dataset
merge_df = pd.merge(mouse_metadata, study_results, on="Mouse ID")
# Display the data table for preview
merge_df.head()


# In[2]:


# Checking the number of mice.
len(merge_df["Mouse ID"].unique())


# In[3]:


# Getting the duplicate mice by ID number that shows up for Mouse ID and Timepoint. 
df = pd.DataFrame(merge_df.duplicated(["Mouse ID","Timepoint"]))
# to compare
booleanDictionary = {True: 'TRUE', False: 'FALSE'}
df = df.replace(booleanDictionary)
df = df.loc[df[0] == "TRUE", :]
df.index.values


# In[4]:


dupt_mice = merge_df.iloc[[909, 911, 913, 915, 917],[0]]
dupt_mice


# In[5]:


# Optional: Get all the data for the duplicate mouse ID. 
merge_df.loc[merge_df["Mouse ID"] == "g989", :]


# In[6]:


# Create a clean DataFrame by dropping the duplicate mouse by its ID.
clean = merge_df.drop(merge_df[merge_df["Mouse ID"]=="g989"].index)
clean


# In[7]:


# Checking the number of mice in the clean DataFrame.
len(clean["Mouse ID"].unique())


# ## Summary Statistics

# In[29]:


# Generate a summary statistics table of mean, median, variance, standard deviation, and SEM of the tumor volume for each regimen

# Use groupby and summary statistical methods to calculate the following properties of each drug regimen: 
# mean, median, variance, standard deviation, and SEM of the tumor volume. 

grouped_df = clean.groupby(['Drug Regimen'])
mean = grouped_df["Tumor Volume (mm3)"].mean()
med = grouped_df["Tumor Volume (mm3)"].median()
var = grouped_df["Tumor Volume (mm3)"].var()
std = grouped_df["Tumor Volume (mm3)"].std()
sem = grouped_df["Tumor Volume (mm3)"].sem()
# Assemble the resulting series into a single summary dataframe.
sum1 = pd.DataFrame({"mean":mean,
                    "med":med,
                    "var":var,
                    "std":std,
                    "sem":sem,
                    })
sum1


# In[31]:


# Generate a summary statistics table of mean, median, variance, standard deviation, and SEM of the tumor volume for each regimen

# Using the aggregation method, produce the same summary statistics in a single line
df_agg = clean.groupby(['Drug Regimen']).agg({"Tumor Volume (mm3)":['mean', 'median', 'var', 'std', 'sem']})
df_agg


# ## Bar and Pie Charts

# In[54]:


# Generate a bar plot showing the total number of timepoints for all mice tested for each drug regimen using Pandas.
num_tp = clean.groupby(["Drug Regimen"]).agg({"Timepoint":["count"]})
num_tp['Drug Regimen'] = num_tp.index
myplot = num_tp.plot(kind='bar')
print(myplot)


# In[59]:


# Generate a bar plot showing the total number of timepoints for all mice tested for each drug regimen using pyplot.
plt.bar(num_tp['Drug Regimen'], num_tp['Timepoint'])
plt.show()


# In[82]:


# Generate a pie plot showing the distribution of female versus male mice using Pandas
gender = clean.groupby(['Sex']).nunique("Mouse ID")
myplot_pie = gender.plot.pie(y='Mouse ID', figsize=(5, 5))
print(myplot_pie)


# In[83]:


# Generate a pie plot showing the distribution of female versus male mice using pyplot
plt.pie(gender['Mouse ID'])
plt.show()


# ## Quartiles, Outliers and Boxplots

# In[86]:


# Calculate the final tumor volume of each mouse across four of the treatment regimens:  
# Capomulin, Ramicane, Infubinol, and Ceftamin
treatment = merge_df.loc[(merge_df["Drug Regimen"] == "Capomulin")|(merge_df["Drug Regimen"] == "Ramicane")|(merge_df["Drug Regimen"] == "Infubinol")|(merge_df["Drug Regimen"] == "Ceftamin"),:]

# Keep the last timepoint of each Mouse ID
treatment = treatment.drop_duplicates("Mouse ID", keep = "last")
treatment


# In[87]:


treatment.dtypes


# In[ ]:


# Put treatments into a list for for loop (and later for plot labels)


# Create empty list to fill with tumor vol data (for plotting)


# Calculate the IQR and quantitatively determine if there are any potential outliers. 

    
    # Locate the rows which contain mice on each drug and get the tumor volumes
    
    
    # add subset 
    
    
    # Determine outliers using upper and lower bounds
    


# In[ ]:


# Generate a box plot of the final tumor volume of each mouse across four regimens of interest


# ## Line and Scatter Plots

# In[96]:


# Generate a line plot of tumor volume vs. time point for a mouse treated with Capomulin
capomulin = merge_df.loc[(merge_df["Drug Regimen"] == "Capomulin"),:]
g_c = capomulin.groupby("Timepoint").agg({"Tumor Volume (mm3)":["mean"]})
g_c['Timepoint'] = g_c.index
x_axis = np.array(g_c["Timepoint"])
y = np.array(g_c["Tumor Volume (mm3)"])
capomulin_p = plt.plot(x_axis, y, marker="s", color="Red", linewidth=1, label="Capomulin")


# In[110]:


# Generate a scatter plot of average tumor volume vs. mouse weight for the Capomulin regimen
g_c_w = capomulin.groupby("Weight (g)").agg({"Tumor Volume (mm3)":["mean"]})
g_c_w['Weight (g)'] = g_c_w.index
x = np.array(g_c_w["Weight (g)"])
y = np.array(g_c_w["Tumor Volume (mm3)"])
plt.scatter(x, y, marker="o", facecolors="yellow", edgecolors="blue", alpha=0.5)
plt.show()


# ## Correlation and Regression

# In[122]:


# Calculate the correlation coefficient and linear regression model 
# for mouse weight and average tumor volume for the Capomulin regimen

y_values = g_c_w["Tumor Volume (mm3)"]
y_values['Weight (g)'] = y_values.index
x_values = y_values["Weight (g)"]
y_values = y_values["mean"]
(slope, intercept, rvalue, pvalue, stderr) = linregress(x_values, y_values)
regress_values = x_values * slope + intercept
line_eq = "y = " + str(round(slope,2)) + "x + " + str(round(intercept,2))
plt.scatter(x_values,y_values)
plt.plot(x_values,regress_values,"r-")
plt.annotate(line_eq,(0,50),fontsize=15,color="red")
plt.xlabel('Weight (g)')
plt.ylabel('Tumor Volume (mm3)')
plt.show()

