import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
sns.set_theme(style="ticks", palette="Dark2")

all_data = pd.read_csv('D:/Alexis_Stats/Alexis_Stats_4_subjects.csv',delimiter=',')
#print(data.head())
hoka = all_data[all_data["shoe"] == 0]
hoka_level = hoka[hoka["incline"] == 0]
hoka_uphill = hoka[hoka["incline"] == 1]
hoka_downhill = hoka[hoka["incline"] == 2]
barefoot = all_data[all_data["shoe"] == 1]
print(hoka.head())

#--------------------------   Hoka  Uphill  -----------------------------
ax = sns.boxplot(x="speed", y="add_mom",
            hue='subject',
            data=hoka_uphill)
sns.despine(offset=10, trim=True)
ax.set_title("Hoka Uphill by Speed")        # Add title and labels
ax.set_xticklabels(['0.5', '0.8', '1.2'])  # Custom
ax.set_xlabel("Walking Speed (m/s)")
ax.set_ylabel("Adduction Moment (%BWxHt)")
plt.show()
#--------------------------   Hoka  Downhill  -----------------------------
ax = sns.boxplot(x="speed", y="add_mom",
            hue='subject',
            data=hoka_downhill)
sns.despine(offset=10, trim=True)
ax.set_title("Hoka Downhill by Speed")        # Add title and labels
ax.set_xticklabels(['0.5', '0.8', '1.2'])  # Custom
ax.set_xlabel("Walking Speed (m/s)")
ax.set_ylabel("Adduction Moment (%BWxHt)")
plt.show()