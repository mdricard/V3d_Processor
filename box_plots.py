import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
sns.set_theme(style="ticks", palette="Dark2")
all_data = pd.read_csv('D:/Alexis_Stats/Alexis_Stats_4_subjects.csv',delimiter=',')
#all_data = pd.read_csv('D:/Alexis_Stats/Alexis_Stats.csv',delimiter=',')
#print(data.head())
hoka = all_data[all_data["shoe"] == 0]
hoka_level = hoka[hoka["incline"] == 0]
hoka_uphill = hoka[hoka["incline"] == 1]
hoka_downhill = hoka[hoka["incline"] == 2]
print(hoka_downhill.head())
barefoot = all_data[all_data["shoe"] == 1]
barefoot_level = barefoot[barefoot["incline"] == 0]
barefoot_uphill = barefoot[barefoot["incline"] == 1]
barefoot_downhill = barefoot[barefoot["incline"] == 2]
print(hoka.head())

#--------------------------   Hoka  Level  -----------------------------
ax = sns.boxplot(x="speed", y="add_mom",
            hue='subject', legend=False,
            data=hoka_level)
sns.despine(offset=10, trim=True)
ax.set_title("Hoka Neutral by Speed")        # Add title and labels
ax.set_xticklabels(['0.5', '0.8', '1.2'])  # Custom
ax.set_xlabel("Walking Speed (m/s)")
ax.set_ylabel("Adduction Moment (%BWxHt)")
plt.show()
#--------------------------   Hoka  Uphill  -----------------------------
ax = sns.boxplot(x="speed", y="add_mom",
            hue='subject', legend=False,
            data=hoka_uphill)
sns.despine(offset=10, trim=True)
ax.set_title("Hoka Uphill by Speed")        # Add title and labels
ax.set_xticklabels(['0.5', '0.8', '1.2'])  # Custom
ax.set_xlabel("Walking Speed (m/s)")
ax.set_ylabel("Adduction Moment (%BWxHt)")
plt.show()
#--------------------------   Hoka  Downhill  -----------------------------
ax = sns.boxplot(x="speed", y="add_mom",
            hue='subject', legend=False,
            data=hoka_downhill)
sns.despine(offset=10, trim=True)
ax.set_title("Hoka Downhill by Speed")        # Add title and labels
ax.set_xticklabels(['0.5', '0.8', '1.2'])  # Custom
ax.set_xlabel("Walking Speed (m/s)")
ax.set_ylabel("Adduction Moment (%BWxHt)")
plt.show()


#--------------------------   Barefoot  Level  -----------------------------
ax = sns.boxplot(x="speed", y="add_mom",
            hue='subject', legend=False,
            data=barefoot_level)
sns.despine(offset=10, trim=True)
ax.set_title("Barefoot Neutral by Speed")        # Add title and labels
ax.set_xticklabels(['0.5', '0.8', '1.2'])  # Custom
ax.set_xlabel("Walking Speed (m/s)")
ax.set_ylabel("Adduction Moment (%BWxHt)")
plt.show()
#--------------------------   Barefoot  Uphill  -----------------------------
ax = sns.boxplot(x="speed", y="add_mom",
            hue='subject', legend=False,
            data=barefoot_uphill)
sns.despine(offset=10, trim=True)
ax.set_title("Barefoot Uphill by Speed")        # Add title and labels
ax.set_xticklabels(['0.5', '0.8', '1.2'])  # Custom
ax.set_xlabel("Walking Speed (m/s)")
ax.set_ylabel("Adduction Moment (%BWxHt)")
plt.show()
#--------------------------   Barefoot  Downhill  -----------------------------
ax = sns.boxplot(x="speed", y="add_mom",
            hue='subject', legend=False,
            data=barefoot_downhill)
sns.despine(offset=10, trim=True)
ax.set_title("Barefoot Downhill by Speed")        # Add title and labels
ax.set_xticklabels(['0.5', '0.8', '1.2'])  # Custom
ax.set_xlabel("Walking Speed (m/s)")
ax.set_ylabel("Adduction Moment (%BWxHt)")
plt.show()