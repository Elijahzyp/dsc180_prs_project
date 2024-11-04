
import matplotlib.pyplot as plt
import seaborn as sns
import os
from scipy import stats
import pandas as pd
main_dir = os.getcwd()[:-3]


def plot_prs_dist(file_name, trait, person):
    
    prs_data = pd.read_csv(main_dir + "PRS/" + file_name + ".profile", delim_whitespace=True)
    

    sns.set(style="whitegrid")
    plt.figure(figsize=(10, 6))
    sns.histplot(prs_data['SCORE'], kde=True, color="skyblue", bins=30)
    plt.title(f'Distribution of PRS for {trait}')
    plt.xlabel('Polygenic Risk Score')
    plt.ylabel('Frequency')


    person_prs = prs_data[prs_data['IID'] == person]['SCORE'].iloc[0]
    plt.axvline(person_prs, color='red', linestyle='dashed', linewidth=2)
    plt.text(person_prs, plt.ylim()[1]*0.9, person, color = 'red')

    plt.legend(['1KG Individuals', person])
    plt.show()

