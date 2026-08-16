import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt
import seaborn as sns
from adjustText import adjust_text


# -- VOLCANO PLOT -- 

def __map_color(results_df, lfc_threshold=1, nlog10_threshold=2):
    '''
    recall you were goind to turn this into a for loop based (or maybe numpy broadcasting) so that we can add a threshold to l2fc and nlog10

    '''
    pass
    
    color = []
    for log2FoldChange, nlog10 in zip(results_df['log2FoldChange'], results_df['nlog10']):

        if abs(log2FoldChange) < lfc_threshold or nlog10 < nlog10_threshold:
            color.append("Low DE")
        elif abs(log2FoldChange) >= lfc_threshold and nlog10 > nlog10_threshold:
            if log2FoldChange >=lfc_threshold:
                color.append("Overexpressed")
            elif log2FoldChange <=-lfc_threshold:
                color.append("Underexpressed")
        return color

def volcano_plot(results_df, ax, lfc_threshold=1, nlog10_threshold=2):
    '''
    'results_df' must have columns: 'padj', 'nlog10'. To show gene symbols at each point the index of 'results_df' must be differentially expressed gene symbols.
    '''
    sigs = results_df.copy()

    sigs["nlog10"] = -np.log10(sigs.padj)
    sigs['color'] = __map_color(results_df, lfc_threshold=lfc_threshold, nlog10_threshold=nlog10_threshold)


    pallete = {"Overexpressed": "#C0392B", "Underexpressed" :"#2980B9", "Low DE": "#95A5A6" }

    sns.scatterplot(data=sigs, x="log2FoldChange", y="nlog10", hue="color", palette=pallete, ax=ax)

    texts = []
    for x, y, gene in zip(sigs["log2FoldChange"], sigs["nlog10"], sigs.index):
        if abs(x) >= 4 and y > 75:
            texts.append(plt.text(x, y, gene))

    adjust_text(texts)
    ax.legend(title="")
    ax.axhline(nlog10_threshold, zorder=0, c="k", lw=1, ls="--")
    ax.axvline(lfc_threshold, zorder=0, c="k", lw=1, ls="--")
    ax.axvline(-lfc_threshold, zorder=0, c="k", lw=1, ls="--")

