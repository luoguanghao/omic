import numpy as np
import pandas as pd
import matplotlib.pyplot as plt



def imshow_with_edge(df):
    '''
    有边框的热图
    '''
    fig, ax = plt.subplots(figsize=(10,4))

    im = plt.imshow(df, cmap='Blues', interpolation='none', vmin=0, vmax=1, aspect='equal')

    ax.set_xticks(np.arange(df.shape[1]))
    ax.set_yticks(np.arange(df.shape[0]))

    ax.set_xticklabels(df.columns)
    ax.set_yticklabels(df.index)

    plt.setp(ax.get_xticklabels(), rotation=-85, ha="left",
            rotation_mode="anchor")


    def rect(pos):
        r = plt.Rectangle(pos-0.5, 1,1, facecolor="none", edgecolor="r", linewidth=2)
        plt.gca().add_patch(r)

    x,y = np.meshgrid(np.arange(df.shape[0]),np.arange(df.shape[1]))
    m = np.c_[y.T[df.astype(bool)],x.T[df.astype(bool)]]
    for pos in m:
        rect(pos)

    plt.show()



from matplotlib import gridspec
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
def combine_graph():
    '''
    just a demo
    '''

    index = np.random.randint(0,20,size=200)

    fig = plt.figure(figsize=(5, 5)) # 定义图尺寸 figsize=(横长×竖长)
    gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1]) # 定义布局，几行、几列，图的高度的比值
    ax0 = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[1])

    ax0.hist(log2_rpkm_expr_mat['CJX20200428'])










