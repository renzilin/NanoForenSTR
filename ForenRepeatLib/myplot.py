import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams.update({'font.size': 11,
                 'axes.titlesize': 15,     # title font size
                 'axes.labelsize': 13,     # xlabel/ylabel font size
                 'axes.linewidth': 1.5,    # axis width
                 'xtick.labelsize': 11,    # xticks font size
                 'ytick.labelsize': 11,    # yticks font size
                 'xtick.major.size': 6,    # xticks length
                 'xtick.major.width': 1.5, # xticks width
                 'ytick.major.size': 6,    # yticks length
                 'ytick.major.width': 1.5, # yticks width
                 })

## plot the histogram of repeat copy number distribution in each locus
def func_save_plot(count_lst, str_name, OUTDIR):
    fig = plt.figure(figsize = (6.2, 4.2))
    plt.hist(count_lst, bins=100)

    ax = plt.gca()
    x_major_locator = plt.MultipleLocator(1)
    ax.xaxis.set_major_locator(x_major_locator)

    plt.xlabel('copy number')
    plt.xlim(0)
    plt.xticks(rotation=-45)
    plt.ylabel('frequency')
    plt.title('copy number from reads: %s' % str_name)
    fig.savefig('%s/copy_number_%s.png' % (OUTDIR, str_name), dpi=300)
    plt.close(fig)
    
    return fig