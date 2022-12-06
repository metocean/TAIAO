from matplotlib import pyplot as plt
import numpy as np

import matplotlib.gridspec as gridspec
from matplotlib.widgets import SpanSelector,RectangleSelector
from matplotlib.dates import date2num


def check_timeseries(y):
    
    nonans = np.nonzero(~np.isnan(y))[0]
    before = nonans-1
    before[before<0] = 0
    after = nonans+1
    after[after>len(y)-1] = len(y)-1
    y_befores = y[before]
    y_afters = y[after]

    inter = len(np.nonzero(np.logical_and(np.isnan(y_befores),
                                          np.isnan(y_afters)))[0])
    tot = len(nonans)
    if tot == 0:
        return False
    else:
        return inter/tot >= 0.5
        

def get_freq(df):
    dt = np.abs(np.diff(df.index))/1000000000
    return dt


def display_info(df,):
    freq = get_freq(df)
    print('min timestep: %i' % freq.min() )
    print('max timestep: %i' % freq.max() )
    print('start time: %s' % df.index[0].strftime('%Y-%m-%d %H:%M:%S ') )
    print('end time: %s' % df.index[-1].strftime('%Y-%m-%d %H:%M:%S ')  )

    info = {
        'Min interval': freq.min(),
        'Max interval': freq.max(),
        'Start time': df.index[0].strftime('%Y-%m-%d %H:%M:%S '),
        'End time':df.index[0].strftime('%Y-%m-%d %H:%M:%S ')
    }
    return info


def plot_graph(dfs,
               variables_to_plot,
               titles,
               ylim=None,
               xlim=None,
               save=None,
               eventname=None,
               show=True):
    """
    Plot variables "variables_to_plot" for the "dfs" dataframes an
    interactive graph and allow the user to manually select data points
    they want to remove. 
    Returns the coordinate of the selected points.

    save (str): path to which to save the plot.
    eventname (str): 
    """
    answer = {} 
    custom_lim = np.ndarray((2,1))
    def printspanselector(xmin, xmax):
        if isinstance(xmin, np.datetime64):
                xmin = date2num(xmin)
                xmax = date2num(xmax)

        ind = np.logical_and(date2num(df.index) >= xmin,
                             date2num(df.index) <= xmax)
        axs[-1].plot(df.index[ind],
                     df[v][ind],
                     'r+',
                     gid='selected')
        if 'delete' not in answer:
            answer['delete'] = []
        answer['delete'].append([xmin, xmax])

    def onclick(event, ax):
        # Only clicks inside this axis are valid.
        try: # use try/except in case we are not using Qt backend
            zooming_panning = ( fig.canvas.cursor().shape() != 0 ) # 0 is the arrow, which means we are not zooming or panning.
        except:
            zooming_panning = False
        if zooming_panning: 
            return

        if eventname == 'choose':
            plt.close('all')
            answer['axes'] = event.inaxes.get_title()
             
    
    plt.close('all')
    fig = plt.figure(figsize = (11.69*(2/3), 8.27/2))
    gs1 = gridspec.GridSpec(len(dfs), 1)
    axs=[]
    for n, df in enumerate(dfs):
        if n == 0:
            axs.append(fig.add_subplot(gs1[n],) )
        else:
            axs.append(fig.add_subplot(gs1[n], sharex=axs[0]) )

        plt.title(titles[n])

        for i,v in enumerate(variables_to_plot):
            if v in df:
                if not check_timeseries(df[v].values):
                    typ = '-'
                else:
                    typ = '+'

                axs[-1].plot(df.index,
                             df[v],
                             typ,
                             label=v)

    if xlim:
        plt.xlim(xlim[0], xlim[1])

    if ylim:
        plt.ylim(ylim[0], ylim[1])

    plt.legend(bbox_to_anchor=(0.0, -.5, 1, 0),
               loc="lower left",
               mode="expand", ncol=2)
    
    fig.autofmt_xdate()

    if save:
        plt.savefig(save,
                    format='png',
                    bbox_inches="tight")

    if eventname=='choose':
        cid = fig.canvas.mpl_connect('button_press_event', lambda event: onclick(event, axs[-1]))
        plt.show()
    elif eventname=='clean':
        selector = SpanSelector(axs[-1],
                                onselect=printspanselector,
                                direction='horizontal',
                                useblit=True,
                                span_stays=False,
                                button=1,
                                rectprops={'facecolor': 'grey', 'alpha': 0.3})
        fig.canvas.mpl_connect('key_press_event',
                               selector)
        plt.show()
    else:
        if show:
            plt.show()

    return answer
