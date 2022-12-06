import numpy as np
from scipy.signal import find_peaks
from numpy.matlib import repmat

from toto.core.toolbox import (
    dir_interval,
    get_increment,
    get_number_of_loops,
    degToCompass
)

import wafo
wafo.plotbackend.plotbackend.interactive(False)
import wafo.stats as ws

import matplotlib.pyplot as plt


def do_fitting(mag,
               fitting,
               method,
               loc=None):

    """ This function does the main fitting for the various distribution and method.
        I had issues with Gumbel and it seems to work when I switch the order of the parameter,
        Looks like GenPareto shape needs to be inverted. not sure why
    """

    if loc is None:
        loc=np.nanmin(mag)*.999

    if fitting.lower() == 'weibull':
        phat = ws.weibull_min.fit2(mag,
                                   floc=loc,
                                   method=method.lower(),
                                   alpha=0.05)
        scale = phat.par[-1]
        shape = phat.par[0]
    elif fitting.lower() == 'gumbel':
        phat = ws.gumbel_r.fit2(mag,
                                method=method,
                                alpha=0.05)
        scale = phat.par[0]
        shape = phat.par[-1]
    elif fitting.lower() == 'gpd':
        phat = ws.genpareto.fit2(mag,
                                 floc=loc,
                                 method=method.lower(),
                                 alpha=0.05)
        scale = phat.par[-1]
        shape = phat.par[0]*-1

    elif fitting.lower() == 'gev':
        phat = ws.genextreme.fit2(mag,
                                  floc=loc,
                                  method=method.lower(),
                                  alpha=0.05)
        scale = phat.par[-1]
        shape = phat.par[0]

    else:
        assert 'Fitting %s not recognize' % fitting

    return phat, scale, shape


def get_mag_for_water_elevation(el_res,
                                et,
                                threshold):
        """
        Here we are re-creating a timeseries by joining the residual and the tide together.
        We do so by doing a joint prob btw the 2 timeseries
        """
        
        ntime = np.ceil( len(et) / len(el_res) )
        el_resJP = repmat(el_res, int(ntime), 1).ravel()
        el_resJP = el_resJP[:len(et)]

        #joint prob tide/surge
        d = 0.01;
        jp, X_interval, Y_interval = joint_prob_annual(el_resJP, et, d, d);
        X = np.array([])
        Y = np.array([])
        n = 0 #%transforms the 2D joint prob into a 1D prob in term of total sea level
        for i in range(0, len(X_interval)-1):
            for j in range(0, len(Y_interval)-1):
                if jp[j,i] != 0:
                    X = np.append(X, X_interval[i] + Y_interval[j] + d)
                    Y = np.append(Y, jp[j,i])
                    n = n+1

        ss = np.arange(0, max(X)+d, d)
        f = np.ones((len(ss)-1,))*np.nan
        for i in range(0, len(ss)-1):
            f[i] = np.sum(Y[np.logical_and( X >= ss[i], X < ss[i+1] )])


        ss = ss[:-1] + d/2.
        F = np.cumsum(f) / np.sum(f)

        ind = (F > threshold/100).nonzero()[0][0]
        ss = ss[ind:]
        F = F[ind:]
        f = np.diff(F)
        multiplicator = 1./np.min(f[f != 0])  #n_events/(1-F(1));
        mag = np.array([])
        for i in range(0, len(ss)-1):# create
            mag = np.append(mag,
                            ss[i+1] - d/2 + d*np.random.rand( int( np.round( np.sqrt(f[i] * multiplicator) ) ), 1) )

        return mag, el_resJP, et
        


def joint_prob_annual(Xdata,
                      Ydata,
                      dx,
                      dy):
    """
    This will do a joint probablilty btw Xdata and Ydata.
    It will go from min to max every dx or dy.
    the JP result is between 0 and 1.
    """

    X_interval = np.arange(np.round(np.nanmin(Xdata), 2) - dx,
                           np.round(np.nanmax(Xdata), 2) + dx,
                           dx)
    Y_interval = np.arange(np.round(np.nanmin(Ydata), 2) - dy,
                           np.round(np.nanmax(Ydata), 2) + dy,
                           dy)

    occurrence = np.zeros((len(Y_interval)-1, len(X_interval)-1))
    for k in range(0, len(X_interval)-1):
        index1 = np.logical_and(Xdata > X_interval[k],
                                Xdata <= X_interval[k+1])
        for m in range(0, len(Y_interval)-1):
            index2=(np.logical_and(Ydata[index1] > Y_interval[m],
                                   Ydata[index1] <= Y_interval[m+1])).nonzero()[0]
            occurrence[m,k] = len(index2) / len(Xdata)

    return occurrence, X_interval, Y_interval


def do_EVA_water(peaks_index,
                 dfout,
                 mag,
                 tide,
                 fitting,
                 method,
                 time_blocking,
                 threshold,
                 tide_join_prob=False):
    """ 
    do the distribution on either Tide+Storm surg or 
    just storm surge depending on the time blocking (yearly,monthly ....)
    """

    eva_stats = {}
    el_res = dfout[mag]
    et = dfout[tide].values

    number_of_loops, identifiers, month_identifier = get_number_of_loops(time_blocking)
    months = dfout.index.month
    idx = np.arange(0, len(months))
    mag_values = dfout[mag]
    for j in range(0, number_of_loops):
    #Pull out relevant indices for particular month/months
        if identifiers[j] in peaks_index:
            peak = peaks_index[identifiers[j]]['Omni']
            index = np.in1d(months, month_identifier[j])

            if identifiers[j] not in eva_stats:
                eva_stats[identifiers[j]] = {}
                eva_stats[identifiers[j]]['Omni'] = {}

            if mag not in eva_stats[identifiers[j]]['Omni']:
                eva_stats[identifiers[j]]['Omni'][mag] = {}

            if tide_join_prob:
                y, el_resJP, et = get_mag_for_water_elevation(dfout[mag][index].values,
                                                              dfout[tide].values,
                                                              threshold)
            else:
                y = dfout[mag][peak].values
                y = y[~np.isnan(y)]
                el_resJP = y

            loc = min(y) - 0.01;
            phat, scale, shape = do_fitting(y, fitting, method, loc=loc)
            eva_stats[identifiers[j]]['Omni'][mag]['phat'] = phat
            eva_stats[identifiers[j]]['Omni'][mag]['scale'] = scale
            eva_stats[identifiers[j]]['Omni'][mag]['shape'] = shape
            eva_stats[identifiers[j]]['Omni'][mag]['el_resJP'] = el_resJP
            eva_stats[identifiers[j]]['Omni'][mag]['et'] = et

    return eva_stats


def get_peaks(dfout,
              mag,
              drr=None,
              time_blocking='Annual',
              directional_interval=[0, 360],
              peaks_options={},
              min_peak=30):
    """
    get peaks in a timeseries and save them in a dict depending on
    the time blocking and directional bins if any
    """

    peaks_index= {}
    if drr == None or drr == '':
        drr = 'direction_optional'
        dfout['direction_optional'] = np.ones((len(dfout[mag].values),)) * 20

    
    number_of_loops,identifiers,month_identifier =\
        get_number_of_loops(time_blocking)
    months = dfout.index.month
    idx = np.arange(0,len(months))
    drr_values = dfout[drr]
    mag_values = dfout[mag]

    for j in range(0, number_of_loops):
    #Pull out relevant indices for particular month/months
        index1 = np.in1d(months,
                         month_identifier[j])
        peaks_index[identifiers[j]] = {}

        for jj in range(0, len(directional_interval)):

            if jj == len(directional_interval) - 1:
                index2 = drr_values > -1
                dir_label = 'Omni'
            else:
                dir_label = degToCompass([directional_interval[jj],
                                          directional_interval[jj+1]])
                if directional_interval[jj+1] <= directional_interval[jj]:
                    index2 = np.logical_or(drr_values > directional_interval[jj],
                                           drr_values <= directional_interval[jj+1])
                else:
                    index2 = (drr_values>directional_interval[jj]) &\
                              (drr_values<=directional_interval[jj+1])
            index = np.logical_and(index1, index2)
            if np.any(index):
                tmp = mag_values.values.copy()
                tmp[~index] = 0
                pk_idx = find_peaks(tmp, **peaks_options)[0]

                if len(pk_idx) > min_peak:
                    peaks_index[identifiers[j]][dir_label] = pk_idx
    return peaks_index
