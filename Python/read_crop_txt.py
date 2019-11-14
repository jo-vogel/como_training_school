#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 14:20:24 2019

@author: christoph
"""

import numpy as np
import glob
import matplotlib.pyplot as plt

K = 273.16
def open_txt_Damocles(var):
    # by Karin
    # open file
    filename = var
    file = open(filename)
    lines = file.readlines()
    file.close()

    # create array
    slab = np.full((1600,len(lines[0].split())),np.nan)

    # loop over all ensemble members
    for i in range(0,1600):
        line = lines[i+1].split()[1:]
        data = np.array(line).astype('float')
        slab[i,:] = data

    # done
    return slab


def plot_mean_with_std(x, gr1_alpha=0.3, gr1_color='b', plt_std=True):
    plt.plot(range(0,x.shape[1]), np.mean(x, axis=0),
             linewidth=2, color=gr1_color)
    if plt_std:
        plt.fill_between(
            range(0,x.shape[1]),
            np.mean(x, axis=0) - np.std(x, axis=0),
            np.mean(x, axis=0) + np.std(x, axis=0),
            alpha=gr1_alpha, color=gr1_color, linewidth=0, label='_nolegend_')
    plt.xticks([1,  31,  62,  93, 121, 152, 182, 213, 243, 274, 305],
               ['N','D', 'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S'])
    for pos in ['right','top']:
        plt.gca().spines[pos].set_visible(False)


def running_sum(x, running_days=30):
    x1 = np.zeros(x.shape)
    x1[:,:] = np.nan
    for year in range(0,x.shape[0]):
        for day in range(running_days,x.shape[1]):
            x1[year,day] = np.sum(x[year,day-running_days:day])
    return x1
            

# ============================================================================
# countries: FR, US, CN         
country = "US"
figsave = 1

path = '/Users/christoph/Desktop/DAMOCLES_training_school/WorkingGroup1/txt_data/'
txtfiles = glob.glob(path + country + "/*.txt")
txtfiles.sort()



# read files
crop = open_txt_Damocles(txtfiles[0])
dps = open_txt_Damocles(txtfiles[1])
pr = open_txt_Damocles(txtfiles[2])
wind = open_txt_Damocles(txtfiles[6])
rsds = open_txt_Damocles(txtfiles[3])
tmax = open_txt_Damocles(txtfiles[4])
tmin = open_txt_Damocles(txtfiles[5])

dpd = tmax-dps

#colors
gr1_orange = '#d95f02' # '#FFD39B'
gr1_blue =  '#7570b3' # '#98F5FF'
gr1_green = '#1b9e77'# '#7FFF00'
gr1_alpha = 0.3
gr1_lw = 2

# Positive extremes: yield value > 1 - alpha percentiles
alpha = 0.05
start_date = {'FR':304, 'US':276, 'CN':272} # 304 is 31st Oct
sowingdate = start_date[country]


months=[31, 30, 31, 31, 28, 31, 30 ,31, 30, 31, 31, 30]

#TODO implement this for US
xticks_country = {'US':
    ([28,  59,  89, 120, 151, 179, 210, 240, 271, 301, 332, 363, 393],
     ['N','D', 'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S']),
     'FR':
    ([1,  31,  62,  93, 121, 152, 182, 213, 243, 274, 305],
     ['N','D', 'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S'])
    }


pos_extreme_indices = crop[:,0] > np.quantile(crop[:,0], 1 - alpha)
neg_extreme_indices = crop[:,0] < np.quantile(crop[:,0], alpha)
avr_indices = np.logical_and((pos_extreme_indices == False),
                             (neg_extreme_indices == False))

tmax_pos_extr = tmax[pos_extreme_indices]
tmax_neg_extr = tmax[neg_extreme_indices]
tmax_av = tmax[avr_indices]

tmin_pos_extr = tmin[pos_extreme_indices]
tmin_neg_extr = tmin[neg_extreme_indices]
tmin_av = tmin[avr_indices]

# =======================
# MAX AND MIN TEMPERATURE
y_limits=[-5,40]
fig1, (ax1, ax2) = plt.subplots(1, 2)
ax1.plot(np.mean(tmax_neg_extr - K, axis=0), color=gr1_orange, linewidth=gr1_lw)
ax1.fill_between(range(0,tmax_neg_extr.shape[1]),
                 np.mean(tmax_neg_extr-K, axis=0)
                 - np.std(tmax_neg_extr, axis=0),
                 np.mean(tmax_neg_extr-K, axis=0)
                 + np.std(tmax_neg_extr, axis=0),
                 alpha=gr1_alpha, color=gr1_orange, linewidth=0)
ax1.plot(np.mean(tmax_pos_extr - K, axis=0), color=gr1_green, linewidth=gr1_lw)
ax1.fill_between(range(0,tmax_pos_extr.shape[1]),
                 np.mean(tmax_pos_extr-K, axis=0)
                 - np.std(tmax_pos_extr, axis=0),
                 np.mean(tmax_pos_extr-K, axis=0)
                 + np.std(tmax_pos_extr, axis=0),
                 alpha=gr1_alpha, color=gr1_green, linewidth=0)
ax1.plot(np.mean(tmax_av - K, axis=0), color=gr1_blue, linewidth=gr1_lw)
ax1.set_xticks([1,  31,  62,  93, 121, 152, 182, 213, 243, 274, 305])
ax1.set_xticklabels(['N','D', 'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S'])
ax1.set_ylabel(r'Temperature [$^{\circ}$C]')
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.set_ylim(y_limits)
ax1.set_title('Maximum Temperature')

ax2.plot(np.mean(tmin_neg_extr - K, axis=0), color=gr1_orange, linewidth=gr1_lw)
ax2.fill_between(range(0,tmin_neg_extr.shape[1]),
                 np.mean(tmin_neg_extr-K, axis=0)
                 - np.std(tmin_neg_extr, axis=0),
                 np.mean(tmin_neg_extr-K, axis=0)
                 + np.std(tmin_neg_extr, axis=0),
                 alpha=gr1_alpha, color=gr1_orange, linewidth=0)
ax2.plot(np.mean(tmin_pos_extr - K, axis=0), color=gr1_green, linewidth=gr1_lw)
ax2.fill_between(range(0,tmin_pos_extr.shape[1]),
                 np.mean(tmin_pos_extr-K, axis=0)
                 - np.std(tmin_pos_extr, axis=0),
                 np.mean(tmin_pos_extr-K, axis=0)
                 + np.std(tmin_pos_extr, axis=0),
                 alpha=gr1_alpha, color=gr1_green, linewidth=0)
ax2.plot(np.mean(tmin_av - K, axis=0), color=gr1_blue, linewidth=gr1_lw)
ax2.set_ylim(y_limits)
ax2.set_title('Minimum Temperature')
ax2.legend(('5pct yrs lowest yield',
            '5pct yrs highest yield',
            '5 - 95pct av. yield',),
            frameon=False)
ax2.set_xticks([1,  31,  62,  93, 121, 152, 182, 213, 243, 274, 305])
ax2.set_xticklabels(['N','D', 'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S'])
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
if figsave:
    fig1.savefig('Max_min_Temp_entire_year_' + country + '.pdf')





dps_pos_extr = dps[pos_extreme_indices]
dps_neg_extr = dps[neg_extreme_indices]
dps_av = dps[avr_indices]

# ===============
# DEWPOINT FIGURE
fig_dewpd = plt.figure()
plot_mean_with_std(tmax_neg_extr-dps_neg_extr, gr1_alpha, gr1_orange)
plot_mean_with_std(tmax_pos_extr-dps_pos_extr, gr1_alpha, gr1_green)
plot_mean_with_std(tmax_av-dps_av, gr1_alpha, gr1_blue, plt_std=False)
plt.xticks([1,  31,  62,  93, 121, 152, 182, 213, 243, 274, 305],
               ['N','D', 'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S'])
plt.ylabel('Dew point depression [K]')
plt.legend(['5pct yrs lowest yield',
            '5pct yrs highest yield',
            '5 - 95pct av. yield'],
            frameon=False, loc=2)
if figsave:
    fig_dewpd.savefig('Dew_point_dp_year_' + country + '.pdf')

pr_pos_extr = running_sum(pr[pos_extreme_indices], 90)
pr_neg_extr = running_sum(pr[neg_extreme_indices], 90)
pr_av = running_sum(pr[avr_indices], 90)

# ====================
# PRECIPITATION FIGURE
fig_precip = plt.figure()
plot_mean_with_std(pr_neg_extr, gr1_alpha, gr1_orange)
plot_mean_with_std(pr_pos_extr, gr1_alpha, gr1_green)
plot_mean_with_std(pr_av, gr1_alpha, gr1_blue, False)
plt.ylabel(r'3-month precipitation sum [mm]')
plt.legend(['5pct yrs lowest yield',
            '5pct yrs highest yield',
            '5 - 95pct av. yield'],
            frameon=False, loc=2)
plt.xlim([80,pr_av.shape[1]])
# plt.ylim([50,280])
if figsave:
    fig_precip.savefig('3_month_sum_precip_' + country + '.pdf')


dpd_p = dpd[pos_extreme_indices]
dpd_n = dpd[neg_extreme_indices]
dpd_a = dpd[avr_indices]

rsds_p = rsds[pos_extreme_indices]
rsds_n = rsds[neg_extreme_indices]
rsds_a = rsds[avr_indices]

# ======================
# SOLAR RADIATION FIGURE
fig_sr = plt.figure()
plot_mean_with_std(rsds_n, gr1_alpha, gr1_orange)
plot_mean_with_std(rsds_p, gr1_alpha, gr1_green)
plot_mean_with_std(rsds_a, gr1_alpha, gr1_blue, plt_std=False)
plt.ylabel(r'Solar radiation [W/m$^2$]')
plt.legend(['5pct yrs lowest yield',
            '5pct yrs highest yield',
            '5 - 95pct av. yield'],
            frameon=False, loc=2)
if figsave:
    fig_sr.savefig('Solar_rad_year_' + country + '.pdf')
    
