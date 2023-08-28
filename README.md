# beh_ecg_phase
This code 1) plots events times of hits, misses, false alarms and correct rejections as well as event times of completed and not completed trials on the circular rose plot
          2) draws a mean direction vector, considering its strength on the plots,
          3) counts event mean phases and performs o-tests for all sessions (date) and all trials together 
          4) performs Rayleigh test and counts d-prime and criterion for all trials. 
          5) Plots means for all sesions on the polar plot


          The code uses circ_mean.m, ig_rose.m and testsim_dprime.m
          I changed ig_rose a bit, cause it didn't work on matlab2014, so I attached this also just in case. 
