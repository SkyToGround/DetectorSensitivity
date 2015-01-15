#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python
# encoding: utf-8

import numpy as np
import pylab as pl
from matplotlib import rc
from matplotlib.ticker import FixedLocator

def PythonPlot(x, y):
	print "Plotting data"
	x = np.array(x)
	y = np.array(y)
	pl.plot(x, y)
	pl.show()
	return 0

def PythonPlotOld2(x, y):
	import numpy as np
	import matplotlib.pyplot as pl
	from matplotlib import rc
	x = np.array(x)
	y = np.array(y)
	
	rc('font', family='serif')
	fig = pl.figure(figsize=(8, 3.5))
	fig.subplots_adjust(left=0.05, right=0.98, bottom=0.13, top=0.98)
	ax = fig.add_subplot(111)
	ax.plot(x, y, lw = 2, c = "k")
	ax.axis([-5, 5, y.min(), None])
	ax.set_xlabel("Tid (s)")
	ax.set_ylabel("Detektorrespons")
	ax.get_yaxis().set_ticks([])
#	ax.set_yscale('log')
	fig.patch.set_alpha(0.0)
	#ax.get_yaxis().set_visible(False)
	fig.savefig("det_resp.pdf", transparent = True)
	
	return 0

def PythonPlotOld3(x, y):
	import numpy as np
	import matplotlib.pyplot as pl
	from matplotlib import rc
	x = np.array(x)
	y = np.array(y)
	
	rc('font', family='serif')
	fig = pl.figure(figsize=(8, 3.5))
	fig.subplots_adjust(left=0.05, right=0.98, bottom=0.13, top=0.98)
	ax = fig.add_subplot(111)
	ax.plot(x, y, lw = 2, c = "k")
	indices = np.where(x >= -0.0)
	x_limit = x[indices]
	y_limit = y[indices]
	
	indices = np.where(x_limit <= 1.0)
	x_limit = x_limit[indices]
	y_limit = y_limit[indices]
	lower = np.zeros(y_limit.shape)
	#ax.fill_between(x_limit, y_limit, lower)
	for i in range(0, 10):
		y_low = y.min()
		y_high = 0.85
		x = i * 2.0 - 5.0
		ax.plot([x, x], [y_low, y_high], "--", lw = 2, c = "k")
	ax.axis([-4, 4, y.min(), None])
	ax.set_xlabel("Time (s)")
	ax.set_ylabel("Detector response")
	ax.get_yaxis().set_ticks([])
	#	ax.set_yscale('log')
	fig.patch.set_alpha(0.0)
	#ax.get_yaxis().set_visible(False)
	#fig.savefig("det_resp_with_marks_and_fill_2.pdf", transparent = True)
	fig.savefig("time_long.pdf", transparent = True)
	return 0

def PythonPlotOld4(x, y):
	import numpy as np
	import matplotlib.pyplot as pl
	from matplotlib import rc
	x = np.array(x)
	y = np.array(y)
	
	rc('font', family='serif')
	fig = pl.figure(figsize=(8, 3.5))
	fig.subplots_adjust(left=0.05, right=0.98, bottom=0.13, top=0.98)
	ax = fig.add_subplot(111)
	ax.plot(x, y, lw = 2, c = "k")
#	indices = np.where(x >= 0)
#	x_limit = x[indices]
#	y_limit = y[indices]
#	
#	indices = np.where(x_limit <= 1.0)
#	x_limit = x_limit[indices]
#	y_limit = y_limit[indices]
#	lower = np.zeros(y_limit.shape)
#	ax.fill_between(x_limit, y_limit, lower)
	for i in range(0, 25):
		y_low = 0.6
		y_high = 0.85
		x = i * 2.0 - 5.0
		ax.plot([x, x], [y_low, y_high], "--", lw = 2, c = "k")
	ax.axis([-5, 5, y.min(), None])
	ax.set_xlabel("Tid (s)")
	ax.set_ylabel("Detektorrespons")
	ax.get_yaxis().set_ticks([])
	#	ax.set_yscale('log')
	fig.patch.set_alpha(0.0)
	#ax.get_yaxis().set_visible(False)
	fig.savefig("time_alt_1.pdf", transparent = True)
	
	return 0

def PythonPlotNew(x, y):
	import locale
	locale.setlocale(locale.LC_ALL, "sv_SE")
	import numpy as np
	import matplotlib.pyplot as pl
	from matplotlib import rc, rcParams
	x = np.array(x)
	y = np.array(y)
	
	rc('font', family='serif')
	rcParams['axes.formatter.use_locale'] = True
	fig = pl.figure(figsize=(8, 3.5))
	fig.subplots_adjust(left=0.05, right=0.98, bottom=0.13, top=0.98)
	ax = fig.add_subplot(111)
	ax.plot(x, y, lw = 2, c = "k")
	ax.axis([0.5, 3, y.min(), None])
	ax.set_xlabel("Integrationstid (s)")
	ax.set_ylabel("Detekterbar aktivitet")
	ax.get_yaxis().set_ticks([])
#	ax.get_xaxis().set_ticklabels(ticks)
	#	ax.set_yscale('log')
	fig.patch.set_alpha(0.0)
	#ax.get_yaxis().set_visible(False)
	fig.savefig("act_vs_int_time.pdf", transparent = True)
	
	return 0

def MinimisePlot(x, y):
	x = np.array(x)
	y = np.array(y)
	y = y / y.max()
	rc('font', family='sans-serif')
	#rc('text', usetex=True)
	fac = 3.0
	fig = pl.figure(figsize=(1.61803 * fac, fac))
	ax = fig.add_subplot(111)
	fig.subplots_adjust(left=0.09, right=0.96, bottom=0.10, top=0.97)
	
	axis_zeros = np.zeros(y.shape)
	
	ax.plot(x, y, lw = 3, color = "k")
	
	min_index = y.argmin()
	x_at_min = x[min_index]
	min_value = y[min_index]
	
	ax.set_xlabel("Integration time, $\Delta t$")
	ax.set_ylabel("Minimum emission rate, $A$")
	ax.set_xticklabels([])
	ax.set_yticklabels([])
	ax.axis([None, None, 0.725, 1.0])
	ax.annotate('Minimum emission rate,\nOptimal integration time', xy=(x_at_min, min_value),  xycoords='data',
					xytext=(30, -25.0), textcoords='offset points',
					arrowprops=dict(arrowstyle="->",
										 connectionstyle="arc,angleA=0,armA=20,angleB=-90,armB=15,rad=7"),
					)
	pl.savefig("Minimise.pdf")
	pl.show()
	return 0

def ExamplePlot(x, y):
	x = np.array(x)
	y = np.array(y)
	rc('font', family='sans-serif')
	fac = 3.0
	fig = pl.figure(figsize=(1.61803 * fac, fac))
	ax = fig.add_subplot(111)
	fig.subplots_adjust(left=0.14, right=0.96, bottom=0.155, top=0.97)
	#fig.subplots_adjust(left=0.07, right=0.96, bottom=0.155, top=0.97)
	
	axis_zeros = np.zeros(y.shape)
	
	ax.plot(x, y, lw = 3, color = "k")
	
	min_index = y.argmin()
	x_at_min = x[min_index]
	min_value = y[min_index]
	
	ax.set_xlabel("Integration time, $\Delta t$ (s)")
	ax.set_ylabel("Min. emission rate, A (Mn/s)")
	ax.axis([1.5, 10, 0.8, 1.2])
	ax.annotate('Minimum emission rate,\nOptimal integration time', xy=(x_at_min, min_value),  xycoords='data',
					xytext=(10, -25.0), textcoords='offset points',
					arrowprops=dict(arrowstyle="->",
										 connectionstyle="arc,angleA=0,armA=20,angleB=-90,armB=15,rad=7"),
					)
	pl.savefig("Example.pdf", transparent = True)
	#pl.show()
	return 0

def TimeDescPlot(x, y):
	x = np.array(x)
	y = np.array(y)
	y = y / y.max()
	rc('font', family='sans-serif')
	#rc('text', usetex=True)
	fac = 3.0
	fig = pl.figure(figsize=(1.61803 * fac, fac))
	ax = fig.add_subplot(111)
	fig.subplots_adjust(left=0.065, right=0.96, bottom=0.19, top=0.97)
	
	axis_zeros = np.zeros(y.shape)
	
	ax.plot(x, y, lw = 3, color = "k")
	
	zero_index = (x**2).argmin()
	tm_index_p = ((x - 1.0)**2).argmin()
	tm2_index_m = ((x + 0.5)**2).argmin()
	tm2_index_p = ((x - 0.5)**2).argmin()
	
	
	fill1 = ax.fill_between(x[zero_index:tm_index_p], axis_zeros[zero_index:tm_index_p], y[zero_index:tm_index_p], facecolor='gray')
	
	used_hatch = ":"
	fill2 = ax.fill_between(x[tm2_index_m:tm2_index_p], axis_zeros[tm2_index_m:tm2_index_p], y[tm2_index_m:tm2_index_p], color="none", hatch = used_hatch, edgecolor="k")
	
	ax.axis([-1.2, 1.2, None, None])

	ax.set_xlabel("Time (s)")
	ax.set_ylabel("Detector response, $S(t)$")
	#x_ticks = ["$-t_m$", "$-t_m/2$", "$0$", "$t_m/2$", "$t_m$"]
	x_ticks = ["$-\Delta t$", "$-\Delta t/2$", "$0$", "$\Delta t/2$", "$\Delta t$"]
	x_tick_locs = [-1, -0.5, 0, 0.5, 1.0]
	ax.set_xticklabels(x_ticks)
	ax.set_yticklabels([])
	ax.xaxis.set_major_locator(FixedLocator(x_tick_locs))
	rec_best = pl.Rectangle((0, 0), 1, 1, fc="white", hatch = used_hatch, alpha = 1.0, fill = True)
	rec_worst = pl.Rectangle((0, 0), 1, 1, fc="gray", alpha = 1.0)
	ax.legend([rec_best, rec_worst], ["$S_\mathrm{best}$", "$S_\mathrm{worst}$"])
	#ax.legend([rec_best, rec_worst], ["Best case", "Worst case"])
	
	
	pl.savefig("TimeDesc.pdf")
	#pl.show()
	return 0
