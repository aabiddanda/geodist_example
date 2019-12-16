import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib as mpl

import matplotlib.font_manager as fm
fm.findSystemFonts(fontpaths=['/home/abiddanda/.fonts/'], fontext='ttf')

mpl.rcParams['font.family'] = "Arial"
mpl.rcParams['font.sans-serif'] = "Arial"
plt.rc('font', family="Arial")
mpl.font_manager._rebuild()

# import click

# Hex-Colors for GeoDist Mappings
# Blues_HSL = [188,45,81]
# Blues_Hex [#B9DFE4, #061784, #FFFFFF]

# cmap_test = mpl.colors.LinearSegmentedColormap.from_list('Custom Blues', ['#B9DFE4', '#061784', '#FFFFFF'], 3)

# # define the bins and normalize
# bounds = np.linspace(0, 20, 21)
# norm = mpl.colors.BoundaryNorm(bounds, cmap.N)


class GeoDistPlot:

  def __init__(self):
    """ Initializing plotting object """
    # Defining general parameters
    self.fontsize = 16
    self.title_fontsize = 12
    self.h_orient ='center'
    self.v_orient = 'center'
    self.bar_color = 'gray'
    self.border_color = 'gray'
    self.line_weight = 0.5
    self.alpha = 0.5
    self.x_lbl_fontsize = 16
    self.y_lbl_fontsize = 12
    # Original Data for GeoDist
    self.orig_geodist = None
    self.orig_ngeodist = None
    self.orig_fgeodist = None
    self.orig_npops = None
    self.orig_ncat = None
    # Actual plotting data for GeoDist
    self.geodist = None
    self.ngeodist = None
    self.fgeodist = None
    self.npops = None
    self.poplist = None
    self.ncat = None
    # Colormap parameters
    self.colors = None
    self.str_labels = None
    self.lbl_colors = None

  def __str__(self):
    """ print all active parameters for plotting """
    test_str = 'Fontsize : %d\n' % self.fontsize
    test_str += 'Bar Color : %s\n' % self.bar_color
    test_str += 'Border Color : %s\n' % self.border_color
    test_str += 'Line Weight : %0.2f\n' % self.line_weight
    test_str += 'Alpha : %0.2f\n' % self.alpha
    # Printing data level parameters
    if self.ngeodist is not None:
      test_str += 'Number of SNPS : %d\n' % np.sum(self.ngeodist)
      test_str += 'Number of Populations : %d\n' % self.npops
      test_str += 'Number of Categories : %d\n' % self.ncat
    else:
      test_str += 'Number of SNPS : 0\n'
      test_str += 'Number of Populations : NA\n'
      test_str += 'Number of Categories : NA\n'

    return(test_str)

  def _add_text_data(self, inputfile):
    """ Add/replace data for a GeoDistPlot object """
    df = np.loadtxt(inputfile, dtype=str)
    geodist = df[:,0]
    ngeodist = df[:,1].astype(int)
    assert(geodist.size == ngeodist.size)
    npops = np.array([len(x) for x in geodist])
    ncat = np.array([max(list(x)) for x in geodist], dtype=int)
    assert(len(np.unique(npops)) <= 1)
    self.orig_npops = npops[0]
    self.orig_geodist = geodist
    self.orig_ngeodist = ngeodist
    self.orig_ncat = max(ncat) + 1
    self.orig_fgeodist = (self.orig_ngeodist / np.sum(self.orig_ngeodist))
    # Setting all the plotting variables
    self.npops = npops[0]
    self.geodist = geodist
    self.ngeodist = ngeodist
    self.ncat = max(ncat) + 1
    self.fgeodist = self.orig_fgeodist

  def _add_data_jsfs(self, jsfs):
    """ Add data from a joint SFS via numpy   """
    npops = len(jsfs.shape)
    ncats = np.array(jsfs.shape)
    # Assert that all of the populations have the same categories
    assert(np.all(ncats == ncats[0]))
    # iterate through them all...
    geo_codes = []
    ngeo_codes = []
    for i,x in np.ndenumerate(jsfs):
      # generate category
      cat = ''.join([str(v) for v in list(i)])
      geo_codes.append(cat)
      ngeo_codes.append(x)
    geodist = np.array(geo_codes)
    ngeodist = np.array(ngeo_codes)
    assert(geodist.size == ngeodist.size)
    self.orig_npops = npops
    self.orig_geodist = geodist
    self.orig_ngeodist = ngeodist
    self.orig_ncat = ncats[0] + 1
    self.orig_fgeodist = (self.orig_ngeodist / np.sum(self.orig_ngeodist))
    # Setting all the plotting variables
    self.npops = npops
    self.geodist = geodist
    self.ngeodist = ngeodist
    self.ncat = ncats[0] + 1
    self.fgeodist = self.orig_fgeodist

  def _filter_data(self, max_freq=0.01, rare=False):
    """ Filter geodist data for easier plotting """
    assert(self.orig_geodist is not None)
    assert(self.orig_ngeodist is not None)
    assert(self.orig_fgeodist is not None)
    if rare:
      idx = np.where(self.orig_fgeodist < max_freq)[0]
    else:
      idx = np.where(self.orig_fgeodist >= max_freq)[0]
    assert(idx[0].size != 0)
    # TODO : have this create a new variable and not reset
    self.ngeodist = self.orig_ngeodist[idx]
    self.geodist = self.orig_geodist[idx]
    self.fgeodist = self.ngeodist / np.sum(self.ngeodist)
    # sort to try to order based on frequency (high to low)
    sorted_idx = np.argsort(self.fgeodist)[::-1]
    self.ngeodist = self.ngeodist[sorted_idx]
    self.geodist = self.geodist[sorted_idx]
    self.fgeodist = self.fgeodist[sorted_idx]

  def _add_cmap(self, base_cmap='Blues',
                str_labels=['U', 'R', 'C'],
                lbl_colors=['black', 'black', 'white']):
    """ Create a colormap object to use """
    assert(self.ncat is not None)
    assert(len(str_labels) == len(lbl_colors))
    assert(len(str_labels) == self.ncat)
    # Generating a discrete mapping
    base = plt.cm.get_cmap(base_cmap)
    color_list = base(np.linspace(0, 1, self.ncat))
    cmap_name = base.name + str(self.ncat)
    test_cmap = base.from_list(cmap_name, color_list, self.ncat)
    # normalizing to 1
    norm = mpl.colors.Normalize(vmin=0, vmax=self.ncat)
    colors = [mpl.colors.rgb2hex(test_cmap(norm(i))) for i in range(self.ncat)]
    self.colors = colors
    self.str_labels = str_labels
    self.lbl_colors = lbl_colors

  def _set_colors(self, colors):
    """ Method to add custom hex colors during plotting """
    assert(self.ncat is not None)
    assert(self.colors is not None)
    assert(len(colors) == len(self.colors))
    self.colors = colors


  def _add_poplabels(self, popfile):
    """ Adding population labels here """
    assert(self.geodist is not None)
    assert(self.ngeodist is not None)
    assert(self.npops is not None)
    assert(self.ncat is not None)
    # Reading the appropriate population file
    pops = np.loadtxt(popfile, dtype=str)
    assert(pops.size == self.npops)
    self.poplist = pops

  def _add_poplabels_manual(self, poplabels):
    """ Adding population labels here """
    assert(self.geodist is not None)
    assert(self.ngeodist is not None)
    assert(self.npops is not None)
    assert(self.ncat is not None)
    assert(poplabels.size == self.npops)
    self.poplist = poplabels

  def plot_geodist(self, ax):
    """ Final function to call to generate a geodist plot on a particular axis """
    # Starting assertions to make sure we can call this
    assert(self.geodist is not None)
    assert(self.ngeodist is not None)
    assert(self.npops is not None)
    assert(self.poplist is not None)
    assert(self.ncat is not None)
    assert(self.colors is not None)
    assert(self.str_labels is not None)
    assert(self.lbl_colors is not None)
    assert(self.poplist is not None)
    # Setting up the codes here
    x_limits = ax.get_xlim()
    xbar_pts = np.linspace(x_limits[0], x_limits[1], num=self.npops + 1)
    xpts_shifted = (xbar_pts[1:] + xbar_pts[:-1])/2.0
    ax.set_xticks(xpts_shifted)
    # setting up the vertical lines
    for x in xbar_pts:
      ax.axvline(x=x, color=self.bar_color, lw=self.line_weight, alpha=self.alpha);

    # changing the border color
    for spine in ax.spines.values():
      spine.set_edgecolor(self.border_color);

    # Plotting the horizontal bars in this case
    ylims = ax.get_ylim()
    y_pts = np.cumsum(self.fgeodist) / (ylims[1] - ylims[0])
    n_y = y_pts.size
    cur_y = ylims[0]
    nsnps = np.sum(self.orig_ngeodist)
    for i in range(n_y):
      y = y_pts[i]
      y_dist = y - cur_y
      ax.axhline(y = y, color=self.bar_color, lw=self.line_weight, alpha=self.alpha);
      cur_code = list(self.geodist[i])
      for j in range(self.npops):
        # Defining the current category
        cur_cat = int(cur_code[j])
        # Drawing in the rectangle
        cur_xy = xbar_pts[j], cur_y
        rect = patches.Rectangle(xy = cur_xy,
                                 width=(xbar_pts[j+1] - cur_xy[0]),
                                 height=y_dist,
                                 facecolor=self.colors[cur_cat]);
        ax.add_patch(rect);
        # Drawing in the text
        fontscale = min(1.25, self.fontsize*y_dist/(ylims[1] - ylims[0]))
        ax.text(x=xpts_shifted[j],
                y=(cur_y + y_dist/2.0),
                ha=self.h_orient, va=self.v_orient,
                s=self.str_labels[cur_cat],
                color=self.lbl_colors[cur_cat],
                fontsize=self.fontsize*fontscale);
      cur_y = y
    ax.set_ylim(ylims)
    return(ax, nsnps, y_pts)
  
  def plot_percentages(self, ax):
    """ Generates a plot with the percentages  """
    ns = self.ngeodist
    fracs = ns/np.sum(ns)
    cum_frac = np.cumsum(fracs)
    # Setting the border here
    for spine in ax.spines.values():
      spine.set_edgecolor(self.border_color);
    prev = 0.0
    for i in range(cum_frac.size):
      ax.axhline(cum_frac[i], lw=0.5, color=self.bar_color)
      # Get the midpoint 
      ydist = (cum_frac[i] - prev)/2.
      fontscale = min(1.5, self.fontsize*fracs[i])
      ax.text(x=0.01, y=prev+ydist, s = '%d (%0.2f%%)' % (ns[i], fracs[i]*100), fontsize=self.fontsize*fontscale)
      prev = cum_frac[i]
    return(ax)
      