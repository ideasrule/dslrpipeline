#!/usr/bin/python -u

from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot
from scipy import interpolate, array
from subprocess import Popen, PIPE
from os.path import basename, splitext
from os.path import join as join_paths
from optparse import OptionParser
from exceptions import Exception
import numpy as np

def read_spfit(fname, columns) :
    """Read the columns with names listed in columns from the given file."""

    f=open(fname, 'r')
    result=[[] for c in columns]
    lineno=1
    for l in f :
        if l[0]=='#' : 
            colnames=l.strip('#').strip().split()
            colnums=map(colnames.index, columns)
            continue
        entries=l.split()
        for col, colnum in zip(result, colnums) :
            try :
                col.append(eval(entries[colnum]))
            except Exception, ex :
                print 'failed to parse column %d (%s) in %s:%d'%(
                    colnum, colnames[colnum], fname, lineno)
                raise
        lineno+=1
    return result

def spline_smooth(x, y, z, weights=None, smooth=None) :
    """ Returns a smoothed version of z(x,y) evaluated at the input points.
    
    Arguments:
    x -- The x coordinates of the points at which the function is known
    y -- The y coordinates of the points at which the function is known
    z -- The noisy values of the function at the points (x[i], y[i])
    weights -- weights to give to the point at (x[i], y[i]) when fitting the
               smoothing spline, if None, equal weight are used for all
               points.
    smooth -- The amount of smoothing to apply. See the s argument to
              scipy.interpolate.bisplrep.
    
    Return value:
    A list of exactly the same length as x, y, z and weights containing the
    values of the smoothed function at the points (x[i], y[i])."""

    tck=interpolate.bisplrep(x, y, z, w=weights, s=smooth)
    return map(lambda xi, yi: interpolate.bisplev(xi, yi, tck), x, y)

def derive_poly_trans(x, y, s, d, k, weights, trans_fname, order=5,
                      iterations=20, sigma=3) :
    """ Compute and store in a file the best fit S, D, K polynomials.

    Arguments:
    x -- The x coordinates of the stars for which S, D, K was fit.
    y -- The y coordinates of the stars for which S, D, K was fit.
    s -- The best fit S value for each star.
    d -- The best fit D value for each star.
    k -- The best fit K value for each star.
    weights -- The weight to give each star in the fit.
    trans_fname -- The filename where to store the deriva transformation.
    order -- The polynomial order of the fit.
    iterations -- How many outlier rejection/re-fitting iterations to do.
    sigma -- How many sigma away from the fit to reject sources during the
             iterations.

    Always returns None. """

    command=['grtrans',  '--input', '-', '--col-xy', '1,2',
             '--col-fit', '3,4,5', '--col-weight', '6', '--order',
             str(order), '--iterations', str(iterations), '--sigma',
             str(sigma), '--offset', '2048,2048', '--scale', '2048',
             '--output-transformation', trans_fname]
    grtrans=Popen(command, stdin=PIPE)
    if weights is None : weights=[1 for xi in x]
    for xysdkw in zip(x, y, s, d, k, weights) :
        grtrans.stdin.write(' '.join(map(repr, xysdkw))+'\n')
    grtrans.stdin.close()
    grtrans.wait()

def poly_smooth(x, y, trans) :
    """ Returns smoothed S(x,y), D(x,y) and K(x, y) at the input points.
    
    Arguments:
    x -- The x coordinates of the points at which the function is known
    y -- The y coordinates of the points at which the function is known
    trans -- The transformation file containing the smoothed polynomial
             coefficients.
   
    Return value:
    Three lists of exactly the same length as x and y containing the values
    of the smoothed S, D and K at the points (x[i], y[i])."""

    command=['grtrans', '--input', '-', '--col-xy', '1,2', '--col-out',
             '3,4,5', '--input-transformation', trans]
    grtrans=Popen(command, stdin=PIPE, stdout=PIPE)
    smooth_s, smooth_d, smooth_k=[], [], []
    for xy in zip(x, y) :
        grtrans.stdin.write(' '.join(map(repr, xy))+'\n')
        line=grtrans.stdout.readline()
        smooth_s.append(eval(line.split()[2]))
        smooth_d.append(eval(line.split()[3]))
        smooth_k.append(eval(line.split()[4]))
    grtrans.stdin.close()
    grtrans.stdout.close()
    grtrans.wait()
    return smooth_s, smooth_d, smooth_k

def plot_residuals(x, y, z, smooth_z, weights, plot_label, pdf,
                   min_symbol_size=3, max_symbol_size=60, alpha=0.5,
                   imwidth=2601, imheight=1734) :
    """ Generates a plot of the residuals of the smoothing function.

    Arguments:
    x -- The x coordinates of the points at which the function is known
    y -- The y coordinates of the points at which the function is known
    z -- The noisy values of the function at the points (x[i], y[i])
    smooth_z -- The smoothed values of the function at the points (x[i],y[i])
    filename -- The filename to save the plots to (multi page PDF format).

    Always returns None."""

    residuals=z-smooth_z
    pyplot.clf()
    pyplot.title(plot_label)
    binsize = 1./50.
    pyplot.hist(residuals, bins=100, normed=True)
    pyplot.axvline()
#    pyplot.xlim(-1,1)
    pyplot.xlabel('Smoothing residuals')
    pyplot.ylabel('Frequency')
    pdf.savefig()

    pyplot.clf()
    pyplot.title(plot_label)
    pyplot.xlabel('Distance from center [pixels]')
    pyplot.ylabel('Smoothing residuals')
    min_weight=min(weights)
    sizes=(max_symbol_size*(array(weights)-min_weight)/
           (max(weights)-min_weight))+min_symbol_size
    pyplot.scatter((x**2+y**2)**0.5, residuals, s=sizes, marker='o', c=sizes,
                   alpha=alpha, linewidths=(0,))
    pyplot.axhline()
    pyplot.ylim(-1,1)
    pdf.savefig()

    nstripes=10
    for stripe in range(nstripes) :
        xmin=float(imwidth)*stripe/nstripes
        xmax=float(imwidth)*(stripe+1)/nstripes
        ystripe, res_stripe, weight_stripe=[], [], []
        for i in range(len(x)) :
            if x[i]>=xmin and x[i]<xmax :
                ystripe.append(y[i])
                res_stripe.append(residuals[i])
                weight_stripe.append(weights[i])
        min_weight=min(weight_stripe)
        sizes=(max_symbol_size*(array(weight_stripe)-min_weight)/
                        (max(weight_stripe)-min_weight))+min_symbol_size
        pyplot.clf()
        pyplot.scatter(ystripe, res_stripe, s=sizes, marker='o', c=sizes,
                       alpha=alpha, linewidths=(0,))
        pyplot.axhline()
        pyplot.ylim(-1,1)
        pyplot.xlabel('y [pix]')
        pyplot.ylabel('Smoothing residuals for %g<x<%g'%(xmin, xmax))
        pdf.savefig()

    for stripe in range(nstripes) :
        ymin=float(imheight)*stripe/nstripes
        ymax=float(imheight)*(stripe+1)/nstripes
        xstripe, res_stripe, weight_stripe=[], [], []
        for i in range(len(y)) :
            if y[i]>=ymin and y[i]<ymax :
                xstripe.append(x[i])
                res_stripe.append(residuals[i])
                weight_stripe.append(weights[i])
        min_weight=min(weight_stripe)
        sizes=(max_symbol_size*(array(weight_stripe)-min_weight)/
                        (max(weight_stripe)-min_weight))+min_symbol_size
        pyplot.clf()
        pyplot.scatter(xstripe, res_stripe, s=sizes, marker='o', c=sizes,
                       alpha=alpha, linewidths=(0,))
        pyplot.axhline()
        pyplot.ylim(-1,1)
        pyplot.xlabel('x [pix]')
        pyplot.ylabel('Smoothing residuals for %g<y<%g'%(ymin, ymax))
        pdf.savefig()

    pyplot.xlim(0, 4096)
    pyplot.ylim(0, 4096)

def parse_command_line() :
    """Returns a structre with all the information from the command line."""

    parser=OptionParser(usage='%prog <FitPSF_file1> [<FitPSF_file2> [...]] '
    '[options]', description='Compute and store in a file the best fit S, '
    'D, K polynomials to a field of stars for whith FitPSF was run with '
    '--order=-1.')
    parser.add_option('-o', '--output', type='string', default='.',
                      help='The name of directory to place the best fit '
                      'transformation files. The files are named exactly '
                      'like the input frames, but with extension sdkfit. '
                      'Default: %default.')
    parser.add_option('-d', '--diagnostic', type='string', default='',
                      help='If this option is passed, a PDF file with the '
                      'given name is generated containing many pages of '
                      'diagnostic plots of the derived fit.')
    try :
        (options,args)=parser.parse_args()
    except Error.CommandLine, error:
        parser.print_help()
        print error
        exit(1)
    if len(args)==0 :
        parser.print_help()
        print "Expected at least one input file."
        exit(1)
    options.input_files=args
    return options

if __name__== '__main__' :
    options=parse_command_line()
    if options.diagnostic : pdf=PdfPages(options.diagnostic)
    for spsdk in options.input_files :
        x, y, s, d, k, sn=map(array, 
                              read_spfit(spsdk,
                                         ['x', 'y', 'S', 'D', 'K', 'S/N']))
        weights=array(sn)**2
        trans_fname=join_paths(options.output, 
                               splitext(basename(spsdk))[0]+'.sdkfit')
        derive_poly_trans(x, y, s, d, k, weights, trans_fname)
        if options.diagnostic :
            smooth_s_poly, smooth_d_poly, smooth_k_poly=poly_smooth(
                x, y, trans_fname)
            print 'len(s)=%d, len(smooth_s_poly)=%d'%(len(s),
                                                      len(smooth_s_poly))
            print 'len(d)=%d, len(smooth_d_poly)=%d'%(len(d),
                                                      len(smooth_d_poly))
            print 'len(k)=%d, len(smooth_k_poly)=%d'%(len(k),
                                                      len(smooth_k_poly))
            plot_residuals(x, y, s, smooth_s_poly, weights, spsdk+' S', pdf)
            plot_residuals(x, y, d, smooth_d_poly, weights, spsdk+' D', pdf)
            plot_residuals(x, y, k, smooth_k_poly, weights, spsdk+' K', pdf)

    """smooth_s_spline=spline_smooth(x, y, s, weights=weights, smooth=1000)
    smooth_d_spline=spline_smooth(x, y, d, weights=weights, smooth=1000)
    smooth_k_spline=spline_smooth(x, y, k, weights=weights, smooth=1000)"""


    """plot_residuals(x, y, s, smooth_s_spline, weights, 'S', pdf)
    plot_residuals(x, y, d, smooth_d_spline, weights, 'D', pdf)
    plot_residuals(x, y, k, smooth_k_spline, weights, 'K', pdf)"""

    if options.diagnostic : pdf.close()
