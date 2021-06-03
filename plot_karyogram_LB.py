# Adopted from a script originally written by Alicia Martin
# https://github.com/armartin/ancestry_pipeline

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pylab
from matplotlib.path import Path
import matplotlib.patches as patches
import matplotlib.collections as mcol
import brewer2mpl
import os

CENTROMERE_FILE = "/datasets/cs284s-sp20-public/ps2/rfmix/centromeres_hg19.bed"
centro = open(CENTROMERE_FILE)
centromeres = {}
for line in centro:
    line = line.strip().split()
    centromeres[line[0]] = line
    
def plot_rects(anc, chr, start, stop, hap, pop_order, colors, ax):  
    centro_coords = [float(item) for item in centromeres[str(chr)]] 
    if len(centro_coords) == 3: #acrocentric chromosome
        val = 0.2
        mask = [
        (centro_coords[1]+2,chr-val), #add +/- 2 at the end of either end
        (centro_coords[2]-2,chr-val),
        (centro_coords[2]+2,chr),
        (centro_coords[2]-2,chr+val),
        (centro_coords[1]+2,chr+val),
        (centro_coords[1]-2,chr),
        (centro_coords[1]+2,chr-val)
        ]
        
        mask_codes = [
        Path.MOVETO,
        Path.LINETO,
        Path.CURVE3,
        Path.LINETO,
        Path.LINETO,
        Path.CURVE3,
        Path.LINETO,
        ]
        clip_mask = Path(vertices=mask, codes=mask_codes)
    
    else: #need to write more complicated clipping mask with centromere masked out
        mask = [
        (centro_coords[1]+2,chr-val), #add +/- 2 at the end of either end
        (centro_coords[2]-2,chr-val),
        (centro_coords[2]+2,chr+val),
        (centro_coords[3]-2,chr+val),
        (centro_coords[3]+2,chr),
        (centro_coords[3]-2,chr-val),
        (centro_coords[2]+2,chr-val),
        (centro_coords[2]-2,chr+val),
        (centro_coords[1]+2,chr+val),
        (centro_coords[1]-2,chr),
        (centro_coords[1]+2,chr-val)
        ]
        
        mask_codes = [
        Path.MOVETO,
        Path.LINETO,
        Path.LINETO,
        Path.LINETO,
        Path.CURVE3,
        Path.LINETO,
        Path.LINETO,
        Path.LINETO,
        Path.LINETO,
        Path.CURVE3,
        Path.LINETO,
        ]
        clip_mask = Path(vertices=mask, codes=mask_codes)
        
    if hap == 'A': #bed_a ancestry goes on top
        verts = [
            (float(start), chr + 0.025), #left, bottom # Lauryn Changed
            (float(start), chr + 0.5), #left, top # Lauryn Changed
            (float(stop), chr + 0.5), #right, top # Lauryn Changed
            (float(stop), chr + 0.025), #right, bottom # Lauryn Changed
            (0, 0), #ignored
        ]
    else: #bed_b ancestry goes on bottom
        verts = [
            (float(start), chr - 0.45), #left, bottom # Lauryn Changed
            (float(start), chr - 0.025), #left, top # Lauryn Changed
            (float(stop), chr - 0.025), #right, top # Lauryn Changed
            (float(stop), chr - 0.45), #right, bottom # Lauryn Changed
            (0, 0), #ignored
        ]

    codes = [
        Path.MOVETO,
        Path.LINETO,
        Path.LINETO,
        Path.LINETO,
        Path.CLOSEPOLY,
    ]
    
    clip_path = Path(verts, codes)
    if anc in pop_order:
        col=mcol.PathCollection([clip_path],facecolor=colors[pop_order.index(anc)], linewidths=0)
    else:
        col=mcol.PathCollection([clip_path],facecolor=colors[-1], linewidths=0)
    if 'clip_mask' in locals():
        col.set_clip_path(clip_mask, ax.transData)
    ax.add_collection(col)
    return ax

def splitstr(option, opt, value, parser):
    return(setattr(parser.values, option.dest, value.split(',')))

def hex_to_rgb(value):
    value = value.lstrip('#')
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))

def PlotKaryogram(bed_a, bed_b, ind, pop_order, outprefix, xaxis_length, yaxis_range):
    """ Plot admixture karyogram
    
    Parameters
    ----------
    bed_a : str
       Input file of segments for the first chromosome copy
    bed_b : str
       Input file of segments for the first chromosome copy
    ind : str
       Sample ID of the individual you want to plot
    pop_order : list of str
       List of population labels [pop1, pop2...]
       Corresponds to same order used to label ancestral groups with RFMix.
    """
    # Load bed files
    bed_a = open(bed_a)
    bed_b = open(bed_b)

    # define plotting space

    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111)
    #ax.set_xlim(-5,300)
    #ax.set_ylim(23,0)
    ax.set_xlim(-5, xaxis_length) # Lauryn Changed
    ax.set_ylim(yaxis_range[0], yaxis_range[1]) # Lauryn Changed
    #plt.xlabel('Genetic position (cM)')
    
    plt.xlabel('Genetic position (cM)', fontsize=16) # Lauryn Changed
    plt.ylabel('Chromosome', fontsize=16) 
    plt.title(ind, fontsize=16)
    plt.yticks(range(yaxis_range[0], yaxis_range[1])) # Lauryn Changed
    plt.xticks(fontsize=20)
    
    ax.get_yaxis().set_visible(False) # Lauryn Added
    # TODO: Hide yaxis line?
    
    # Define colors
    bmap = brewer2mpl.get_map('Set1', 'qualitative', 4)
    colors=bmap.mpl_colors
    colors.append((0,0,0))
    
    # Define centromeres
    centro = open(CENTROMERE_FILE)
    centromeres = {}
    for line in centro:
        line = line.strip().split()
        centromeres[line[0]] = line

    # Plot rectangles
    for line in bed_a:
        line = line.strip().split()
        ax = plot_rects(line[3], int(line[0]), line[4], line[5], 'A', pop_order, colors, ax)
    for line in bed_b:
        line = line.strip().split()
        ax = plot_rects(line[3], int(line[0]), line[4], line[5], 'B', pop_order, colors, ax)

    # Write a legend
    p = []
    for i in range(len(pop_order)):
        p.append(plt.Rectangle((0, 0), 1, 1, color=colors[i]))
    p.append(plt.Rectangle((0, 0), 1, 1, color='k'))
    labs = list(pop_order)
    labs.append('UNK')
    leg = ax.legend(p, labs, loc=4, fancybox=True, prop={"size": 20})
    leg.get_frame().set_alpha(0)

    # Get rid of annoying plot features
    spines_to_remove = ['top', 'right']
    for spine in spines_to_remove:
        ax.spines[spine].set_visible(False)
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('none')
    
    fig.savefig(outprefix)

    # Return the axis
    return ax
