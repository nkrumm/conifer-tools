import pandas
import numpy as np
import sys
import matplotlib
import matplotlib.pyplot as plt
import pylab as P
import rpkm
from conifertools import CallTable
import operator
import locale

def _contiguous_regions(condition):
    """Finds contiguous True regions of the boolean array "condition". Returns
    a 2D array where the first column is the start index of the region and the
    second column is the end index."""
    
    # Find the indicies of changes in "condition"
    d = np.diff(condition)
    idx, = d.nonzero() 
    
    # We need to start things after the change in "condition". Therefore, 
    # we'll shift the index by 1 to the right.
    idx += 1
    
    if condition[0]:
        # If the start of condition is True prepend a 0
        idx = np.r_[0, idx]
    
    if condition[-1]:
        # If the end of condition is True, append the length of the array
        idx = np.r_[idx, condition.size - 1]
    
    # Reshape the result into two columns
    idx.shape = (-1,2)
    return idx

class ConiferPlotTrack(object):
    def __init__(self, plotter, data_in, name, collapse=True, position='auto', data_field=None, collapsed_linespacing=0.125, **args):
        self.name = name
        self.p = plotter
        self.position = position
        self.collapsed = collapse
        # Default style
        self.style = {"color":'r',"linewidth":5, "linestyle":'-',"alpha":0.5,"solid_capstyle":'butt'}
        self.collapsed_linespacing = collapsed_linespacing
        for key in args:
            if key in ["color", "linewidth", "linestyle","alpha"]:
                self.style[key] = args[key]

        self.fields = ["chromosome", "start", "stop"]
        if data_field is not None:
            self.fields.extend(data_field)
            self.data_field = data_field

        self._data = pandas.DataFrame()
        
        if isinstance(data_in, CallTable):
            self._import(data_in.calls[self.fields], type='calltable')
        else:
            try:
                # read in bedfile
                bedfile = pandas.read_csv(data_in,sep="\t",index_col=False,names=["chromosome","start","stop"])
            except IOError:
                # Cannot read file
                raise IOError("Cannot read bed file %s" % data_in)
            self._import(bedfile, type='bed')
 
    def _import(self, bedfile, type):
        bedfile_chrs = set(bedfile["chromosome"].values)
        for probe_tbl in self.p.r.h5file.root.probes:
            probe_array = pandas.DataFrame(probe_tbl.read()[["start","stop"]])
            probe_array_chr = int(probe_tbl.title.replace("chr",""))
            probe_array["chr"] = probe_array_chr

            # semi-hack to convert chrX/Y to chr23/chr24
            if type == 'calltable':
                chr_string = probe_array_chr
            elif type == 'bed':
                if probe_array_chr not in bedfile_chrs:
                    if probe_array_chr == 23:
                        chr_string = "chrX"
                    elif probe_array_chr == 24:
                        chr_string = "chrY"
                    else:
                        chr_string = probe_tbl.title
                else:
                    chr_string = probe_tbl.title
            else:
                print "Indicate either calltable or bedfile for type!"

            bed_array = bedfile[bedfile["chromosome"] == chr_string]
            if len(bed_array) == 0:
                continue
            # this line returns if the intervals of the probe_array are contained or overlapping the bed_array.
            # It works by using the trick that the searchsorted function will return differing indices if a start or stop
            # from bed_array falls between one of the probes, in that case the searchsorted result will not be the same.    
            if self.collapsed:
                probe_array["feature"] = np.searchsorted(np.sort(bed_array["start"].values), probe_array["stop"].values) != np.searchsorted(np.sort(bed_array["stop"].values), probe_array["start"].values)
                self._data = self._data.append(probe_array,ignore_index=True)
            else:
                bed_array = bed_array.sort(["start","stop"])
                bed_array["chromosome"] = probe_array_chr
                bed_array["exon_start"] = np.searchsorted(probe_array["start"].values, bed_array["start"].values)
                bed_array["exon_stop"] = np.searchsorted(probe_array["stop"].values, bed_array["stop"].values)
                self._data = self._data.append(bed_array, ignore_index=True)

    def plot(self, ax, position, chr, data_range):
        
        if self.collapsed:
            for start,stop in _contiguous_regions(self._data[self._data.chr==chr]["feature"].values[data_range[0]:data_range[-1]]):
                _ = ax.add_line(matplotlib.lines.Line2D([start-0.5,stop-0.5],[position,position],**self.style))
        else:
            range_mask = (self._data["chromosome"] == chr) & (self._data.exon_stop>=data_range[0]) & (self._data.exon_start<=data_range[-1])
            curr_pos = position
            for ix, row in self._data[range_mask].iterrows():
                start = row["exon_start"] - data_range[0]
                stop = row["exon_stop"] - data_range[0]
                if start == stop:
                    continue
                curr_pos -= self.collapsed_linespacing
                _ = ax.add_line(matplotlib.lines.Line2D([start,stop],[curr_pos,curr_pos],**self.style))


class ConiferPlotter():
    def __init__(self, conifer_file):
        try:
            _ = locale.setlocale(locale.LC_ALL, '')
        except:
            print "Warning: Could not set locale"

        if conifer_file is not None:
            self.h5file = conifer_file
            # TODO test if valid file
            try:
                self.r = rpkm.rpkm_reader(self.h5file)
            except:
                print "Error in reading or opening HDF5 analysis file"
                return None
            # get contig names
            self.contigs = self.r.getContigList()

            # get sample names
            self.samples = self.r.getSampleList()

            # set up empty track list
            self.tracks = []
        else:
            pass

    def _plotRawData(self, axis, rpkm_data, color='r', linewidth=0.9,label=None):
        zero_stack = np.zeros(len(rpkm_data))
        positions = np.repeat(np.arange(0, len(rpkm_data)), 3)
        logr = np.vstack([zero_stack, rpkm_data.flatten(), zero_stack]).transpose().ravel()
        axis.plot(positions, logr, color=color, marker=None, linewidth=linewidth,label=label)

    def _plotGenes2(self, axis, rpkm_data, levels=5,y_pos=-2, data_range=None,text_pos='right',line_spacing=0.1,text_offset=0.25,linewidth=5,fontsize=6, **kwargs):
        def drawGene(start, stop, axis, prev_gene, prev_gene_strand, prev_gene_isPPG, counter, y_pos, line_spacing,text_offset,linewidth,fontsize):
            
            if prev_gene_isPPG:
                gene_color = (0,0,0,0.3)
            else:
                gene_color = (102/255.,33/255.,168/255.,0.6)
            
            dir_color = (1,1,1,0.5)
            
            gene_y_pos = y_pos - (counter * line_spacing)
            
            x_points = np.arange(start-0.5,stop-0.5,3)
            y_points = np.repeat(gene_y_pos,len(x_points))
            
            _ = axis.add_line(matplotlib.lines.Line2D([start-0.5,stop-0.5],[gene_y_pos,gene_y_pos] ,color=gene_color,linewidth=linewidth,linestyle='-',solid_capstyle='butt'))
            
            if prev_gene_strand == "+":
                _ = axis.add_line(matplotlib.lines.Line2D(x_points[1:], y_points[1:],color=dir_color,linewidth=4,linestyle='>',marker=">",markersize=5.5,markeredgewidth=0,markerfacecolor=dir_color,alpha=0.5, solid_capstyle='butt'))
            elif prev_gene_strand == "-":
                _ = axis.add_line(matplotlib.lines.Line2D(x_points[1:], y_points[1:],color=dir_color,linewidth=4,linestyle='<',marker="<",markersize=5.5,markeredgewidth=0,markerfacecolor=dir_color,alpha=0.5, solid_capstyle='butt'))
            else:
                _ = axis.add_line(matplotlib.lines.Line2D(x_points[1:], y_points[1:],color=dir_color,linewidth=4,linestyle='-',alpha=0.5, solid_capstyle='butt'))
            
            if prev_gene != 'None':
                _ = axis.text(stop+text_offset,gene_y_pos, prev_gene, ha='left',va='center',fontsize=fontsize)
        
        
        if data_range is not None:
            exon_set = rpkm_data.exons[data_range]
        else:
            exon_set = rpkm_data.exons
        
        if "isPPG" in exon_set[0]:
            prev_gene_isPPG =  exon_set[0]["isPPG"]
        else:
            prev_gene_isPPG = False
        
        if "strand" in exon_set[0]:
            prev_gene_strand =  exon_set[0]["strand"]
        else:
            prev_gene_strand = ""
        
        prev_gene = exon_set[0]["name"]
        prev_stop = -1
        i = -1
        counter = -1
        for exon in exon_set:
            i += 1
            gene = exon["name"]
            if exon["name"] == prev_gene:
                stop = i
                continue
            
            start = prev_stop 
            stop = i
            
            drawGene(start, stop, axis, prev_gene, prev_gene_strand, prev_gene_isPPG, counter, y_pos,  line_spacing,text_offset,linewidth,fontsize)
            
            counter +=1
            prev_gene = exon["name"]
            try:
                prev_gene_isPPG = exon["isPPG"]
            except IndexError:
                prev_gene_isPPG = False
            try:
                prev_gene_strand = exon["strand"]
            except IndexError:
                prev_gene_strand = ""
            
            prev_stop = stop
            if counter > levels:
                counter = 0
        
        drawGene(prev_stop, stop, axis, prev_gene, prev_gene_strand, prev_gene_isPPG, counter, y_pos,  line_spacing,text_offset,linewidth,fontsize)

    def _plotSegDups(self, axis, rpkm_data, y_pos=2.5,data_range=None):
        if data_range is not None:
            segdups = rpkm_data.exons[data_range]["isSegDup"].flatten()
        else:
            segdups = rpkm_data.exons["isSegDup"].flatten()
        
        for start,stop in _contiguous_regions(segdups):
            _ = axis.add_line(matplotlib.lines.Line2D([start-0.5,stop-0.5],[y_pos,y_pos],color='y',linewidth=5,linestyle='-',alpha=0.5,solid_capstyle='butt'))
    
    def _plotGenomicCoords(self, plt, rpkm_data, interval=50,fontsize=10,rotation=0,data_range=None):
        if data_range is not None:
            exon_set = rpkm_data.exons[data_range]
        else:
            exon_set = rpkm_data.exons
        genomic_coords = np.array(map(operator.itemgetter("start"),exon_set))
        ticks = range(0,len(exon_set)+interval,interval)
        ticks[-1] -= 1 # the last tick is going to be off the chart, so we estimate it as the second to last genomic coord.
        labels = [locale.format("%d", genomic_coords[i], grouping=True) for i in ticks if i < len(genomic_coords)]
        if rotation != 0:
            ha = "right"
        else:
            ha = "center"
        _ = plt.xticks(ticks,labels,fontsize=fontsize,rotation=rotation,ha=ha)

    def _chrInt2Str(self, chromosome_int):
        if int(chromosome_int) == 23:
            return 'chrX'
        elif int(chromosome_int) == 24:
            return 'chrY' 
        else:
            return 'chr' + str(chromosome_int)

    def _chrStr2Int(self, chromosome_str):
        chr = chromosome_str.replace('chr','')
        if chr == 'X':
            return 23
        elif chr == 'Y':
            return 24 
        else:
            return int(chr)
    
    def add_track(self, track):
        if track.position == 'auto':
            if len(self.tracks) > 0:
                pos = min(map(operator.itemgetter("position"), self.tracks)) - 0.075
            else:
                pos = 2.2
        else:
            pos = track.position
        
        self.tracks.append({"track":track,
                            "name": track.name,
                            "position": pos})

    def basicPlot(self, call, window=50, outdir=None):
        chromosome = int(call["chromosome"])
        start = int(call["start"])
        stop = int(call["stop"])
        sampleID = str(call["sampleID"])
        exons = self.r.getExonIDs(chromosome,start,stop)
        
        window_start = max(exons[0]-window,0)
        window_stop = exons[-1]+window
        
        data = self.r.getExonValuesByExons(chromosome,window_start, window_stop,sampleList=[sampleID])
        #_ = data.smooth()
        
        fig, ax = plt.subplots(figsize=(10,5))
        #ax.plot(data.rpkm, linewidth = 0.3, c='k')
        #ax.plot(data.getSample([sampleID]), linewidth = 1, c='r', label = sampleID)

        self._plotRawData(ax, data.getSample([sampleID]),color='r', label=sampleID)
        
        plt.legend(prop={'size':10},frameon=False)
        
        self._plotGenes2(ax,data, fontsize=8, linewidth=6, line_spacing=0.2, levels=3)
        label_interval = (window_stop-window_start)/5
        self._plotGenomicCoords(plt,data,interval=label_interval)
        #self._plotSegDups(ax, data,y_pos=2.5)

        for t in self.tracks:
            t["track"].plot(ax, position=t["position"], chr=chromosome, data_range=xrange(window_start, window_stop))

        exon_start = np.where(data.exons["start"] == start)[0][0]
        exon_stop = np.where(data.exons["stop"] == stop)[0][0]
        _ = ax.add_line(matplotlib.lines.Line2D([exon_start,exon_stop],[2,2],color='k',lw=6,linestyle='-',alpha=1,solid_capstyle='butt'))
        
        _ = plt.xlim(0,data.shape[1])
        _ = plt.ylim(-3,3)
        plt.xlabel("Position")
        plt.ylabel("SVD-ZRPKM Values")
        
        if "cnvrID" in call:
            plt.title("%s; %s: %s - %s; cnvr: %d; frequency: %d" % (sampleID, self._chrInt2Str(chromosome),locale.format("%d",start, grouping=True),locale.format("%d",stop, grouping=True), call["cnvrID"], call["cnvr_frequency"]))
            outfile = "%d_%s_%d_%d_%s.png" % (call["cnvrID"], self._chrInt2Str(chromosome), start, stop, sampleID)
        else:
            plt.title("%s; %s: %s - %s" % (sampleID, self._chrInt2Str(chromosome),locale.format("%d",start, grouping=True),locale.format("%d",stop, grouping=True)))
            outfile = "%s_%d_%d_%s.png" % (self._chrInt2Str(chromosome), start, stop, sampleID)


        
        if outfile is not None:
            plt.savefig(outdir + "/" + outfile)
            P.close(fig)
            plt.clf()
            del fig
            del ax            
        else:
            plt.show()


def QC_Chromosome_Plot(calls, title=None, outfile=None):
    x = np.array(calls.calls["num_probes"].values)
    y = np.array(calls.calls["median_svdzrpkm"].values)
    
    fig = plt.figure(1, figsize=(9,9))
    
    from mpl_toolkits.axes_grid import make_axes_locatable
    
    axScatter = plt.subplot(111)
    divider = make_axes_locatable(axScatter)
    
    # create a new axes with a height of 1.2 inch above the axScatter
    axHistx = divider.new_vertical(1.2, pad=0.1, sharex=axScatter)
    
    # create a new axes with a width of 1.2 inch on the right side of the
    # axScatter
    axHisty = divider.new_horizontal(1.2, pad=0.1, sharey=axScatter)
    
    fig.add_axes(axHistx)
    fig.add_axes(axHisty)
    
    # make some labels invisible
    plt.setp(axHistx.get_xticklabels() + axHisty.get_yticklabels(),
             visible=False)
    
    # the scatter plot:
    #axScatter.scatter(x[x<0], y[x<0], lw=0, alpha=0.3, color="r")
    axScatter.scatter(x[x>=0], y[x>=0], lw=0, alpha=0.3, color="b")
    #axScatter.set_aspect(1.)
    axScatter.set_xscale("symlog")
    axHistx.set_xscale("symlog")
    
    axScatter.set_xlabel("Size of call (# of probes)")
    axScatter.set_ylabel("Signal Strength (Median SVD-ZRPKM)")
                         
    axHistx.hist(x, bins=np.arange(0,np.max(x)), histtype="stepfilled", lw=0,align='left')
    axHisty.hist(y, bins=200, orientation='horizontal', histtype="stepfilled", lw=0)
    
    for tl in axHistx.get_xticklabels():
        tl.set_visible(False)
    
    axHisty.set_xticklabels(["%d" % i for i in axHisty.get_xticks()], rotation=-90)
    if title is not None:
        axHistx.set_title(title)
    
    if outfile is not None:
        plt.savefig(outfile)

