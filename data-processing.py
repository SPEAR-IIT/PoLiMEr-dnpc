import os
import numpy as np
import pandas as pd
import platform
import argparse
import sys
import datetime
import matplotlib.pyplot as plt
from matplotlib.pyplot import text
from adjustText import adjust_text
import custom_mpl as tx

MONITOR= "PoLiMEr" #"iMon"

tx.latexify(columns=1)

parser = argparse.ArgumentParser()
parser.add_argument('-p','--prefix-path', help='Base path containing the '+MONITOR+' output files. Default: current working directory.', default='.', required=True)
parser.add_argument('-x','--'+MONITOR+'-prefix', help='The value of the "'+MONITOR+'_PREFIX" environment variable used.', default='', required=False)
parser.add_argument('-o','--output-path', help='Location at which plots should be stored. Default: current working directory.', default='.', required=False)
parser.add_argument('-f', '--output-prefix', help='Prefix to be appended to every output file.', default='', required=False)
parser.add_argument('-t','--num-threads', help='Number of threads used for the run. Use only for marking output plots.', default='', required=False)
parser.add_argument('-rpn', '--ranks-per-node', help='Number of reanks per node used for the run. Use only for marking output plots.', default='', required=False)
parser.add_argument('-j', '--threads-per-core', help='Number of threads per core used for the run. Use only for marking output plots.', default='', required=False)

args = parser.parse_args()

slash = '/' #separator for paths, platform dependent

print("Determining platform:", platform.system())
if platform.system() == 'Windows':
    slash = '\\'
    
figformat = '.pdf'
rounding = 1

def collect_files(path):
    print("Collecting files in", path)
    ''' Filenames should be of the following format:
        <iMon_PREFIX>_iMon_<node_id>_<job_id>.txt
        <iMon_PREFIX>_iMon_energy-tags_<node_id>_<job_id>.txt
        <iMon_PREFIX>_iMon_powercap-tags_<node_id>_<job_id>.txt
    '''
    nodes_per_job = {}
    jobfiles = {}
    
    jobid_index = -1
    node_index = -2
    tag_name_index = -3
    
    for root, dirs, files in os.walk(path):
        for file in files:
            if MONITOR in file:
                
                if "CLEAN" in file:
                    continue
                prefix = ""
                filenamecomponents = file.split('_')
                jobid = filenamecomponents[jobid_index].split('.')[0].strip()
                node = filenamecomponents[node_index].strip()
                filetype = filenamecomponents[tag_name_index].strip()
                prefix_end = -3
                if filetype == "energy-tags" or filetype == "powercap-tags":
                    prefix_end = -4
                if filenamecomponents[0] != MONITOR:
                    if MONITOR in filenamecomponents[0]:
                        clean = filenamecomponents[0].split(MONITOR)
                        pfxnames = []
                        for c in clean:
                            if c != MONITOR:
                                pfxnames.append(c.strip())
                        prefix = '-'.join(pfxnames)
                    else:
                        if MONITOR in filenamecomponents[prefix_end]:
                            clean = filenamecomponents[prefix_end].split(MONITOR)
                            pfxnames = []
                            for c in clean:
                                if c != MONITOR and len(c) > 0:
                                    pfxnames.append(c.strip())
                            prefix = '-'.join(filenamecomponents[:prefix_end] + pfxnames)
                            
                        elif filenamecomponents[0] == filenamecomponents[prefix_end]:
                            prefix = filenamecomponents[0].strip()
                        else:
                            prefix = '-'.join(filenamecomponents[:prefix_end])
                    prefix += '-'
                prefix += jobid
                if node.isdigit():
                    node = int(node)
                if prefix not in nodes_per_job:
                    nodes_per_job[prefix] = {}
                    jobfiles[prefix] = {}
                
                nodefiles = jobfiles[prefix]
                files_per_node = nodes_per_job[prefix]   

                if node not in files_per_node:
                    files_per_node[node] = [-1,-1,-1]
                    nodefiles[node] = [-1, -1, -1]
                    
                if "energy-tags" in file:
                    try:
                        files_per_node[node][1] = load_etag_file(file)
                        nodefiles[node][1] = file
                    except Exception as e:
                        print("File", file, "couldn't be loaded:", sys.exc_info())
                        sys.exit(1)
                elif "powercap-tags" in file:
                    if "CLEAN_" not in file:
                        nodefiles[node][2] = file
                        file = clean_pcaptag_file(file)
                    try:
                        files_per_node[node][2] = load_pcaptag_file(file)
                        
                    except Exception as e:
                        print("File", file, "couldn't be loaded")                    
                else:
                    if "CLEAN_" not in file:
                        nodefiles[node][0] = file
                        file = clean_poll_file(file) #for now
                    try:
                        files_per_node[node][0] = load_poll_file(file)
                    except Exception as e:
                        print("File", file, "couldn't be loaded")
                        print("ERROR", str(e))
                nodes_per_job[prefix] = files_per_node
                jobfiles[prefix] = nodefiles
    
    return nodes_per_job, jobfiles


def collect_energy_tags (path, jobfiles):
    for job, nodes in jobfiles.items():
        for node, files in nodes.items():
            etagfile = load_etag_file(files[1])
            for tag_name in etagfile['Tag name']:
                if tag_name != "application_summary":
                    loag_etag_from_name(tag_name, files[0])

def clean_poll_file (file):
    newfile = "CLEAN_"+file
    with open(os.path.join(args.prefix_path, file), 'r') as infile:
        with open(os.path.join(args.prefix_path, newfile), 'w') as ofile:
            lines = infile.readlines()
            for line in lines:
                if '---' in line or '***' in line:
                    continue
                else:
                    if len(line.split('\t')) > 38 and "energy-tags" not in file and "powercap-tags" not in file:
                        newline = []
                        for field in line.split('\t')[2:]:
                            if '.' in field:
                                newline.append(field)
                        newline = '\t'.join(line.split('\t')[:2]+newline)
                        if len(newline.split('\t')) > 38:
                        #    print("something is wrong with this file:", file)
                            continue
                        ofile.write(newline)
                    else:
                        ofile.write(line)
                    
    return newfile

def clean_pcaptag_file (file):
    newfile = "CLEAN_"+file
    with open(os.path.join(args.prefix_path, file), 'r') as infile:
        with open(os.path.join(args.prefix_path, newfile), 'w') as ofile:
            lines = infile.readlines()
            header = lines[0]
            ofile.write(header)
            headerfields = header.split('\t')
            headerlen = len(headerfields)
            headerfields = [hf.strip(' \n') for hf in headerfields]
            try:
                listbegin_index = headerfields.index("Emon tag list")
            except Exception as e:
                print("Warning: This file may not be formatted correctly")
            for line in lines[1:]:
                linefields = line.split('\t')
                if len(linefields) > headerlen:
                    etaglist = ','.join(linefields[listbegin_index:])
                    ofile.write('\t'.join(linefields[:listbegin_index]) + '\t'+etaglist)
                else:
                    ofile.write(line)
            
    return newfile

def load_poll_file (file):
    path = os.path.join(args.prefix_path, file)
    return pd.read_csv(path, sep='\t', header=0, parse_dates=[1], index_col = 0, na_values=[-1,'-1.000000'])

def load_etag_file (file):
    path = os.path.join(args.prefix_path, file)
    dateparse = lambda x: pd.datetime.strptime(x, '%Y-%m-%d %H:%M:%S')
    df = pd.read_csv(path, sep="\t", delimiter="\t", index_col=False, parse_dates=['Timestamp'], header=0, na_values=[-1,'-1.000000'])
    df = df.set_index('Tag Name')
    return df

def load_pcaptag_file (file):
    path = os.path.join(args.prefix_path, file)
    return pd.read_csv(path, sep='\t', header=0, dtype={'Emon tag list': np.str}, parse_dates=[2], index_col = 0, na_values=[-1,'-1.000000'])

def get_pcap_markers(df):
    tcol = 'Time since start (s)'
    flag = 'PCAP FLAG'
    pcol = 'Power Cap Long (W)'
    try:
        udf = df[df[flag] > 0]
    except Exception as e:
        return [], []
    return udf[tcol].values, udf[pcol].values
    

def get_etag_markers(df):
    try:
        udf = df.drop(['application_summary'])
    except Exception as e:
        print(e)
        return [],[],[]
    start = 'Start Time (s)'
    end = 'End Time (s)'    
    return udf[start].values, udf[end].values, udf.index.values
    
'''
    i = 0
    with open(os.path.join(path, pollfile), 'r') as infile:
        lines = infile.readlines()
        current_line = -1
        pair = [-1, -1]
        start_time = 0
        for line in lines:
            current_line+=1
            if line.split('\t')[0].strip() == "0":
                start_time = float(line.split('\t')[2].strip())
            if '---' in line:
                etag_name = line.split(':')[1].strip()
                if etag_name != "application_summary":
                    nextline = lines[current_line + 1]
                    j=2
                    while '---' in nextline or '***' in nextline:
                        nextline = lines[current_line + j]
                        j+=1
                    time = nextline.split('\t')[2].strip()
                    time = float(time)
                    if 'START' in line:
                        pair[0] = time
                    if 'END' in line:
                        pair[1] = time
                        if pair[0] == -1:
                            print("Warning: End of tag found without start. Will set to start time")
                            pair[0] = start_time
                        etag_markers.append(pair)
                        etag_names.append(etag_name)
                        pair = [-1, -1]
                        i+=1
            
    return etag_markers, etag_names
'''



def plot_aggregate_poll(path, nodes_per_job, jobfiles, columns, xcolumn, kind="line", xlabel="", ylabel="", etag_markers = False, pcap_markers = False):
    print("Plotting aggregate poll information")
    success = 0
    for job, files_per_node in nodes_per_job.items():
        fig, ax = plt.subplots()
        totalset = []
        for node, dataset in files_per_node.items():
            polldata = dataset[0]
            try:
                totalset.append(polldata)
                success = 1
            except Exception as e:
                print("Node", node, "does not contain a valid data set")
        if success == 1:
            totalset = pd.concat(totalset, axis=0)
            if xcolumn == "Time since start (s)":
                totalset[xcolumn] = totalset[xcolumn].round(rounding)
            totalset = totalset[totalset[xcolumn] < 310]
            for column in columns:
                totalset = totalset[totalset[column] < 10000]
            grouped = totalset.groupby(xcolumn)
            summed = grouped[columns].aggregate(sum)
            summed = summed.dropna()
            #print(summed[summed["RAPL dram P (W)"] > 10000])
            maxy = summed.max(axis=1).max()
            maxx = summed.index.max()
            print("Energy values per tag for job:", args.output_prefix + "_" +job)
            for column in columns:
                ax = summed.reset_index().plot(ax=ax, x=xcolumn, y=column, kind=kind, label=column)
                if '(W)' in column:
                    print(column, " - Total energy:", np.trapz(y=summed.as_matrix([column]).flatten(), x=summed.reset_index().as_matrix([xcolumn]).flatten()))  
            if etag_markers == True:
                texts = []
                etag_stimes, etag_etimes, etag_names = get_etag_markers(nodes_per_job[job][node][1])
                i = 0
                for start, end, name in zip(etag_stimes, etag_etimes, etag_names):
                    ax.axvline(start, color='gray')
                    ax.axvline(end, color='gray')
                    tagname = etag_names[i]
                    if '_' in tagname:
                        tagnameparts = tagname.split('_')
                        tagname = '\_'.join(tagnameparts)
                    texts.append(text(start, maxy, tagname, verticalalignment='center', fontdict={'family':'monospace', 'size':6}))
                    texts.append(text(end, maxy, tagname+" end", verticalalignment='center', fontdict={'family':'monospace', 'size':6}))
                    i+=1
                adjust_text(texts, arrowprops=dict(arrowstyle="->"))
            if pcap_markers == True:
                texts = []
                pcap_times, pcap_labels = get_pcap_markers(nodes_per_job[job][node][2])
                if len(pcap_times) > 0 and len(pcap_labels) > 0:
                    for time, label in zip(pcap_times, pcap_labels):
                        if time > maxx - 0.1*maxx:
                            continue
                        ax.axvline(time, color="gray")
                        text(time, maxy+maxy*0.1, label, verticalalignment="center", rotation=90)
                    
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            ax.set_ylim([0,9000])
            box = ax.get_position()
            #ax.set_position([box.x0, box.y0 - box.height * 0.1, box.width, box.height])
            legend = ax.legend(loc="upper center", bbox_to_anchor=(0.5,1.4), ncol=round(len(columns)/2))
            legend.get_frame().set_linewidth(0.5)
            ax.grid(False)
            ax = tx.format_axes(ax)
            fig.savefig(os.path.join(args.prefix_path, args.output_prefix+"aggregate_"+job+"_"+str(node)+"_" +str(len(files_per_node))+figformat), bbox_extra_artists=(legend,), bbox_inches='tight')
            plt.close()
        

def plot_boxplot(nodes_per_job, groupcol, column, xlabel="", ylabel="", user_prefix="", offset=10):
    for job, files_per_node in nodes_per_job.items():
        fig, ax = plt.subplots()
        totalset = []
        numnodes = len(files_per_node)
        for node, dataset in files_per_node.items():
            polldata = dataset[0]
            try:
                #polldata[groupcol] = polldata[groupcol].round(1)
                #polldata = polldata.iloc[::10,:]
                totalset.append(polldata)
                success = 1
            except Exception as e:
                print("Node", node, "does not contain a valid data set")        
                print(e)
                
        totalset = pd.concat(totalset)
        totalset = totalset[[groupcol, column]]
        totalset[groupcol] = totalset[groupcol].round(1)
        totalset[groupcol] = pd.cut(totalset[groupcol],offset)
        bp = totalset.boxplot(column=column, by=groupcol, ax=ax)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        ax.grid(False)
        ax = tx.format_axes(ax)
        for tick in ax.get_xticklabels():
            tick.set_rotation(45)
            tick.set_fontsize(6)
        bp.get_figure().gca().set_title("")
        fig.texts = []
        fig.suptitle('')
        fig.savefig(os.path.join(args.prefix_path, args.output_prefix+user_prefix+"boxplot_"+get_filetype(0)+"_"+job+"_"+str(node)+"_" +str(numnodes)+figformat), bbox_inches='tight')
        plt.close()
            


def plot_for_multinodes(nodes_per_job, fileindex, xcolumn, ycolumn, kind="line", xlabel="", ylabel="", row=-1, user_prefix=""):
    print("Plotting for multiple nodes")
    success = True
    for job, files_per_node in nodes_per_job.items():
        fig, ax = plt.subplots()
        numnodes = len(files_per_node)
        for node, dataset in files_per_node.items():
            polldata = dataset[fileindex]
            if row != -1 and fileindex == 1:
                try:
                    xval = [polldata.loc[row,xcolumn]]
                    yval = [polldata.loc[row,ycolumn]]
                    ax.plot(xval, yval, 'o', label=node)
                except Exception as e:
                    print("There is no data frame for this file")
                    success = False
            else:
                try:
                    polldata = polldata[polldata[xcolumn] < 300]
                    polldata = polldata[polldata[ycolumn] < 10000]
                    if numnodes > 10:
                        ax = polldata.plot(ax=ax, kind=kind, x=xcolumn, y=ycolumn, legend=False)
                        #polldata = polldata.iloc[::10,:]
                        #ax = polldata.boxplot(ycolumn, xcolumn, ax=ax)
                    else:
                        ax = polldata.plot(ax=ax, kind=kind, x=xcolumn, y=ycolumn, label=node)
                except Exception as e:
                    print("There is no data frame for this file")
                    success = False
        if success == True:
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
            ax.grid(False)
            ax = tx.format_axes(ax)            
            if numnodes <= 10:
                legend = ax.legend(loc="center left", bbox_to_anchor=(1,0.5))
                fig.savefig(os.path.join(args.prefix_path, args.output_prefix+user_prefix+"plot_"+get_filetype(fileindex)+"_"+job+"_"+str(node)+"_" +str(numnodes)+figformat), bbox_extra_artists=(legend,), bbox_inches='tight')
            else:
                fig.savefig(os.path.join(args.prefix_path, args.output_prefix+user_prefix+"plot_"+get_filetype(fileindex)+"_"+job+"_"+str(node)+"_" +str(numnodes)+figformat), bbox_inches='tight')
            plt.close()

    return

def get_filetype(fileindex):
    if fileindex == 0:
        return "poll"
    if fileindex == 1:
        return "etag"
    if fileindex == 2:
        return "ptag"

def plot_energy_distribution (nodes_per_job, column, columnname = 'column', xlabel = "", ylabel = "", row=0, bins=10):
    print("Plotting energy distribution")
    success = True
    vals = []
    for job, files_per_node in nodes_per_job.items():
        fig, ax = plt.subplots()
        for node, dataset in files_per_node.items():
            polldata = dataset[1]
            try:
                vals.append(polldata.iloc[row,column])
            except Exception as e:
                print("There is no data frame for this file")
                success = False
        if success == True:    
            df = pd.DataFrame({columnname:vals})
            df.hist(ax=ax, column=0, bins=bins)
            title = "\#nodes: " + str(len(files_per_node))
            ax.set_title(title)
            ax.set_xlabel(columnname)
            ax.grid(False)
            ax = tx.format_axes(ax)
            fig.savefig(os.path.join(args.prefix_path, args.output_prefix+"distribution_"+job+"_"+str(len(files_per_node))+figformat), bbox_inches='tight')
            plt.close()

def split_label(x):
    if '_' in x:
        x = "\\shortstack{"+x+"}"
    return x


def get_total_energy_from_tags(nodes_per_job, jobfiles, ep_columns, perf_columns=[], plot=True, xlabel = "", ylabel1 = "", ylabel2 = "", exclude_rows=[], user_prefix=""):
    print("Plotting total energy values for tags")
    for job, files_per_node in nodes_per_job.items():
        total_set = []
        fig, ax = plt.subplots()
        for node, dataset in files_per_node.items():
            df = dataset[1][ep_columns + perf_columns]
            #for LAMMPS with 100 steps only
            #looprows = ['force_atom', 'force_kspace', 'force_pair', 'neighbor_setup']
            #df.loc[looprows] = df.loc[looprows] / 100.0
            #df.loc['neighbor_build'] = df.loc['neighbor_build'] / 11.0
            #df.loc['write_output'] = df.loc['write_output'] / 2.0            
            total_set.append(df)
        total_set = pd.concat(total_set, axis=0)
        grouped = total_set.groupby(total_set.index)
        if "Total Time (s)" in ep_columns:
            summed = grouped[ep_columns].aggregate(max)
        else:
            summed = grouped[ep_columns].aggregate(sum)
        if len(exclude_rows) > 0:
            summed.drop(exclude_rows, inplace=True)
         
        summed.index = summed.index.str.replace('_', '\\\\\_')
        summed.index = summed.index.map(split_label)
        summed.plot(ax=ax, kind='bar', position=1)
        
        summed_perf = pd.DataFrame()
        if len(perf_columns) > 0:
            summed_perf = grouped[perf_columns].aggregate(max)
            print(summed_perf)
            summed_perf.index = summed_perf.index.str.replace('_','\_')
            ax2 = ax.twinx()
            ax2.set_ylabel(ylabel2)
            summed_perf.plot(ax=ax2, kind='bar', position=0)
            box2 = ax2.get_position()
            ax2.set_position([box2.x0, box2.y0, box2.width * 0.8, box2.height])
            legend2 = ax2.legend(loc="lower left", bbox_to_anchor=(1,0.5))
        numnodes = len(nodes_per_job[job])
        nodename = ""
        if numnodes == 1:
            nodename = str(node)+"_"
        if plot == True:
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel1)
            box = ax.get_position()
            #ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
            legend = None #ax.legend(loc="best", bbox_to_anchor=(0.6,0.5))
            tx.format_axes(ax)
            fn = "plot-tag-summary_"
            if len(exclude_rows) > 0:
                fn = "user-tags-"+fn
            fig.savefig(os.path.join(args.prefix_path, args.output_prefix+user_prefix+fn+job+"_"+nodename+str(numnodes)+figformat), bbox_extra_artist=(legend,), bbox_inches='tight')
            plt.close()


def get_total_energy(nodes_per_job, jobfiles, column, tag):
    print("Printing Total Energy for tag", tag, "and components", column)
    for job, files_per_node in nodes_per_job.items():
        total_energy = 0
        for node, dataset in files_per_node.items():
            etagdata = dataset[1] 
            se = etagdata[column]
            total_energy += se[tag]
            
        print(job, args.output_prefix, column, "Tag:", tag, "Total energy: ", total_energy)
            
            
    return

def count_tags(nodes_per_job):
    for job, files_per_node in nodes_per_job.items():
        for node, files in files_per_node.items():
            print("JOB:", job)
            print(files[1].groupby(files[1].index).size())
            print("=====================")

def plot_mpirank_vs_energy(nodes_per_job, tag, columns):
    mpiranks = []
    totalset = []
    for job, files_per_node in nodes_per_job.items():
        mpiranks.append(int(job.split('-')[0]))
        for node, files in files_per_node.items():
            totalset.append(files[1])
    totalset = pd.concat(totalset)
    #mpiranks = pd.Series(mpiranks, name="MPI ranks")
    
    grouped = totalset.groupby(totalset.index)[columns].aggregate(lambda x: list(x))
    fig, ax = plt.subplots()
    for column in columns:
        y = grouped.loc[tag, column]
        ax.plot(mpiranks, y, 'o', label=column)
    ax.legend()
    plt.show()
    

nodes_per_job, jobfiles = collect_files(args.prefix_path)

#plot_mpirank_vs_energy(nodes_per_job, "application_summary", ["Total RAPL pkg E (J)", "Total RAPL dram E (J)", "Total Cray node E (J)", "Total Cray cpu E (J)", "Total Cray memory E (J)"])

#get_total_energy_from_tags(nodes_per_job, jobfiles, ["Total RAPL dram E (J)", "Total Cray memory E (J)"], xlabel = "", ylabel1 = "Total Energy (J)", ylabel2 = "Total Time (s)")#, exclude_rows = ['application_summary', 'verlet_run'])

#get_total_energy_from_tags(nodes_per_job, jobfiles, ["Total RAPL pkg P (W)", "Total RAPL dram P (W)", "Total Cray node calc P (W)", "Total Cray cpu calc P (W)", "Total Cray memory calc P (W)"], xlabel = "", ylabel1 = "Total Energy (J)", ylabel2 = "Total Time (s)")#, exclude_rows = ['application_summary', 'verlet_run'])

looprows = ['force_atom', 'force_kspace', 'force_pair', 'neighbor_setup', 'verlet_cleanup', 'write_output', 'application_summary']

get_total_energy_from_tags(nodes_per_job, jobfiles, ["Total RAPL pkg E (J)", "Total RAPL dram E (J)", "Total Cray node E (J)", "Total Cray cpu E (J)", "Total Cray memory E (J)"], xlabel = "", ylabel1 = "Total Energy (J)", user_prefix="big_energy_",exclude_rows = looprows)#['application_summary', 'verlet_run', 'verlet_setup'])
#get_total_energy_from_tags(nodes_per_job, jobfiles, ["Total RAPL pkg P (W)", "Total RAPL dram P (W)", "Total Cray node calc P (W)", "Total Cray cpu calc P (W)", "Total Cray memory calc P (W)"], xlabel = "", ylabel1 = "Total Power (W)", user_prefix="_power_",exclude_rows = ['application_summary'])#, 'verlet_run'])
#get_total_energy_from_tags(nodes_per_job, jobfiles, ["Total Time (s)"], xlabel = "", ylabel1 = "Total Time (s)", user_prefix="_time_",exclude_rows = ['application_summary', 'verlet_run', 'verlet_setup'])


#get_total_energy(nodes_per_job, jobfiles, "Total RAPL pkg E (J)", "read")
#get_total_energy(nodes_per_job, jobfiles, "Total RAPL pkg E (J)", "write")
#get_total_energy(nodes_per_job, jobfiles, "Total RAPL pkg E (J)", "application_summary")

#plot_for_multinodes(nodes_per_job, 0, "Time since start (s)", "RAPL pkg P (W)", xlabel="Time since start (s)", ylabel="RAPL pkg P (W)", user_prefix="-RAPLpkgP-")
#plot_for_multinodes(nodes_per_job, 0, "Time since start (s)", "RAPL pkg E (J)", xlabel="Time since start (s)", ylabel="RAPL pkg E (J)", user_prefix="-RAPLpkgE-")
#plot_for_multinodes(nodes_per_job, 0, "Time since start (s)", "Cray cpu E (J)", xlabel="Time since start (s)", ylabel="Cray cpu E (J)", user_prefix="-CrayCPUE-")
#plot_for_multinodes(nodes_per_job, 0, "Time since start (s)", "Cray cpu P (W", xlabel="Time since start (s)", ylabel="Cray cpu P (W)", user_prefix="-CrayCPUP-")
#plot_for_multinodes(nodes_per_job, 0, "Time since start (s)", "Cray node P (W)", xlabel="Time since start (s)", ylabel="Cray node P (W)", user_prefix="-CrayNodeP-")
#plot_for_multinodes(nodes_per_job, 0, "Time since start (s)", "Cray node E (J)", xlabel="Time since start (s)", ylabel="Cray node E (J)", user_prefix="-CrayNodeE-")


#plot_for_multinodes(nodes_per_job, 1, 3, 8, kind="line", xlabel="Total Time (s)", ylabel=" Total RAPL pkg P (W)", row=0)
#plot_for_multinodes(nodes_per_job, 1, 'Total Time (s)', 'Total Cray node E (J)', xlabel="Total Time (s)", ylabel="Total Node Energy (J)", row='application_summary')  
#plot_energy_distribution(nodes_per_job, 8, columnname = "Total Cray node E (J)", bins=15)

#plot_aggregate_poll(args.prefix_path, nodes_per_job, jobfiles, ["RAPL pkg P (W)", "RAPL dram P (W)", "Cray node P (W)", "Cray cpu P (W)", "Cray memory P (W)"], "Time since start (s)", xlabel="Time since start (s)", ylabel="Total Power (P)", etag_markers=True, pcap_markers=False)

#plot_aggregate_poll(args.prefix_path, nodes_per_job, jobfiles, ["RAPL pkg E (J)", "RAPL dram E (J)", "Cray node E (J)", "Cray cpu E (J)", "Cray memory E (J)"], "Time since start (s)", xlabel="Time since start (s)", ylabel="Total Energy (J)", etag_markers=True)
#plot_aggregate_poll(args.prefix_path, nodes_per_job, jobfiles, ["RAPL dram E (J)"], "Time since start (s)", xlabel="Time since start (s)", ylabel="Total Energy (J)", etag_markers=True)

#plot_boxplot(nodes_per_job, "Time since start (s)", "RAPL pkg P (W)", xlabel="Time since start (s)", ylabel="RAPL Package Power (W)")


#count_tags(nodes_per_job)


'''
STREAM TAGS DDR
texts.append(text(start - 2.5, maxy+0.2*maxy, tagname, verticalalignment='center', fontdict={'family':'monospace', 'size':6}))
texts.append(text(end - 6, maxy+0.15*maxy, tagname+" end", verticalalignment='center', fontdict={'family':'monospace', 'size':6}))

STREAM TAGS MCDRAM
texts.append(text(start-0.5, maxy+0.1*maxy, tagname, verticalalignment='center', fontdict={'family':'monospace', 'size':6}))
if tagname == 'copy':
    texts.append(text(end - 1.5, maxy+0.05*maxy, tagname+" end", verticalalignment='center', fontdict={'family':'monospace', 'size':6}))
else:
    texts.append(text(end - 1, maxy+0.05*maxy, tagname+" end", verticalalignment='center', fontdict={'family':'monospace', 'size':6}))
'''

