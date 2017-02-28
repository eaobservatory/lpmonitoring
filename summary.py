# Copyright (C) 2015-2016 East Asian Observatory
# All Rights Reserved.
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation; either version 2 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful,but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# this program; if not, write to the Free Software Foundation, Inc.,51 Franklin
# Street, Fifth Floor, Boston, MA  02110-1301, USA

from __future__ import absolute_import, division, print_function
from collections import namedtuple, OrderedDict
import datetime
import os

import matplotlib
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib import gridspec
import numpy as np
import StringIO
from io import BytesIO

from flask import send_file
from omp.db.db import OMPDB
from omp.siteconfig import get_omp_siteconfig
from omp.obs.state import OMPState

# Object to hold info about a project
ProgInfo = namedtuple('ProgInfo', 'projectid name info instruments url color allocdict')

class Bands:
    """ class for weather info"""
    BAND1 = [0,0.05]
    BAND2 = [0.05, 0.08]
    BAND3 = [0.08, 0.12]
    BAND4 = [0.12, 0.2]
    BAND5 = [0.2, 100]


faultstatus = {
    0: 'Open',
    1: 'Closed',
    2: 'WorksForMe',
    3: 'NotAFault',
    4: 'WontBeFixed',
    5: 'Duplicate',
    6: 'OpenWillBeFixed',
    7: 'Suspended',
}

def get_bands(tau):
    """
    Turn lists of 225 GHz taus into JCMT bands.
    """

    tau = np.asarray(tau)
    bands = np.zeros(tau.shape)

    band1 = tau < 0.05
    band2 = (tau >= 0.05) & (tau < 0.08)
    band3 = (tau >= 0.08) & (tau < 0.12)
    band4 = (tau >= 0.12) & (tau < 0.2)
    band5 = (tau >= 0.2)

    bands[band1] = 1
    bands[band2] = 2
    bands[band3] = 3
    bands[band4] = 4
    bands[band5] = 5

    return bands

class LP:
    """
    class for handling information about an LP
    """
    TRANSIENT = ProgInfo('M16AL001',
                        'Transient',
                        'Protostar variability: observations of 8 sources, repeated ~ every 21 days',
                        'SCUBA-2',
                        'https://www.eaobservatory.org/jcmt/science/large-programs/transient/',
                         '#4D4D4D',
                         {'SCUBA-2':{1:50,2:50,3:50}})
    S2COSMOS = ProgInfo('M16AL002',
                        'S2COSMOS',
                        'Cosmology: 2 square degree uniform deep 850um map',
                        'SCUBA-2',
                        'https://www.eaobservatory.org/jcmt/science/large-programs/s2-cosmos/',
                        '#5DA5DA',
                        {'SCUBA-2':{2:111,3:112}})
    SCOPE = ProgInfo('M16AL003',
                     'SCOPE',
                     'Early protostars: survey of 1000 densest Planck clumps',
                     'SCUBA-2',
                     'https://www.eaobservatory.org/jcmt/science/large-programs/scope/',
                     '#FAA43A',
                     {'SCUBA-2':{3:150,4:150}})
    GBBFIELD = ProgInfo('M16AL004',
                        'BISTRO',
                        'POL-2 observations of dense GBS regions',
                        'SCUBA-2',
                        'https://www.eaobservatory.org/jcmt/science/large-programs/gb_bfields/',
                        '#60BD68',
                        {'SCUBA-2':{2:224}})
    JINGLE = ProgInfo('M16AL005',
                      'JINGLE',
                      'Continuum maps of 190 IR galaxies and CO(2-1) fluxes for subset of 75 galaxies',
                      'SCUBA-2, RXA3M',
                      'https://www.eaobservatory.org/jcmt/science/large-programs/jingle/',
                      '#F17CB0',
                      {'SCUBA-2':{2:57,3:123,4:75},'RxA3m':{4:125, 5:400}})
    STUDIES = ProgInfo('M16AL006',
                       'STUDIES',
                       'Ultra deep confusion limited daisy at 450um',
                       'SCUBA-2',
                       'https://www.eaobservatory.org/jcmt/science/large-programs/studies/',
                       '#B2912F',
                       {'SCUBA-2':{1:330}})
    MALATANG = ProgInfo('M16AL007',
                        'MALATANG',
                        'Harp HCN/HCO+ jiggles and lines from cores of a sample of nearby galaxies',
                        'HARP',
                        'https://www.eaobservatory.org/jcmt/science/large-programs/malatang/',
                        '#F15854',
                        {'HARP':{2:40,3:100,4:250}})

    PROJECTS = [TRANSIENT, S2COSMOS, SCOPE, GBBFIELD, JINGLE, STUDIES, MALATANG]

    @classmethod
    def lookup_code(cls, code):
        for i in cls.PROJECTS:
            if i.projectid == code:
                return i
        raise Exception('Unkown LP project code {0}'.format(code))


def create_summary(ompdb):
    """
    Get the summary information for the JCMT large programmes.

    """

    # Basic info for web page.
    webinfo = {}
    webinfo['title'] = 'Summary'

    # Ensure the LP.PROJECTS list is available.
    webinfo['LP'] = LP
    ompallocs = ompdb.get_allocation_project('M16AL%', like=True)
    timecharged = get_time_charged(ompdb)

    projsum = namedtuple('projsum',
        'summary msbs faults number_observations time_observations charged unconfirmed allocation percentcomplete')
    faultsum = namedtuple('faultsum', 'project numopen total')
    chargeinfo = namedtuple('summary', 'percent timecharged timeallocated')

    summinfo = ompdb.get_summary_obs_info('M16AL%')
    msbinfo = ompdb.get_summary_msb_info('M16AL%')
    faultinfo = ompdb.get_fault_summary('M16AL%')

    projdict = {}

    for p in LP.PROJECTS:
        # Time charged information
        charged = sum([i.timespent for i in timecharged[p.projectid]]) / (60.0*60.0)
        unconfirmed = charged - sum([i.timespent  for i in timecharged[p.projectid] if i.confirmed]) / (60.0*60.0)
        allocation = ompallocs[p.projectid].allocated/ (60.0*60.0)
        percentcomplete = 100.0 * charged/allocation

        # Observation information (all obs and bands)
        summary =  [i for i in summinfo if i.project == p.projectid]
        number_observations = {}
        number_observations['good'] = sum([ i.number for i in summary
                                            if i.status==OMPState.GOOD])
        number_observations['bad'] = sum([ i.number for i in summary
                                           if i.status==OMPState.BAD])
        number_observations['questionable'] = sum([ i.number for i in summary
                                                    if i.status==OMPState.QUESTIONABLE])
        number_observations['rejected'] = sum([ i.number for i in summary
                                                if i.status==OMPState.REJECTED])
        number_observations['junk'] = sum([ i.number for i in summary
                                            if i.status==OMPState.JUNK])
        time_observations = {}
        time_observations['good'] = sum([ i.totaltime for i in summary
                                            if i.status==OMPState.GOOD])
        time_observations['bad'] = sum([ i.totaltime for i in summary
                                           if i.status==OMPState.BAD])
        time_observations['questionable'] = sum([ i.totaltime for i in summary
                                                    if i.status==OMPState.QUESTIONABLE])
        time_observations['rejected'] = sum([ i.totaltime for i in summary
                                                if i.status==OMPState.REJECTED])
        time_observations['junk'] = sum([ i.totaltime for i in summary
                                            if i.status==OMPState.JUNK])
        time_observations['uncharged'] = sum([ i.totaltime for i in summary
                                              if i.status in (OMPState.JUNK, OMPState.REJECTED, OMPState.BAD)
                                           ])

        # Instrument and band.
        instruments = [i.instrument for i in summary]
        msbs = [i for i in msbinfo if i.project == p.projectid]
        faults = [i for i in faultinfo if i.project.lower() == p.projectid.lower()]
        f = faultsum(p.projectid, len([i for i in faults if i.status in [0,6]]),
                     len(faults))
        projdict[p.projectid] = projsum(summary, msbs, f, number_observations, time_observations,
                                        charged, unconfirmed, allocation, percentcomplete)

    webinfo['projdict'] = projdict
    return webinfo


def get_project_information(ompdb, projectid, instrument=None):

    """
    Get detailed information about each obs in a project (slow).
    """

    obs = ompdb.get_observations_from_project(projectid, instrument=instrument)
    instdict = {}
    wvmtaus = [(i.wvmtaust +  i.wvmtauen)/2.0 for i in obs.values()]
    bands = get_bands(wvmtaus)
    times = np.asarray([i.duration for i in obs.values()])
    instruments = np.asarray([i.instrument for i in obs.values()])
    status = np.asarray([i.status for i in obs.values()])
    source = np.asarray([i.object for i in obs.values()])
    for inst in set(instruments):
        instdict[inst] = {}
        for b in set(bands):
            instdict[inst][b] = {}
            instdict[inst][b]['good'] = times[(instruments==inst)&(bands==b)&(status==OMPState.GOOD)]
            instdict[inst][b]['bad'] = times[(instruments==inst)&(bands==b)&(status==OMPState.BAD)]
            instdict[inst][b]['questionable'] = times[(instruments==inst)&(bands==b)&(status==OMPState.QUESTIONABLE)]
            instdict[inst][b]['rejected'] = times[(instruments==inst)&(bands==b)&(status==OMPState.REJECTED)]
        instdict[inst]['questionable'] = [i for i in obs.values() if i.status == OMPState.QUESTIONABLE]
        instdict[inst]['bad'] = [i for i in obs.values() if i.status == OMPState.BAD]
        instdict[inst]['rejected'] = [i for i in obs.values() if i.status == OMPState.REJECTED]
        instdict[inst]['junk'] = [i for i in obs.values() if i.status == OMPState.JUNK]
        instdict[inst]['obs'] = [i for i in obs.values() if i.instrument == inst]
    return instdict


def get_omp_database():
    """
    Construct a read-only OMP database access object.

    Access comes from the omp site_config file.
    """

    config = get_omp_siteconfig()
    credentials = 'hdr_database'

    ompdb = OMPDB(
            server=config.get(credentials, 'server'),
            user=config.get(credentials, 'user'),
            password=config.get(credentials, 'password'),
            read_only=True)

    return ompdb

def get_time_charged(ompdb):
    timecharged = OrderedDict()
    for proj in LP.PROJECTS:
        p = proj.projectid
        timecharged[p] = ompdb.get_time_charged_project_info(p)
    return timecharged


def get_observations(ompdb, project):
    observations = OrderedDict()
    for proj in LP.PROJECTS:
        p = proj.projectid
        observations[p] = ompdb.get_observations_from_project(p)
    return observations

def prepare_completion_chart(ompdb):
    """
    Create a plot showing the percentage completion of the large projects.

    Return sendfile object of mime-type image/png
    """

    # Create matplotlib figure
    figure = Figure(figsize=(15,5))
    fig = make_completion_chart(ompdb, figure)

    # Create sendfile object
    canvas = FigureCanvas(fig)
    img = StringIO.StringIO()
    canvas.print_png(img)
    img.seek(0)
    return send_file(img, mimetype='image/png')




def make_completion_chart(ompdb, fig):
    """
    Create plot showing completion of projects.
    Return a matplotlib figure
    """
    ompallocs = ompdb.get_allocation_project('M16AL%', like=True)
    timecharged = get_time_charged(ompdb)

    # Date range to plot over.
    sdate = datetime.datetime(2015,11,25,0,0)
    edate = datetime.datetime(2017,1,1,0,0)

    # Current date for title
    currdate = datetime.datetime.utcnow().strftime('%Y-%m-%d %H:%M')

    # Set up two plots on one figure.
    gs = gridspec.GridSpec(1,2, width_ratios=[5,2])
    ax = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])


    # Plot each projs time charged per day as cumulative percentage observed
    for p in LP.PROJECTS:
        dates = [i.date for i in timecharged[p.projectid]]
        times = [i.timespent for i in timecharged[p.projectid]]
        cumultimes = [sum(times[0:i+1])for i,j in enumerate(times)]
        allocation = ompallocs[p.projectid].allocated
        color = p.color
        label = p.name

        ax.plot_date(dates, [100* i/allocation for i in cumultimes],linestyle='solid', marker='x',color=color, label=label,
                linewidth=2.0, markeredgewidth=2.0)

    # set up plot limits and labels
    ax.set_ylim(0, 100)
    #ax.set_xlim(sdate,edate)
    ax.set_ylabel('Percentage completion')
    fig.patch.set_visible(False)
    fig.autofmt_xdate()
    ax.tick_params(axis='x', which='minor', bottom='on', top='on')
    ax.xaxis.set_minor_locator(matplotlib.dates.MonthLocator())
    ax.set_title('Completion as of ' + currdate)
    ax.yaxis.grid()

    # Get times of LAP blocks
    lapfile = 'LAP-UT-blocks.txt'
    if os.path.isfile(lapfile):
        print('found LAP file')
        f = open('LAP-UT-blocks.txt', 'r')
        blocks = f.readlines()
        f.close()
        blocks = [i.strip() for i in blocks]
        blocks = [i.split() for i in blocks]
        for start, end, project in blocks:
            start = datetime.datetime(int(start[0:4]), int(start[4:6]), int(start[6:8]), 0, 0)
            end =  datetime.datetime(int(end[0:4]), int(end[4:6]), int(end[6:8]), 23, 59)
            p = LP.lookup_code(project.upper())
            ax.axvspan(start, end, alpha=0.25, color=p.color, edgecolor=None)


    # Set up second plot: bar chart of percentage observed.
    observedtimes = [sum([i.timespent for i in t]) for t in timecharged.values()]
    allocations=[i.allocated for i in ompallocs.values()]
    total = [100 for i in observedtimes]
    percent_observed = [100*j/allocations[i]for i,j in enumerate(observedtimes)]

    # Plot a bar for each project.
    ax2.barh([i+0.1 for i in range(0,7)],total, label='allocated',
             alpha=0.3, color=[i.color for i in LP.PROJECTS],
             edgecolor='white', height=0.8)
    ax2.barh([i+0.1 for i in range(0,7)],percent_observed, label='observed',
             alpha=1.0, color=[i.color for i in LP.PROJECTS], edgecolor='white',
             height=0.8)

    # Label each bar
    for i in range(len(observedtimes)):
        ax2.text(90, i+0.5, LP.PROJECTS[i].name + ' %.1F/%.1F hrs'%(observedtimes[i]/(60.0*60.0), allocations[i]/(60.0*60.0)), color='black', ha='right', alpha=1.0, va='center', weight='bold')

    # Finish up figure.
    ax2.set_axis_off()
    ax2.invert_yaxis()
    fig.set_tight_layout(tight=True)

    return fig




def prepare_project_chart(project, summary):
    """
    Create chart for project obs.
    """
    try:
        pinfo = LP.lookup_code(project)
    except:
        pinfo = None

    # Number of charts
    instruments = list(set([i.instrument for i in summary]))

    figsize = (5*len(instruments), 5)

    # Set up chart on one figure.
    fig = Figure(figsize=figsize)

    # Go through each instrument.
    for instno in range(len(instruments)):
        instrument = instruments[instno]
        ax = fig.add_subplot(1, len(instruments), instno+1)
        data = [i for i in summary if i.instrument == instrument]
        gooddata = 5*[0]
        questdata = 5*[0]
        for b in range(0,5):
            band = b+1
            gooddata[b] = sum([i.totaltime/(60*60.0) if (int(i.band)==band)&(i.status==OMPState.GOOD) else 0.0 for i in data ])
            questdata[b] = sum([i.totaltime/(60*60.0) if (int(i.band)==band)&(i.status==OMPState.QUESTIONABLE) else 0.0 for i in data ])

        if pinfo:
            color = pinfo.color
        else:
            color = None
        ax.bar([i+0.05 for i in range(0,5)],
               gooddata, color=color,
                edgecolor='white', width=0.9)
        ax.bar([i+0.05 for i in range(0,5)],
                questdata,
                bottom=gooddata,
                color = 'orange',
                edgecolor='white', width=0.9)

        ax.set_xticks([i+0.5 for i in range(0,5)])
        xlabels = []
        for i in range(0,5):
            label = 'B%i'%(i+1)
            xlabels.append(label)

        ax.set_xticklabels(xlabels, weight='bold')

        if pinfo:
            for i in range(0,5):
                allocation = pinfo.allocdict.get(instrument, {}).get(i+1,0)
                label = '%.1f' %(gooddata[i])
                if questdata[i] != 0:
                    label += ' G \n %.1F Q \n' %(questdata[i])
                else:
                    label += '\n'
                label += ' /%i hr' % allocation
                if allocation > 0:
                    label += '\n'
                    label += '%.1f' %(100.0*(gooddata[i] + questdata[i])/allocation)
                    label += '%'
                ax.text(i+0.5, 3, label,
                        color='black', ha='center', alpha=1.0, va='bottom', weight='bold', rotation='horizontal')

        #ax.set_ylim(0, max(ax.get_xlim()[1], max(pinfo.allocdict.get(instrument,{0:0}).values())))
        ax.set_title('{}: {:.1F} hrs observations on sky'.format(instrument, sum(gooddata) + sum(questdata))
                     , weight='bold')
        ax.set_ylabel('Hrs on sky')

    fig.patch.set_visible(False)
    fig.set_tight_layout(tight=True)

    # Create sendfile object
    canvas = FigureCanvas(fig)
    img = StringIO.StringIO()
    canvas.print_png(img)
    img.seek(0)
    return send_file(img, mimetype='image/png')

