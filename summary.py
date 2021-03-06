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
from flask import url_for
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
from project import create_msb_image

# Object to hold info about a project
ProgInfo = namedtuple('ProgInfo', 'projectid name info instruments url color allocdict')
MsbCustomInfo = namedtuple('MSBCustomInfo', 'project instrument target title coordstype ra2000 dec2000 remaining timeest taumin taumax wavelength')




class Bands:
    """ class for weather info"""
    BAND1 = [0,0.05]
    BAND2 = [0.05, 0.08]
    BAND3 = [0.08, 0.12]
    BAND4 = [0.12, 0.2]
    BAND5 = [0.2, 100]

    @classmethod
    def get_band(self, tau):
        """ Return the band a given tau is in."""
        if tau < self.BAND1[1]:
            return 1
        elif tau < self.BAND2[1]:
            return 2
        elif tau < self.BAND3[1]:
            return 3
        elif tau < self.BAND4[1]:
            return 4
        elif tau >= self.BAND5[0]:
            return 5
        else:
            return np.nan
    @classmethod
    def get_fromname(self, name):
        if int(name) == 1:
            return self.BAND1
        elif int(name) == 2:
            return self.BAND2
        elif int(name) == 3:
            return self.BAND3
        elif int(name) == 4:
            return self.BAND4
        elif int(name) == 5:
            return self.BAND5
        else:
            return 'unknown'


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





allocations = {}
allocations['M16AL001'] = {'SCUBA-2':{1:50,2:50,3:50}}
allocations['M16AL002'] = {'SCUBA-2':{2:111,3:112}}
allocations['M16AL003'] = {'SCUBA-2':{3:150,4:150}}
allocations['M16AL004'] = {'SCUBA-2':{2:224}}
allocations['M16AL005'] = {'SCUBA-2':{2:57,3:123,4:75},'RxA3m':{4:125, 5:400}}
allocations['M16AL006'] = {'SCUBA-2':{1:330}}
allocations['M16AL007'] = {'HARP':{2:40,3:100,4:250}}

allocations['M17BL005'] = {'SCUBA-2':{2:221}, 'HARP':{3:55.3}}
allocations['M17BL011'] = {'POL-2':{1:56, 2:168}}
allocations['M17BL009'] = {'SCUBA-2':{1:319}}
allocations['M17BL004'] = {'HARP':{1:85.5,2:218.4,4:50,5:50}}

allocations['M17BL002'] = {'SCUBA-2':{1:10, 2:125}, 'HARP':{3:18,4:106},
                           'RxA3m':{4:99,5:157}}
allocations['M17BL001'] = {'SCUBA-2':{2:436.5,3:436.5}}
allocations['M17BL010'] = {'RxA3m':{4:285,5:169}}
allocations['M17BL006'] = {'SCUBA-2':{2:152.6, 3:152.6}}
#allocations['M17BL007'] = {'SCUBA-2':{3:400}}


# In order of priority.
LP_1 = ['M16AL001', 'M16AL002', 'M16AL003', 'M16AL004', 'M16AL005','M16AL006',
               'M16AL007']
LP_2A = ['M17BL005','M17BL011', 'M17BL009', 'M17BL004', 'M17BL002']
LP_2B = ['M17BL001','M17BL010', 'M17BL006', 'M17BL007']


def create_summary(ompdb, semester='LAP', queue=None, patternmatch=None, projects=None, exclude_projects=None,
                   exclude_done=True,
                   details=True, blocks=True, allprojectmsb=True, queryparams=None):
    """
    Get the summary information for the JCMT large programmes.

    """

    originalprojects = projects[:]

    # Basic info for web page.
    webinfo = {}
    title = 'Summary of projects: '
    elements = []
    if semester:
        elements += ['semester={}'.format(semester)]
    if queue:
        elements += ['queue={}'.format(queue)]
    if patternmatch:
        elements += ['matching {}'.format(patternmatch)]
    if projects:
        elements += ['projects={}'.format(', '.join(projects))]
    title += ' and '.join(elements)
    webinfo['fulltitle'] = title
    webinfo['title'] = 'Summary'
    webinfo['allocations'] = allocations

    #Put original information into webinfo, for query_form.
    webinfo['semester'] = semester
    webinfo['queue'] = queue
    webinfo['patternmatch'] = patternmatch
    webinfo['requested_projects'] = originalprojects
    webinfo['exclude_projects'] = exclude_projects
    webinfo['details'] = bool(details)
    webinfo['exclude_done'] = bool(exclude_done)
    webinfo['blocks'] = bool(blocks)
    webinfo['allprojectmsb'] = bool(allprojectmsb)

    ompallocs = ompdb.get_allocations(semester=semester, queue=queue, projects=projects,
                                      patternmatch=patternmatch, telescope='JCMT')
    projects = list(ompallocs.keys())

    if exclude_projects:
        print('excluding projects!', exclude_projects)
        for p in exclude_projects:
            try:
                projects.remove(p)
            except ValueError:
                pass

    timecharged = ompdb.get_time_charged_group(projects=projects, telescope='JCMT')

    # Get the overall completion of this set of projects
    overallallocation = 0
    overallobserved = 0
    for p in projects:
        allocation = ompallocs[p].allocated
        observed = (allocation - ompallocs[p].remaining) + ompallocs[p].pending
        overallallocation += allocation
        overallobserved += observed
    webinfo['overallallocation_hrs'] = overallallocation/(60.0*60.0)
    webinfo['overallobserved_hrs'] = overallobserved/(60.0*60.0)

    projsum = namedtuple('projsum',
        'summary msbs faults number_observations time_observations charged unconfirmed allocation percentcomplete')

    faultsum = namedtuple('faultsum', 'project numopen total')
    chargeinfo = namedtuple('summary', 'percent timecharged timeallocated')

    projdict = {}
    projinfodict = {}

    webinfo['ompallocations'] = ompallocs
    webinfo['timecharged'] = timecharged

    # msbinfo needs to be gotten either for details==True or for
    # allprojectmsbs==True. Don't want to do it twice though.
    msbinfo = None
    if details:
        summinfo = ompdb.get_summary_obs_info_group(projects=projects)
        end3 = datetime.datetime.now()

        msbinfo = ompdb.get_summary_msb_info_group(projects=projects)
        faultinfo = ompdb.get_fault_summary_group(projects=projects)
        for p in projects:
            projinfodict[p] = ompdb.get_project_info(p)
        webinfo['projinfodict'] = projinfodict



        for p in projects:

            # Time charged information
            if p in timecharged:
                charged = sum([i.timespent for i in timecharged[p]]) / (60.0*60.0)
                unconfirmed = charged - sum([i.timespent  for i in timecharged[p] if i.confirmed]) / (60.0*60.0)
            else:
                charged = 0.0
                unconfirmed = 0.0

            allocation = ompallocs[p].allocated/ (60.0*60.0)
            try:
                percentcomplete = 100.0 * charged/allocation
            except ZeroDivisionError:
                print('Warning: project {} has 0 hours allocated'.format(p))
                percentcomplete = np.nan


            # Observation information (all obs and bands)
            summary =  [i for i in summinfo if i.project == p]
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
            msbs = [i for i in msbinfo if i.project == p]
            faults = [i for i in faultinfo if i.project.lower() == p.lower()]
            f = faultsum(p, len([i for i in faults if i.status in [0,6]]),
                         len(faults))
            projdict[p] = projsum(summary, msbs, f, number_observations, time_observations,
                                            charged, unconfirmed, allocation, percentcomplete)
        webinfo['projdict'] = projdict

    else:
        webinfo['details'] = False

    projects = sorted(ompallocs, key = lambda name: ompallocs[name].priority)
    webinfo['projects'] = projects


    print('Creating completion pie chart')
    fig = Figure(figsize=(8,8))
    fig.patch.set_visible(False)
    axpie = fig.add_subplot(111)

    allocationspie = np.array([ompallocs[p].allocated for p in projects])
    observedtimes = np.array([ompallocs[p].allocated - ompallocs[p].remaining for p in projects])
    axpie = make_completion_summary_piechart(axpie, projects, allocationspie, observedtimes, radius=0.8)
    fig.tight_layout()
    canvas = FigureCanvas(fig)
    img = StringIO.StringIO()
    canvas.print_png(img)
    img.seek(0)
    webinfo['piechart'] = img.getvalue().encode('base64')

    print('creating completion charts')
    figs = OrderedDict()
    if (semester == 'LAP' or queue == 'LAP'):
        completionchartprojs = [LP_1, LP_2A, LP_2B]
    else:
        completionchartprojs = [projects]

    for projs in completionchartprojs:
        # Sort project by priority

        figure = Figure(figsize=(15,5))
        if exclude_done:
            allocs = { p: ompallocs[p] for p in projs  if p in ompallocs and ompallocs[p].remaining > 0}
        else:
            allocs = { p: ompallocs[p] for p in projs  if p in ompallocs}

        # Ensure projects are sorted in priority order.
        projs = sorted(allocs, key = lambda name: allocs[name].priority)

        # Exclude finished projects.
        charges = { p: timecharged[p] for p in projs if p in timecharged}

        fig = make_completion_chart(ompdb, figure, projs, allocs, charges,
                                    blocks=blocks)
        fig_url = url_for('completion_chart', proj=projs, blocks=int(blocks), inprog=int(exclude_done))
        # Create sendfile object
        canvas = FigureCanvas(fig)
        img = StringIO.StringIO()
        canvas.print_png(img)
        img.seek(0)
        figs[fig_url] = img.getvalue().encode('base64')
    webinfo['figs'] = figs


    if details:
        unknown_observations = {}
        print('creating project weatherband plots')
        projectfigs = {}
        for proj in projects:
            summary = [i for i in summinfo if i.project==proj]
            summary_unknown = [i for i in summary if i.band == 'unknown']
            if summary_unknown:
                unknown_observations[proj] = summary_unknown
            summary = [i  for i in summary if i.band != 'unknown']
            if summary != []:
                allocdict = allocations.get(proj, None)
                img = create_project_chart(proj, summary, allocdict=allocdict)
                projectfigs[proj] = img.getvalue().encode('base64')
        webinfo['projectfigs'] = projectfigs
        webinfo['unknownwvmobs'] = unknown_observations

    else:
        webinfo['projectfigs'] = {}


    if allprojectmsb:

        #If we already got msbinfo earlier don't redo the query.
        if not msbinfo:
            msbinfo = ompdb.get_summary_msb_info_group(projects=projects)

        # Remove all msbs where there is no time left in the project.
        msbinfo = [i for i in msbinfo if (i.project in ompallocs and ompallocs[i.project].remaining > 0)]
        remainingmsbs = {}
        for p in msbinfo:
            remainingmsbs[p.project] = ompdb.get_remaining_msb_info(p.project)

        # Mangle the msbinformation into a useful form.
        fullmsblist = []
        fullmsbdict = {}
        for p in remainingmsbs:
            projmsblist = []
            for m in remainingmsbs[p]:
                if m.pol == 1:
                    instrument = m.instrument + '-POL'
                else:
                    instrument = m.instrument
                taumin = max(m.taumin, ompallocs[p].taumin)
                taumax = min(m.taumax, ompallocs[p].taumax)
                # Tau values that are clearly meant to be in the range are not quite there...
                taumin = float('{:.2f}'.format(taumin))
                taumax = float('{:.2f}'.format(taumax))
                msb = MsbCustomInfo(p, instrument, m.target, m.title, m.coordstype, m.ra2000, m.dec2000,
                                    m.remaining, m.timeest, taumin, taumax, m.wavelength)
                fullmsblist.append(msb)
                projmsblist.append(msb)
            fullmsbdict[p] = projmsblist

        # Have to replace tau with actually available tau: set of allocation and msb setting
        instruments = set(i.instrument for i in fullmsblist)
        bdict = {}
        # Get all observations for given band
        bdict['1'] = [i for i in fullmsblist if i.taumin < Bands.BAND1[1]]
        bdict['2'] = [i for i in fullmsblist if (i.taumin < Bands.BAND2[1] and  i.taumax > Bands.BAND2[0])]
        bdict['3'] = [i for i in fullmsblist if (i.taumin < Bands.BAND3[1] and i.taumax > Bands.BAND3[0])]
        bdict['4'] = [i for i in fullmsblist if (i.taumin < Bands.BAND4[1] and i.taumax > Bands.BAND4[0])]
        bdict['5'] = [i for i in fullmsblist if i.taumax > Bands.BAND5[0]]
        bdict['extra'] = [i for i in fullmsblist if i not in bdict['1'] + bdict['2'] + bdict['3']+bdict['4'] + bdict['5']]
        webinfo['msbs'] = bdict

        # Create images.
        weatherfigs = {}
        for band in ['1','2','3','4','5','extra']:
            weatherfigs[band] = {}
            for instrument in instruments:
                msbs = [i for i in bdict[band] if i.instrument == instrument]
                if msbs:
                    fig = create_msb_image(msbs, utdate=datetime.date.today(),
                                           semesterdates=[datetime.date.today(),
                                                          datetime.date.today() + datetime.timedelta(days=6*31)],
                                           multiproject=True)
                    #fig.suptitle('Band {}: {}'.format(band, instrument))
                    canvas = FigureCanvas(fig)
                    img = StringIO.StringIO()
                    canvas.print_png(img)
                    img.seek(0)
                    figureobject = img.getvalue().encode('base64')
                    weatherfigs[band][instrument] = figureobject

            # If there are no values, remove.
            if set(weatherfigs[band].values()) == set([None]):
                print(band)
                weatherfigs[band] = None

        webinfo['msbfigs']=weatherfigs



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







def prepare_completion_chart(ompdb, semester='LAP', queue=None, patternmatch=None,
                             projects=None, exclude_projects=None, telescope='JCMT',
                             blocks=True,
                             exclude_done=True):
    """
    Create a plot showing the percentage completion of the large projects.

    Return sendfile object of mime-type image/png
    """
    if semester=='LP_1':
        semester = None
        projects = LP_1
    elif semester == 'LP_2A':
        print('Using LP_2A')
        semester = None
        projects = LP_2A
    elif semester == 'LP_2B':
        semester=None
        projects = LP_2B

    ompallocations = ompdb.get_allocations(semester=semester, queue=queue, projects=projects,
                                          patternmatch=patternmatch,
                                          telescope=telescope)

    timecharged = ompdb.get_time_charged_group(semester=semester, queue=queue, projects=projects,
                                               patternmatch=patternmatch, telescope=telescope)

    # Remove projects with no time remaining if requested.
    if exclude_done:
        for p in ompallocations:
            if ompallocations[p].remaining  <= 0:
                ompallocations.pop(p)
                timecharged.pop(p)

    # Remove excluded projects if provided.
    if exclude_projects:
        for p in exclude_projects:
            if p in ompallocations:
                ompallocations.pop(p)
            if p in timecharged:
                timecharged.pop(p)

    # Create matplotlib figure.
    figure = Figure(figsize=(15,5))
    fig = make_completion_chart(ompdb, figure, ompallocations.keys(),
                                ompallocations, timecharged,
                                blocks=blocks)

    # Create sendfile object
    canvas = FigureCanvas(fig)
    img = StringIO.StringIO()
    canvas.print_png(img)
    img.seek(0)
    return send_file(img, mimetype='image/png')




def make_completion_chart(ompdb, fig, projects, ompallocations, timecharged,
                          blocks=False, piechart=False):
    """
    Create plot showing completion of projects.

    Return a matplotlib figure

    Semester can be real semester, or LP1, LP2A or LP2B,
    Queue is the same as 'country': usually EC, DDT, PI, UH, LAP

    projects and exclude_projects can be lists.

    patternmatch needs to include the '%' wildcard symbols.
    """

    timechargedorig = timecharged
    ompallocsorig = ompallocations
    print('MAKE COMPLETION CHART: blocks=',blocks)

    # Find the last date being plotted.
    if timechargedorig:
        maxdate = max([item.date for sublist in timechargedorig.values() for item in sublist])
        mindate = min([item.date for sublist in timechargedorig.values() for item in sublist])
    else:
        maxdate = None
        mindate = None


    # Find out how many things we're cycling through, and don't use more than that number per figure.
    colorlen = len(matplotlib.rcParams['axes.prop_cycle'])
    if len(projects) > colorlen:
        # Need to put projects on multiple axes.
        import math
        number_plots = int(math.ceil(len(projects)/ colorlen))

        figsize = fig.get_size_inches()
        figsize_new = (figsize[0], figsize[1]*number_plots)
        fig.set_size_inches(figsize_new)
        # Need to divide up the projects.
        project_superlist = [projects[x:x+colorlen] for x in xrange(0, len(projects), colorlen)]

    else:
        number_plots = 1
        project_superlist = [projects]


    # Current date for title
    currdate = datetime.datetime.utcnow().strftime('%Y-%m-%d %H:%M')

    # Set up two plots on one figure.
    gs = gridspec.GridSpec(number_plots, 2, width_ratios=[4.5,2.5])

    for count in range(number_plots):

        # Set up axes, timecharged and ompallocs for this plot.
        ax = fig.add_subplot(gs[count*2])
        ax2 = fig.add_subplot(gs[(count*2)+1])

        projects = project_superlist[count]
        timecharged = { proj: timechargedorig[proj] for proj in projects if proj in timechargedorig}
        ompallocs = { proj: ompallocsorig[proj] for proj in projects }

        # Plot each projects time charged per day as cumulative percentage observed
        colordict = {}
        for p in projects:
            if p in timecharged:

                dates = [i.date for i in timecharged[p]]
                times = [i.timespent for i in timecharged[p]]
                cumultimes = [sum(times[0:i+1])for i,j in enumerate(times)]
                markevery=None

                # Add on a fake value so the graph continues to the
                # edge. BUT don't plot a marker as its not a real
                # observation.
                if dates[-1] != maxdate:
                    dates += [maxdate]
                    cumultimes += [cumultimes[-1]]
                    markevery = list(np.arange(len(cumultimes) - 1))
                allocation = ompallocs[p].allocated
                label = p

                lines = ax.plot(dates, [100* i/allocation for i in cumultimes],
                                label=label, markevery=markevery,
                                linewidth=2.0, markeredgewidth=2.0, marker='x')
                colordict[p] = lines[0].get_color()
            else:
                color = next(ax._get_lines.prop_cycler)['color']
                colordict[p] = color

        # set up plot limits and labels
        ylim=ax.get_ylim()
        ax.set_ylim(0, max(ylim[1], 100))
        #ax.set_xlim(sdate,edate)
        ax.set_ylabel('Percentage completion')
        fig.patch.set_visible(False)
        xlim = ax.get_xlim()
        ax.set_xlim(mindate, maxdate)

        ax.xaxis.set_minor_locator(matplotlib.dates.MonthLocator())

        ax.yaxis.grid()

        # Get times of scheduled blocks.
        if blocks:
            lapfile = 'LAP-UT-blocks-new.txt'
            if os.path.isfile(lapfile):
                print('found LAP file')
                f = open('LAP-UT-blocks-new.txt', 'r')
                blocks = f.readlines()
                f.close()
                blocks = [i.strip() for i in blocks]
                blocks = [i.split() for i in blocks]
                for start, end, project in blocks:
                    if project.upper() in projects:
                        start = datetime.datetime(int(start[0:4]), int(start[4:6]), int(start[6:8]), 0, 0)
                        end =  datetime.datetime(int(end[0:4]), int(end[4:6]), int(end[6:8]), 23, 59)
                        ax.axvspan(start, end, alpha=0.25, facecolor=colordict[project.upper()], edgecolor='None')


        # Set up second plot: bar chart of percentage observed.
        observedtimes = []
        allocations = []
        for p in projects:
            if p in timecharged:
                t = timecharged[p]
                observedtimes.append(sum([i.timespent for i in t]))
            else:
                observedtimes.append(0.0)
            allocations.append(ompallocs[p].allocated)


        total = [100 for i in observedtimes]
        percent_observed = np.array([100*j/allocations[i] if allocations[i] != 0 else np.nan
                            for i,j  in enumerate(observedtimes) ])

        labels = []
        for i in range(len(projects)):
            enabled = bool(ompallocs[projects[i]].enabled)
            if not enabled:
                symbol = u'\u2717'
            else:
                symbol = ' '
            hourspent = observedtimes[i]/(60.0*60.0)
            allochrs = allocations[i]/(60.0*60.)
            if allochrs != 0:
                percent = 100 * hourspent/allochrs
            else:
                percent = np.nan
            textstring = symbol + ' {}: {:.1F}/{:.1F} ({:.3G}%)'.format(projects[i],
                                                                           hourspent, allochrs,
                                                                           percent)#, ompallocs[projects[i]].title)
            labels.append(textstring)

        if piechart:
            projlabels = [i.split(':')[0].strip() for i in labels]
            ax2 = make_completion_summary_piechart(ax2, projects, allocations, observedtimes,
                                                   projectlabels=projlabels, radius=0.8, percentages=True)

        else:
            # Plot a bar for each project.
            ax2.barh([i+0.05 for i in range(len(projects))], total, label='allocated',
                     alpha=0.3, color=[colordict[p] for p in projects],
                     edgecolor='white', height=0.9)

            ax2.barh([i+0.05 for i in range(len(projects))],percent_observed, label='observed',
                     alpha=1.0, color=[colordict[p] for p in projects], edgecolor='white',
                     height=0.9)

            # Label each bar.
            #labels = []

            for i in range(len(projects)):
                textstring = labels[i]
                ax2.text(5, i-0.2, textstring,
                         color='black', ha='left', alpha=1.0, va='top', weight='bold', size='medium')
                ax2.text(5, i+0.4, ompallocs[projects[i]].title,
                         color='black', ha='left', alpha=1.0, va='bottom', size='small')
            ax.text(0.1, 0.9, 'Completion as of ' + currdate, transform=ax.transAxes)


            # Finish up plot
            ax2.set_axis_off()
            ax2.invert_yaxis()

    fig.autofmt_xdate()

    if (number_plots > 1):
        ax = fig.axes[0]
        ax.tick_params(axis='x', which='major', top='on', bottom='on')
        ax.tick_params(axis='x', which='major',labeltop='on', labelbottom=False)
        majticklabels_shown = fig.axes[-2].get_xmajorticklabels()
        majticklabels = ax.get_xmajorticklabels()
        angle = majticklabels_shown[0].get_rotation()
        for i in majticklabels:
            i.set_visible(True)
            i.set_rotation(angle)
            i.set_ha('left')
    fig.set_tight_layout(tight=True)

    return fig





def prepare_project_chart(project,summary, allocdict=None):
    img = create_project_chart(project, summary, allocdict=allocdict)
    return send_file(img, mimetype='image/png')

def create_project_chart(project, summary, allocdict=None):
    """
    Create chart for project obs.
    """

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

        color = None
        ax.bar([i+0.05 for i in range(0,5)],
               gooddata, color=color,
                width=0.9)
        ax.bar([i+0.05 for i in range(0,5)],
                questdata,
                bottom=gooddata,
                color = 'orange',
                width=0.9)

        ax.set_xticks([i for i in range(0,5)])
        xlabels = []
        for i in range(0,5):
            label = 'B%i'%(i+1)
            xlabels.append(label)

        ax.set_xticklabels(xlabels, weight='bold')

        if allocdict:
            for i in range(0,5):
                allocation = allocdict.get(instrument, {}).get(i+1,0)
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
                ax.text(i*0.2 + 0.1, 0.1, label, transform=ax.transAxes,
                        color='black', ha='center', alpha=1.0, va='center', weight='bold', rotation='horizontal')

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
    return img



def make_completion_summary_piechart(ax, projects, allocations, observedtimes, projectlabels=None,
                                     colors=None, cmap='viridis_r',
                                     radius=1.0, percentages=None, radial_labels=True, scale=1.5):

    # Label with projects if no specific project labels.
    if not projectlabels:
        projectlabels = projects

    if colors is None:
        cmap = matplotlib.cm.get_cmap(cmap)
        colors = cmap(np.linspace(0., 1., int(scale*len(projects))))
        lower = int(np.floor(0.5*(int(scale*len(projects)) -len(projects))))
        colors = colors[lower:lower+len(projects)]

    percent_observed = 100 * np.array(observedtimes)/np.array(allocations)



    wedges1, texts, autolabels = ax.pie(allocations, labels=projectlabels,
                                         autopct='%.1F', startangle=90.0, counterclock=False,
                                         labeldistance=1.1, pctdistance=0.8*radius, radius=radius, colors=colors)
    wedges2, _ = ax.pie(allocations, startangle=90.0, counterclock=False, radius=radius)

    # Go through wedges, set color, alpha and size appropriately.
    for i in range(len(projects)):
        color = wedges1[i].get_facecolor()
        #texts[i].set_color(color)
        if percentages:
            autolabels[i].set_text('{:.1F}%'.format(percent_observed[i]))
        else:
            autolabels[i].set_visible(False)

        wedges1[i].set_alpha(0.2)
        wedges2[i].set_facecolor(color)
        # Set radius so that fraction of area filled in is proportional to copmletion
        wedges2[i].set_radius(radius*np.sqrt(percent_observed[i]/100.0))

    # Add a vertical black line at start, and an arrow indicating priority order.
    lines = ax.plot([0,0], [0,1.05*radius], color='black', linewidth=2.0)
    posA = (0, 1.05*radius)
    posB =(1.05*radius * np.cos(np.deg2rad(90 - 10.0)), 1.05*radius*np.sin(np.deg2rad(90 - 10.0)))
    arr = ax.annotate('', xy=posB, xycoords='data', xytext=posA, textcoords='data',
                          arrowprops=dict(arrowstyle="simple", fc='black', ec="none",
                                          connectionstyle="arc3, rad=-0.1"))

    # Pie charts should be circular!
    ax.set_aspect('equal')

    # Ensure labels are pointing outwards.
    if radial_labels:
        for t, a in zip(texts, autolabels):
            position = t.get_position()
            angle = np.rad2deg(np.arctan(position[1]/position[0]))
            if position[0] > 0 and position[1] < 0:
                t.set_ha('left')
                t.set_va('baseline')
            if position[0] > 0 and position[1] > 0:
                t.set_ha('left')
                t.set_va('bottom')
            if position[0] < 0 and position[1] < 0:
                t.set_ha('right')
                t.set_va('baseline')
            if position[0] < 0 and position[1] > 0:
                t.set_ha('right')
                t.set_va('bottom')
            if angle > 90:
                angle += 180
            elif angle < -90:
                angle += 180
            t.set_rotation(angle)
            a.set_rotation(angle)

    for t, a in zip(texts, autolabels):
        t.set_size('small')
        a.set_size('small')

    # Add faint black lines around the patches
    for w in ax.patches:
        w.set_edgecolor((0.0,0.0,0.0,0.5))
        w.set_linewidth(1.0)

    return ax
