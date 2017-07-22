# Copyright (C) 2016 East Asian Observatory
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

from flask import render_template
from collections import namedtuple
import numpy as np
from jsa_proc.admin.directories import get_output_dir
from jsa_proc.error import NoRowsError
from omp.error import OMPDBError
import werkzeug.exceptions
from astropy.time import Time as aTime
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy import units as u
import datetime
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from flask import send_file
from io import BytesIO
import matplotlib

import base64

ObsInfo = namedtuple('ObsInfo', 'obsid obsnum '+ \
                     'instrument project time obstype source '\
                     + 'tau transmission seeing '\
                     + 'comment hetinfo jobinfo')
Time = namedtuple('Time',
                  'start end duration')

Source = namedtuple('Source',
                    'name moving_target standard')
Tau = namedtuple('Tau',
                 'average start end start_offset end_offset source')
Transmission = namedtuple('Tranmission',
                          'average start end starttime endtime wvmsource')
Comment = namedtuple('Comment',
                     'status text author datetime')
HetInfo = namedtuple('HetInfo',
                     'restfreq molecule transition velocity velosys bwmode subsysnr obsidss')

JobInfo = namedtuple('JobInfo',
                     'job_id obsid obsidss state preview')

def get_obstype(obs):
    if obs.obs_type != 'science':
        return obs.obs_type.capitalize()
    else:
        if obs.instrume == 'SCUBA-2':
            obstype =  obs.scan_pat
            if obstype=='CV_DAISY':
                return 'Daisy'
            elif obstype=='CURVY_PONG':
                return 'Pong-{}'.format(int(obs.map_wdth))
        else:
            if obs.sam_mode in ['scan', 'raster']:
                return 'Raster'
            else:
                return obs.sam_mode.capitalize() + '-' + obs.sw_mode.capitalize()

    # outputs should be one of
    # pointing, CV_DAISY, CURVY_PONG, #SCUBA-2 or obs
    # If obs_type != 'science': return 'obs_type'
    # else if science & scuba-2: return scan_pat
    # else if ACSIS & science:
    #   if sam_mode=scan, or raster: return scan_pat (but replace DISCRETE_BOUSTROPHEDON with Raster)
    #   if sam_mode = 'jiggle'  or 'grid': return sam_mode

    #
def create_obsinfo(obsinfo, acsisinfo, procjobs, jsaprocdb):
    """Create obsinfo objects for every observation.

    Return a dictionary, with each key being the name of the instrument,
    and the value being the list of observations.

    """
    # Turn acsisinfo into dict
    aobsids = set(a.obsid for a in acsisinfo)
    acsis = {}
    for o in aobsids:
        acsis[o] = [a for a in acsisinfo if a.obsid==o]

    # Get SCUBA-2 jobs if relevant
    scuba_lookupdict={}

    for a in [i.obsid for i in obsinfo if i.instrume=='SCUBA-2']:
        match850 = 'jcmt_{}_reduced-850um_preview_64.png'.format(a.lower())
        match450 = 'jcmt_{}_reduced-450um_preview_64.png'.format(a.lower())
        info850 = [(j, i) for j in procjobs if j.outputs is not None for i in j.outputs if i==match850]
        info450 = [(j, i) for j in procjobs if j.outputs is not None for i in j.outputs if i==match450]
        scuba_lookupdict[a] = {'450': info450, '850': info850}

    # Get acsis jobs if relevant
    acsis_lookupdict = {}
    for a in acsisinfo:

        match1 = 'jcmt_{}'.format(a.obsid.lower())
        match2 = '{}_preview_64.png'.format(a.subsysnr)

        info = [(j, i) for j in procjobs if j.outputs is not None for i in j.outputs if i.startswith(match1) and i.endswith(match2)]

        if not a.obsid in acsis_lookupdict:
            acsis_lookupdict[a.obsid] = {a.subsysnr: info}
        else:
            acsis_lookupdict[a.obsid][a.subsysnr] = info

    for o in aobsids:
        acsis[o] = [a for a in acsisinfo if a.obsid==o]


    output = []
    for o in obsinfo:
        time = Time(o.date_obs, o.date_end, o.date_end - o.date_obs)
        source = Source(o.object, bool(o.moving_target), bool(o.standard))
        if o.wvmtaust and o.wvmtauen and o.wvmdatst and o.wvmdaten:
            tau = Tau(np.mean([o.wvmtaust,o.wvmtauen]), o.wvmtaust, o.wvmtauen,
                o.date_obs - o.wvmdatst, o.date_end - o.wvmdaten, 'WVM')
        elif o.tau225st and o.tau225en and o.taudatst and o.taudaten:
            tau = Tau(np.mean([o.tau225st,o.tau225en]), o.tau225st, o.tau225en,
                      o.date_obs - o.taudatst, o.date_end - o.taudaten, o.tausrc)
        else:
            tau = None
        comment = Comment(o.commentstatus, o.commenttext, o.commentauthor,
                          o.commentdate)
        hetinfos = []
        if o.backend=='ACSIS':
            for a in acsis[o.obsid]:
                obsidss = o.obsid + '_' + str(a.subsysnr)
                hetinfos.append(HetInfo(a.restfreq, a.molecule, a.transition,
                                        a.zsource*3e8*1e-3,
                                        a.doppler, a.bwmode, a.subsysnr, obsidss))
            hetinfos.sort(key=lambda x: x.subsysnr)


        output.append(ObsInfo(o.obsid, o.obsnum, get_instrument(o), o.project, time,
                              get_obstype(o), source, tau, None,
                              None, comment, hetinfos, None))

    outputdict = {}
    instruments = set(o.instrument for o in output)
    for i in instruments:
        values = [o for o in output if o.instrument==i]
        values.sort(key=lambda o: o.time.start)
        outputdict[i] = values

    lookupdict = {}
    lookupdict.update(acsis_lookupdict)
    lookupdict.update(scuba_lookupdict)

    return outputdict, lookupdict

def get_previews_for_observation_from_job(obsid, acsisvals, job_id):
    try:
        outputs = jsaprocdb.get_output_files(job_id)
        allpreview = [a for a in outputs if a.endswith('64.png')]

        preview = [a for a in allpreview if a.startswith('jcmt_{}'.format(obsid.lower())) and 'reduced' in a]

    except:
        preview = None
    return preview

def prepare_observation_page(ompdb, jsaprocdb, projectid, utdatestart=None, utdateend=None,
                             ompstatus=None):
    try:
        projinfo = ompdb.get_project_info(projectid)

    except OMPDBError:
        raise werkzeug.exceptions.InternalServerError(
            'Project {} was not found in the database.'.format(projectid))

    if ompstatus:
        ompstatus = int(ompstatus)
    if utdatestart:
        utdatestart = int(utdatestart)
    if utdateend:
        utdateend = int(utdateend)

    obsinfo = ompdb.get_observations(projectid, utdatestart=utdatestart,
                                     utdateend=utdateend, ompstatus=ompstatus)

    procjobs = jsaprocdb.find_jobs(task=['jcmt-nightly'], obsquery={'project': projectid},
                                   outputs='jcmt\_%\_0%\_%\_reduced-%64.png')

    if obsinfo:

        if 'ACSIS' in set(i.backend for i in obsinfo):
            acsisinfo = ompdb.get_acsis_info(projectid)
        else:
            acsisinfo = []
        obsinfodict, lookupdict = create_obsinfo(
                obsinfo, acsisinfo, procjobs, jsaprocdb)
        obsids = set([o.obsid for o in obsinfo])
        extracss = create_css(obsids)
    else:
        obsinfodict={}
        extracss = ''
        lookupdict={}

    return render_template('observations.html', projinfo=projinfo,
                           title=projectid, obsinfodict=obsinfodict, ompstatus=ompstatus,
                           utdatestart=utdatestart, utdateend=utdateend,
                           lookupdict=lookupdict)

def prepare_project_page(ompdb, jsaprocdb, projectid, observations=False):
    s = datetime.datetime.now()
    try:
        projinfo = ompdb.get_project_info(projectid)

    except OMPDBError:
        raise werkzeug.exceptions.InternalServerError(
            'Project {} was not found in the database.'.format(projectid))

    title = '{}: {}'.format(projinfo.id, projinfo.title)

    if observations:
        obsinfo = ompdb.get_observations(projectid)


        if obsinfo:
            procjobs = jsaprocdb.find_jobs(task=['jcmt-nightly'], obsquery={'project': projectid},
                                   outputs='jcmt\_%\_0%\_%\_reduced-%64.png')
            if 'ACSIS' in set(i.backend for i in obsinfo):
                acsisinfo = ompdb.get_acsis_info(projectid)
            else:
                acsisinfo = []
            obsinfodict, acsis_lookupdict, scuba_lookupdict = create_obsinfo(
                obsinfo, acsisinfo, procjobs, jsaprocdb)

            obsids = set([o.obsid for o in obsinfo])
            extracss = create_css(obsids)
        else:
            obsinfodict={}
            extracss = ''
    else:
        extracss = ''
        obsinfodict = {}


    msbfullinfo = ompdb.get_remaining_msb_info(projectid)
    time_in_msbs = (sum([i.remaining*i.timeest for i in msbfullinfo]))
    # Turn into dictionary by instrument.
    instruments = set((m.instrument, m.pol) for m in msbfullinfo)
    msbresults = {}
    for i in instruments:
        if i[1]==1:
            instrument = i[0] + '-POL'
        else:
            instrument = i[0]
        msbresults[instrument] = [m for m in msbfullinfo if (m.instrument==i[0] and m.pol==i[1])]

    completion_percentage = 100.0 * (projinfo.allocated_hours - projinfo.remaining_hours) / projinfo.allocated_hours
    mtime = datetime.datetime.now()
    print('MSB fullinfo {}'.format((mtime -s).total_seconds()))

    # Create custom css for observation table
    #obssummary = ompdb.get_summary_obs_info(projectid)
    obssummary = None
    otime = datetime.datetime.now()

    #print('Obssummary time {}'.format((otime - mtime).total_seconds()))
    #fullobs = ompdb.get_observations(projectid)
    otime2 = datetime.datetime.now()
    #print('Obs full time {}'.format((otime2 - otime).total_seconds()))

    #project_chart = prepare_project_chart(projectid, obssummary)
    project_chart=None
    fig = create_msb_image(msbfullinfo, datetime.date.today(),
                     (datetime.datetime.strptime('2016-08-01', '%Y-%m-%d').date(),
                      datetime.datetime.strptime('2017-02-01', '%Y-%m-%d').date()))
    canvas = FigureCanvas(fig)
    msbimg= BytesIO()
    canvas.print_png(msbimg)
    msbimg.seek(0)
    msbpng = base64.b64encode(msbimg.getvalue())
    ctime = datetime.datetime.now()
    print('Figures time {}'.format((ctime - otime2).total_seconds()))
 
    return render_template('project.html', title=title,projinfo=projinfo, obsinfodict=obsinfodict,
                           msbfullinfo=msbfullinfo, time_in_msbs=time_in_msbs,
                           completion_percentage=completion_percentage, extracss=extracss,
                           msbresults=msbresults, obssummary=obssummary, project_chart=project_chart,
                           msbpng = msbpng)

def create_css(obsids):
    cssselector = []
    for o in obsids:
        cssselector.append('tr.{} + tr.{} td.het-additional'.format(o,o))
    cssselector = ',  '.join(cssselector)
    cssrule = '{}{{color:gray; font-style:italic;}}'.format(cssselector)
    return cssrule

def get_instrument(o):
    if o.instrume == o.backend:
        instrument = o.instrume
    else:
        instrument = '-'.join([o.instrume, o.backend])

    if instrument == 'SCUBA-2' and o.inbeam and 'pol' in o.inbeam.lower():
        instrument = 'POL-2'
    else:
        if o.inbeam:
            instrument = '-'.join([o.inbeam, instrument])
    return instrument



def create_msb_image(msbs, utdate, semesterdates, multiproject=False):
    """
    Create a plot of observability ala source plot, for all sources in msb list.

    msbs: list of msbinfo objects
    figure: matplotlib figure object.
    utdate: datetime.date object
    semesterdates: start and end tuples for semester.

    Returns a figure object.
    """
    fig = Figure(figsize=(12,5))
    fig.set_facecolor([0.7,0.7,0.7,0.0])
    # Get telescope position.
    jcmt = EarthLocation(lat=19.82283890588*u.degree, lon=-155.4770278387 *u.degree, height=4120.0*u.meter)

    #get time
    utcoffset = -10*u.hour # HST time
    time = utdate.strftime('%Y-%m-%d 0:00:00') # Today

    midnight_hi = aTime(time) - utcoffset
    delta_midnight = np.linspace(-12,12,100)*u.hour
    frame_tonight = AltAz(obstime=midnight_hi + delta_midnight, location=jcmt)

    # semester stuff
    start=aTime(semesterdates[0].strftime('%Y-%m-%d'))
    end = aTime(semesterdates[1].strftime('%Y-%m-%d'))
    delta = end - start
    semtimes = start + np.linspace(0, delta.value-1, delta.value) * u.day
    # Get Coordinate info
    coordstypes = set([i.coordstype for i in msbs])
    plotdict={}
    coorddict={}


    # First plot: observability at requested night.
    ax = fig.add_subplot(121)

    for coord in coordstypes:
        if coord != 'RADEC':
            print('Warning: non-RA-DEC coordinates not yet supported')
        else:
            ra = [i.ra2000 for i in msbs if i.coordstype==coord]
            dec = [i.dec2000 for i in msbs if i.coordstype==coord]
            if not multiproject:
                labels = [i.target for i in msbs if i.coordstype==coord]
            else:
                labels = ['{}: {}'.format(i.project, i.target) for i in msbs if i.coordstype==coord]
                projects = [i.project for i in msbs if i.coordstype==coord]
                projectcolors = {}
                for p in set(projects):
                    projectcolors[p] = next(ax._get_lines.prop_cycler)['color']
                colors = [projectcolors[i.project] for i in msbs if i.coordstype==coord]

            coords = SkyCoord(ra=np.rad2deg(ra)*u.degree, dec=np.rad2deg(dec)*u.degree, frame='fk5')
            coorddict[coord] = coords
            sources_tonight = coords[:, np.newaxis].transform_to(frame_tonight)
            plotdict[coord] = sources_tonight, labels




    for coord, labels in plotdict.values():
        times = np.array([delta_midnight.value]*(len(coord.alt.value))).swapaxes(0,1) * u.hour
        pcoords = coord.alt.value.swapaxes(0,1)

        lines = ax.plot(times, pcoords)

        if multiproject:
            for l, c in zip(lines, colors):
                l.set_color(c)

        peak_alts = coord.alt.value.max(axis=1)
        peak_times = delta_midnight.value[coord.alt.value.argmax(axis=1)]
        for a,t,la, li in zip(peak_alts, peak_times, labels, ax.lines):
            ax.text(t,a,la, color=li.get_color(), zorder=100)


    ax.set_xlim(-12,12)
    xticks = np.array(ax.get_xticks())
    xticks[xticks < 0] = xticks[xticks<0] + 24

    ax.set_xticklabels(['{}'.format(int(i)) for i in xticks])
    ax.set_ylim(0,90)
    ax.grid()
    ax.set_xlabel('Time (HST)')
    ax.set_ylabel('Altitude')
    ax.set_title('Observability at {}'.format(utdate.strftime('%Y-%m-%d')))
    ax.hlines(30.0, -12, 12)
    ax.fill_betweenx([0,90], [18.5-24, 18.5-24.0], [6.50, 6.5], color='0.7', alpha=0.5)
    ax.fill_betweenx([0,90], [6.5, 6.5], [12.0, 12.0], color='0.7', alpha=0.2)

    # Second figure: observability over semester
    if 'RADEC' in coorddict:
        c = coorddict['RADEC']
        c=c[:, np.newaxis]

        semtimeshst  = semtimes - utcoffset
        transits = (24 - (semtimeshst.sidereal_time('mean', longitude=jcmt.lon).value  - c.ra.hourangle )) % 24
        transits[ transits > 12] -= 24
        # Prevent wrapping
        for i in transits:
            i[i.argmax()]=np.nan

        ax2 = fig.add_subplot(122)
        times=np.array([semtimeshst.datetime]*(len(c.ra.value))).swapaxes(0,1)
        ptransits = transits.swapaxes(0,1)
        lines = ax2.plot(ptransits, times)
        if multiproject:
            for l, c in zip(lines, colors):
                l.set_color(c)

        loc=matplotlib.dates.WeekdayLocator(byweekday=1, interval=2)
        ax2.yaxis.set_major_locator(loc)
        ax2.yaxis.set_major_formatter(matplotlib.dates.DateFormatter(fmt='%Y-%m-%d'))
        ax2.yaxis.tick_right()
        ax2.set_ylim(times.max(), times.min())

        ax2.set_xlim(-12, 12)
        xticks = np.array(ax2.get_xticks())
        xticks[xticks < 0] = xticks[xticks<0] + 24
        ax2.set_xticklabels(['{}'.format(int(i)) for i in xticks])

        ax2.grid()
        ax2.minorticks_on()
        ax2.set_xlabel('Time (HST)')
        ax2.set_ylabel('Date')
        ax2.set_title('Time of transits {} to {}'.format(semesterdates[0].strftime('%Y-%m-%d'),
                                                   semesterdates[1].strftime('%Y-%m-%d')))
        ax2.fill_betweenx([times.max(),times.min()], [18.5-24, 18.5-24.0], [6.50, 6.5], color='0.7', alpha=0.5)
        ax2.fill_betweenx([times.max(), times.min()], [6.5, 6.5], [12.0, 12.0], color='0.7', alpha=0.2)
        ax2.hlines(utdate, -12, 12)
    fig.set_tight_layout(True)

    return fig



# Fields we care about:
# instrumnet( instrume, backend, inbeam)
# Date (from date_obs)
# duration (date_end - date_obs)
# source (object)
# project (not needed here, but will be sometimes)
# average tau and method (wvm, cso, warn if time of wvm > 10 minutes from start and end of observation)
# Transmission? -- tau, 
# seeing? if have
# scan mode -- scan_pat, sam_mode, sw-mode
# status
# last comment
# moving_target ? (not right, but should fix)
