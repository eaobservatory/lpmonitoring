# -*- coding: utf-8 -*-
# Copyright (C) 2015 East Asian Observatory
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
import datetime
import functools
import os
from flask import Flask, redirect, url_for, render_template, send_file, request

from summary import create_summary, get_omp_database, prepare_completion_chart, prepare_project_chart

from project import prepare_project_page, create_msb_image, prepare_observation_page

from omp.obs.state import OMPState
from omp.error import OMPError
from jsa_proc.config import get_database as get_jsaproc_database
from jsa_proc.admin.directories import get_output_dir
from jsa_proc.state import JSAProcState
from matplotlib.figure import Figure
import StringIO
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

from tag_allocation_summary import prepare_tagallocation_summary

def web_pages():
    """ Create web pages for large program monitoring

    These are designed to allow easy monitoring by the
    observatory of the completion of these projects.

    Import notes: completion rates, MSBs available,
    state of observations, keeping track of observations observed
    in the wrong weather band.

    This list will  likely expand in general.
    """

    app = Flask('lpmonitoring')
    ompdb = get_omp_database()
    jsaprocdb = get_jsaproc_database()


    # Setup root page to redirect to summary
    @app.route('/')
    def home_page():
        return redirect(url_for('summary', semester='LAP', inprog=0), code=303)

    @app.route('/project/<projectid>/observations')
    def observation(projectid):
        projectid = str(projectid)
        projectid=projectid.replace('-','/')
        utdatestart = request.args.get('utdatestart', None)
        utdateend = request.args.get('utdateend', None)
        ompstatus = request.args.get('ompstatus', None)
        if ompstatus:
            try:
                ompstatus = OMPState.lookup_name(ompstatus)
            except OMPError:
                ompstatus=ompstatus
        return prepare_observation_page(ompdb, jsaprocdb, projectid,
                                        utdatestart=utdatestart,
                                        utdateend=utdateend,
                                        ompstatus=ompstatus)

    @app.route('/project/<projectid>/msbchart')
    def msbchart(projectid):
        msbs = ompdb.get_remaining_msb_info(str(projectid).replace('-','/'))

        utdate = request.args.get('utdate', None)
        if utdate:
            utdate = datetime.datetime.strptime(utdate, '%Y-%m-%d').date()
        else:
            utdate = datetime.date.today()
        semstart = request.args.get('semstart', None)
        if semstart:
            semstart = datetime.datetime.strptime(semstart, '%Y-%m-%d').date()
        else:
            semstart = datetime.datetime.now().date()

        semend = request.args.get('semend', None)
        if semend:
            semend = datetime.datetime.strptime(semend, '%Y-%m-%d').date()
        else:
            semend = semstart + datetime.timedelta(days=6*31)
        fig = create_msb_image(msbs, utdate, (semstart, semend))
        canvas = FigureCanvas(fig)
        img = StringIO.StringIO()
        canvas.print_png(img)
        img.seek(0)
        return send_file(img, mimetype='image/png')

    # Summary page
    @app.route('/summary/')
    def summary():
        semester = request.args.get('semester', None, type=str)
        queue = request.args.get('queue', None, type=str)
        patternmatch = request.args.get('patternmatch', None, type=str)
        if patternmatch:
            patternmatch = '%{}%'.format(patternmatch)
        projects = request.args.getlist('proj', type=str)
        exclude_projects = request.args.getlist('notproj', type=str)
        #telescope = request.args.get('telescope', 'JCMT', type=str)
        blocks = bool(request.args.get('blocks', 1, type=int))
        exclude_done = bool(request.args.get('inprog', 0, type=int))
        details = bool(request.args.get('details', 1, type=int))
        if not semester and not queue and not patternmatch and not projects:
            semester='LAP'
        values = create_summary(ompdb, semester=semester, queue=queue, patternmatch=patternmatch,
                                projects=projects, exclude_projects=exclude_projects, exclude_done=exclude_done,
                                details=details, blocks=blocks)
        return render_template('summary.html', **values)

    # Chart showing overall completion by project.
    @app.route('/summary/completionchart/')
    def completion_chart():
        semester = request.args.get('semester', None, type=str)
        queue = request.args.get('queue', None, type=str)
        patternmatch = request.args.get('patternmatch', None, type=str)
        if patternmatch:
            patternmatch = '%{}%'.format(patternmatch)
        projects = request.args.getlist('proj', type=str)
        exclude_projects = request.args.getlist('notproj', type=str)
        telescope = request.args.get('telescope', 'JCMT', type=str)
        blocks = bool(request.args.get('blocks', 0, type=int))
        exclude_done = bool(request.args.get('inprog', 1, type=int))
        if not patternmatch and not semester and not queue and not projects:
            semester='LAP'
        return prepare_completion_chart(ompdb, semester=semester, queue=queue, patternmatch=patternmatch,
                                        projects=projects, exclude_projects=exclude_projects, telescope=telescope,
                                        blocks=blocks, exclude_done=exclude_done)


    @app.route('/tag-summary/<semester>')
    def tag_summary(semester):
        values =  prepare_tagallocation_summary(semester)
        return render_template('tag-summary-semester.html', **values)


    @app.route('/summary/summary_query', methods=['POST'])
    def summary_query():

        semester = request.form.get('semester', None, type=str)
        if semester:
            semester = semester.upper()
        queue = request.form.get('queue', None, type=str)
        if queue== '':
            queue = None
        if queue:
            queue = queue.upper()
        patternmatch = request.form.get('patternmatch', None, type=str)
        if patternmatch == '':
            patternmatch = None
        projects = request.form.getlist('projects', type=str)
        if projects == [] or projects == ['']:
            projects = None
        else:
            projects = set(projects)
        exclude_projects = request.form.getlist('exclude_projects', type=str)
        if exclude_projects == [] or exclude_projects == ['']:
            exclude_projects = None
        else:
            exclude_projects = set(exclude_projects)
        blocks = request.form.get('blocks', 0, type=int)
        details = request.form.get('details', 0, type=int)
        exclude_done = request.form.get('exclude_done', 0, type=int)
        print(details, blocks, exclude_done)
        fopid = request.form.get('fopid', None, type=str)
        if fopid:
            fopid = fopid.upper()
            projects_fopid = set(ompdb.get_support_projects(fopid, semester))
            print(projects_fopid)
            if projects:
                projects = projects_fopid.intersection(projects)
            else:
                projects = projects_fopid
        if projects:
            projects = list(projects)
        if exclude_projects:
            exclude_projects = list(exclude_projects)
        return redirect(url_for('summary', semester=semester, queue=queue, patternmatch=patternmatch, proj=projects,
                                notproj=exclude_projects, blocks=blocks, inprog=exclude_done, details=details))
    @app.route('/summary/chart/<projectid>/')
    def project_chart(projectid):
        summary = ompdb.get_summary_obs_info(str(projectid))
        return prepare_project_chart(projectid, summary)


    @app.route('/fop')
    def fop():
        return render_template('fop-form.html')

    @app.route('/fop_summary_getres', methods=['POST'])
    def fop_summary_getres():
        fopompid = request.form['fopompid']
        semester = request.form['semester']
        return redirect(url_for('fop_listing', fopuname=fopompid, semester=semester))


    @app.route('/fop/<fopuname>/<semester>')
    def fop_listing(fopuname, semester):
        fopuname = str(fopuname).upper()
        semester = str(semester).upper()
        projects = ompdb.get_support_projects(fopuname, semester)
        summarydict = {}
        for p in projects:
            projinfo = ompdb.get_project_info(p)
            obsinfo = ompdb.get_summary_obs_info(p)
            obsdict = {}
            for s in OMPState.STATE_ALL:
                number = sum([i.number for i in obsinfo if i.status==s])
                totaltime = sum([i.totaltime for i in obsinfo if i.status==s])
                obsdict[OMPState.get_name(s).lower()] = [number, totaltime]


            msbinfo = ompdb.get_summary_msb_info(p)
            faultinfo = ompdb.get_fault_summary(p)
            faultinfo = [len([i for i in faultinfo if i.status in [0,6]]), len(faultinfo)]
            summarydict[p] = projinfo, obsdict, msbinfo, faultinfo
        return render_template('fop-semester-summary.html', summary=summarydict,
                               fopuname=fopuname, semester=semester)

    @app.route('/project/<projectid>')
    def project_page(projectid):
        projectid = str(projectid)
        projectid=projectid.replace('-', '/')
        print('Projectid is {}'.format(projectid))
        return prepare_project_page(ompdb, jsaprocdb, projectid)

    @app.route('/jsaproc_preview/<job_id>/<preview>')
    def job_preview(job_id, preview):
        path = os.path.join(get_output_dir(int(job_id)), preview.lower())
        return send_file(str(path), mimetype='image/png')

    # Get all the bands for which any data is allocated from an
    # allocation dictionary.
    @app.template_filter('allocationbandset')
    def allocationbandset(allocdict):
        return set([item for sublist in allocdict.values() for item in sublist.keys()])

    # A filter to get the set of an attributes values from a list
    @app.template_filter('getset')
    def getset(listofthings, attribute):
        return set([getattr(i,attribute) for i in listofthings])

    @app.template_filter('contactable')
    def contactable(userinfo):
        if not userinfo.contactable:
            return u'<td title="Does not receive notifications.">✗</td>'
        else:
            return u'<td title="' + userinfo.email +u'">✉</td>'
    @app.template_filter('max')
    def mymax(alist):
        themax = max(alist)
        return themax
    @app.template_filter('cadc')
    def cadc(userinfo):
        if not userinfo.cadcuser:
            return u'<td title="No CADC username">✗</td>'
        else:
            return u'<td title="CADC: ' + userinfo.cadcuser +u'">✓</td>'


    @app.template_filter('acsispreview')
    def get_acsis_jobid_preview(hetinfo, o, jobinfo):
        # find matching jobinfo
        print('Hetinfo is ', hetinfo)
        print('job info is ', jobinfo)
        print('obsinfo is ', o)
        jobinfos = [j for j in jobinfo if int(j.obsidss.split('_')[-1])==int(hetinfo.subsysnr)]

        if jobinfos:
            jobinfos = jobinfos[0]
            job_id = jobinfos.job_id
            preview = jobinfos.preview
        else:
            job_id=None
            preview = None
        return job_id, preview


    @app.template_filter('ompstatus')
    def ompstatus(ompstatus):
        if ompstatus:
            return OMPState.get_name(ompstatus)
        else:
            return "Good"
    @app.template_filter('procstatus')
    def jsaprocstatus(status):
        if jsaprocstatus:
            return JSAProcState.get_name(status)
        else:
            return ""

    @app.template_filter('sn')
    def suppress_none(value):
        if value is not None:
            return value
        else:
            return ''
    @app.template_filter('enabled')
    def enabled(value, symbol=False):
        if bool(value) is True:
            if symbol:
                return u'✓'
            return 'enabled'
        elif bool(value) is False:
            if symbol:
                return 'x'
            return 'disabled'
        else:
            return 'unknown'
    @app.template_filter('projenc')
    def encode_project_string(projectid):
        return projectid.replace('/','-')


    @app.template_filter('formatcoords')
    def formatcoords(coord, coordtype, index):
        if coordtype=='RADEC':
            if index == 0:
                return ra_to_str(coord)
            if index == 1:
                return dec_to_str(coord)
        elif coordtype in ['PLANET','CAL','TBD']:
            return coord
        else:
            if not coord:
                return coord
            else:
                return '{:.2f}'.format(coord)

    @app.template_filter('ratostr')
    def ra_to_str(ra):
        """
        Turn floating ra (radians, degrees) to sexagesimal string.
        """
        from math import modf
        # Convert to decimal hours.
        ra = ra * 57.295779513/15
        if ra==0:
            sign=1.0
        else:
            sign = abs(ra)/ra
        frachours, hours = modf(abs(ra))
        fracminutes, minutes = modf(frachours * 60)
        seconds = 60*fracminutes
        rastring = '{:02}:{:02}:{:04.1F}'.format(int(sign*hours), int(minutes), seconds)
        return rastring

    @app.template_filter('dectostr')
    def dec_to_str(dec):
        """
        Turn floating dec (radians, degrees) to sexagesimal string.
        """
        from math import modf
        # Convert to decimal hours.
        dec = dec * 57.295779513
        if dec==0:
            sign=1.0
        else:
            sign = abs(dec)/dec

        fracdegs, degs = modf(abs(dec))
        fracminutes, minutes = modf(fracdegs * 60)
        seconds = 60*fracminutes
        decstring = '{:+02}:{:02}:{:04.1F}'.format(int(sign*degs), int(minutes), seconds)
        return decstring
    @app.template_filter('commentprint')
    def print_comment(comment):
        outputstring = ""
        if comment.author:
            outputstring = comment.author
        if comment.datetime:
            outputstring += ': {}'.format(comment.datetime.strftime('%Y-%m-%dT%H:%M:%S'))
        if comment.text and not comment.text.isspace():
            outputstring+= '\n\n{}'.format(comment.text)
        return outputstring

    @app.template_test('notwhitespace')
    def not_whitespace(text):
        if text and not text.isspace():
            return True
        else:
            return False

    @app.template_filter('ompinst')
    def omp_instrument(instname):
        if 'SCUBA-2' in instname:
            return 'SCUBA-2'
        elif '-ACSIS' in instname:
            return instname.split('-ACSIS')[0]

        else:
            return 'UNKNOWN'

    @app.context_processor
    def add_to_context():
        return{
            'url_for_omp_comment': url_for_omp_comment,
            }
    return app



def url_for_omp_comment(obsid, instrument, obsnum, date_obs):
    ompurl = 'http://omp.eao.hawaii.edu/cgi-bin'
    url = '{}/staffobscomment.pl?oid={}&inst={}&runnr={}&ut={}'.format(
        ompurl, obsid, instrument, obsnum, date_obs)
    return url
