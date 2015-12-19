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

import functools
import os
from flask import Flask, redirect, url_for, render_template

from summary import create_summary, get_omp_database, prepare_completion_chart, prepare_project_chart



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


    # Setup root page to redirect to summary
    @app.route('/')
    def home_page():
        return redirect(url_for('summary'), code=303)

    # Summary page
    @app.route('/summary/')
    def summary():
        values = create_summary(ompdb)
        return render_template('summary.html', **values)

    # Chart showing overall completion by project.
    @app.route('/summary/completionchart/')
    def completion_chart():
        return prepare_completion_chart(ompdb)

    @app.route('/summary/chart/<projectid>/')
    def project_chart(projectid):
        summary = ompdb.get_summary_obs_info(str(projectid))
        return prepare_project_chart(projectid, summary)

    # Get all the bands for which any data is allocated from an
    # allocation dictionary.
    @app.template_filter('allocationbandset')
    def allocationbandset(allocdict):
        return set([item for sublist in allocdict.values() for item in sublist.keys()])

    # A filter to get the set of an attributes values from a list
    @app.template_filter('getset')
    def getset(listofthings, attribute):
        return set([getattr(i,attribute) for i in listofthings])
    return app

