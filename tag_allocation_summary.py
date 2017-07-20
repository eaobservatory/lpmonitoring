import json
import numpy as np
from summary import get_omp_database, make_completion_summary_piechart
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from omp.obs.state import OMPState
from astropy.table import Table, Column
import matplotlib
import StringIO
from omp.obs.state import OMPState

matplotlib.rcParams['savefig.transparent'] = True


affiliations = ['Canada (listed institution)',
               'Canada (national)',
               'China',
               'EAO Staff',
               'Japan',
               'South Korea',
               'Taiwan',
               'United Kingdom (listed institution)']

aff_assignment = affiliations[:]
aff_assignment.remove('EAO Staff')

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


def create_base64_figure(fig):
    fig.patch.set_visible(False)
    canvas = FigureCanvas(fig)
    img = StringIO.StringIO()
    canvas.print_png(img)
    img.seek(0)
    return img.getvalue().encode('base64')

def prepare_tagallocation_summary(semester):
    """
    semester must be 16A, 16B, 17A or 17B.
    """
    table, figsummary, bandcharts, instcharts, regioncharts = create_tagallocation_summary_objects(semester,
                                                                     filter_allocations_hours=0.1)
    webinfo = {}
    webinfo['semester'] = semester

    tab = StringIO.StringIO()
    table.write(tab, format='ascii.html', htmldict = dict(table_class='obstable display',
                                                          table_id='obs-table-scuba2'))
    tab.seek(0)
    table_html = tab.read()
    table_html = table_html.split('\n')
    headerrow_start = [table_html.index(i) for i in table_html if '<th>' in i][0] -1
    headerrow_end = [table_html.index(i) for i in table_html if '</th>' in i][-1] +1
    header = table_html[headerrow_start: headerrow_end+1]
    header2 = header[:]
    header2[0] = header2[0].replace('<tr>', '<tr id="filterrow">')
    [table_html.insert(headerrow_end+1+i, header2[i]) for i in range(len(header2))]
    table_html = '\n'.join(table_html)
    webinfo['table'] = table_html

    webinfo['figsummary'] = create_base64_figure(figsummary)
    for b in bandcharts:
        bandcharts[b] = create_base64_figure(bandcharts[b])
    for r in regioncharts:
        regioncharts[r] = create_base64_figure(regioncharts[r])
    for i in instcharts:
        instcharts[i] = create_base64_figure(instcharts[i])
    webinfo['regioncharts'] = regioncharts
    webinfo['bandcharts'] = bandcharts
    webinfo['instcharts'] = instcharts

    return webinfo



def create_tagallocation_summary_objects(semester, filter_allocations_hours=0.0):
    """
    Create summary's of semester's based on OMP and allocation information.

    Only supports PI queue.

    Filter out allocations less than filter_allocations_hours from pie charts
    """


    # Read JSON information into nexted dictionaries.
    semester_file = 'M_{}_P_info.json'.format(semester)
    with open(semester_file, 'r') as file_:
        proposal_info = json.load(file_)

    # First get all the allocated projects. Ignore all unaccepted proposals.
    allocated_projects = [p for p in proposal_info if proposal_info[p]['state']=='Accepted']


    # Get allocated projects sorted by rating.
    ratings = [proposal_info[p]['rating'] for p in allocated_projects]
    projrat = zip(allocated_projects, ratings)
    projrat.sort(key=lambda x: x[1], reverse=True)
    allocated_projects_rating = np.array(projrat)[:,0]


    ompdb = get_omp_database()
    ompallocs = ompdb.get_allocations(projects=allocated_projects_rating)
    allocations_overall = [60.0*60.0*sum([a['time'] for a in proposal_info[p]['allocation']]) for p in allocated_projects_rating]
    observedtimes = np.array([ompallocs[p].allocated - ompallocs[p].remaining for p in allocated_projects_rating])
    fig_summary = Figure(figsize=[8,8])
    axpie = fig_summary.add_subplot(111)
    axpie = make_completion_summary_piechart(axpie, allocated_projects_rating,
                                             allocations_overall, observedtimes, radius=0.8)



    # Create  table containing the results.
    tableheadings = [
        'project',
        'PI',
        'PI affiliation',
        'Exempt',
        'Rating',
        'Completion %',
        'Hours allocated',
        'Hours charged',
        'Alloc: B1',
        'Alloc: B2',
        'Alloc: B3',
        'Alloc: B4',
        'Alloc: B5',
        'Alloc: HARP',
        'Alloc: RxA3',
        'Alloc: SCUBA-2',
        'Alloc: POL-2',
        'Obs: B1',
        'Obs: B2',
        'Obs: B3',
        'Obs: B4',
        'Obs: B5',
        'Obs: HARP',
        'Obs: RxA3',
        'Obs: SCUBA-2',
        'Obs: POL-2',
        'Tech asses',

        'Affil: Can. List',
        'Affil: Can. Nat',
        'Affil: China',
        'Affil: Japan',
        'Affil: South Korea',
        'Affil: Taiwan',
        'Affil: UK',
        'ToP',
        ]

    def remove_0(value):
        if value == 0:
            return ''
        else:
            return '{:.1f}'.format(value)

    def bool_to_check(value):
        if value == False:
            return u'\u2717'
        else:
            return u'\u2713'

    formatdict = {}
    othervalues = ['project', 'PI', 'PI affiliation', 'Exempt', 'Tech asses', 'ToP']
    for h in tableheadings:
        if h not in othervalues:
            formatdict[h] = remove_0

    formatdict['Exempt'] = bool_to_check
    formatdict['ToP'] = bool_to_check


    table = []
    for project in allocated_projects_rating:
        print(project)
        info = proposal_info[project]
        allocations = info['allocation']
        hours_allocated = sum([a['time'] for a in allocations])
        hours_charged = (ompallocs[project].allocated - ompallocs[project].remaining)/(60.0*60.0)
        percent_complete = 100 * hours_charged / hours_allocated


        allocated_bands = set([a['weather'] for a in allocations])
        allocated_bands = [int(i.split()[1]) for i in allocated_bands]


        obs = ompdb.get_observations_from_project(project)
        instruments_observed, bands_observed, extra_time = get_bands_instruments_observed(allocated_bands,
                                                                                          obs)
        affiliation_values = [100 * info['affiliation_assignment'][e] if  e in info['affiliation_assignment'] else 0.0 for e in aff_assignment ]

        tableline = [
            project,
            info['member_pi']['person_name'],
            info['member_pi']['affiliation_name'],
            info['decision_exempt'],
            info['rating'],
            percent_complete,
            hours_allocated,
            hours_charged,
            sum([a['time'] for a in allocations if a['weather']=='Band 1']),
            sum([a['time'] for a in allocations if a['weather']=='Band 2']),
            sum([a['time'] for a in allocations if a['weather']=='Band 3']),
            sum([a['time'] for a in allocations if a['weather']=='Band 4']),
            sum([a['time'] for a in allocations if a['weather']=='Band 5']),
            sum([a['time'] for a in allocations if a['instrument']=='HARP']),
            sum([a['time'] for a in allocations if 'RXA3' in a['instrument'].upper()]),
            sum([a['time'] for a in allocations
                 if a['instrument']=='SCUBA-2'
                 and a['ancillary'] is None
                 and info['jcmt_options']['polarimetry'] is False]),
            sum([a['time'] for a in allocations
                 if a['instrument']=='SCUBA-2'
                 and (a['ancillary']=='POL-2'
                      or info['jcmt_options']['polarimetry'] is True)]),
            bands_observed.get(1, 0.0)/(60.0*60.0),
            bands_observed.get(2, 0.0)/(60.0*60.0),
            bands_observed.get(3, 0.0)/(60.0*60.0),
            bands_observed.get(4, 0.0)/(60.0*60.0),
            bands_observed.get(5, 0.0)/(60.0*60.0),
            instruments_observed.get('HARP', 0.0)/(60.0*60.0),
            instruments_observed.get('RxA3', 0.0)/(60.0*60.0),
            instruments_observed.get('SCUBA-2', 0.0)/(60.0*60.0),
            instruments_observed.get('POL-2', 0.0)/(60.0*60.0),
            info['review_technical']['review_assessment'],
            ]
        tableline += affiliation_values
        tableline += [info['jcmt_options']['target_of_opp']]
        table.append(tableline)

    table = Table(data=map(list, zip(*table)), names=tableheadings)
    print(table)
    for key in formatdict:
        table.replace_column(key, Column(name=key, data=table[key].data, format=formatdict[key]))

    # Now make piecharts from Table.


    # Get total time.
    totaltime = np.sum(table['Hours charged'])

    regioncharts = {}
    defaulttime = totaltime/(len(affiliations))
    #1. by PI affilation
    for r in affiliations:
        fig = Figure(figsize=(5,5))
        ax = fig.add_subplot(111)
        affmask = table['PI affiliation']==r
        projects = table[affmask]['project']
        allocations = table[affmask]['Hours allocated']
        observed = table[affmask]['Hours charged']
        fractional_time = sum(allocations)/defaulttime
        radius = 0.6 * np.sqrt(fractional_time)

        make_completion_summary_piechart(ax, projects, allocations, observed,
                                         scale=1.5, radius=radius)
        regioncharts[r] = fig


    # Completion by instrument.
    instcharts = {}
    instruments = ['HARP', 'RxA3', 'SCUBA-2', 'POL-2']
    defaulttime = totaltime/len(instruments)
    for i in instruments:
        fig = Figure(figsize=(5,5))
        ax = fig.add_subplot(111)
        imask = (table['Alloc: {}'.format(i)] != 0) & (table['Alloc: {}'.format(i)] > filter_allocations_hours)
        projects = table[imask]['project']
        allocations = table[imask]['Alloc: {}'.format(i)]
        observed = table[imask]['Obs: {}'.format(i)]
        fractional_time = sum(allocations)/defaulttime
        radius = 0.6 * np.sqrt(fractional_time)
        make_completion_summary_piechart(ax, projects, allocations, observed,
                                         scale=1.5, radius=radius)
        instcharts[i] = fig

    # Completion by Band.
    bandcharts = {}
    bands = list(range(1,6))
    defaulttime = totaltime/len(bands)
    for b in bands:
        fig = Figure(figsize=(5,5))
        ax = fig.add_subplot(111)
        bmask = (table['Alloc: B{}'.format(b)] != 0) & (table['Alloc: B{}'.format(b)] > filter_allocations_hours)
        projects = table[bmask]['project']
        allocations = table[bmask]['Alloc: B{}'.format(b)]
        observed = table[bmask]['Obs: B{}'.format(b)]
        fractional_time = sum(allocations)/defaulttime
        radius = 0.6 * np.sqrt(fractional_time)
        make_completion_summary_piechart(ax, projects, allocations, observed,
                                         scale=1.5, radius=radius)
        bandcharts[b] = fig

    return table, fig_summary, bandcharts, instcharts, regioncharts



def get_bands_instruments_observed(allocated_bands, obs):
    extra_time = {}
    bands_observed = {}
    instruments_observed = {}

    # Divvy up the observed time by a. requested BAND if available
    # from req_mintau, req_maxtau header.  If observed bektween
    # requested bands, get the actual band observed in. If
    # observed outside requested bands, use closest band.


    # Get the observations
    observations = obs.values()

    # Skip out all JUNK, REJECTED and BAD observation:
    observations = [o for o in observations if o.status in [OMPState.GOOD, OMPState.QUESTIONABLE, None]]

    for o in observations:

        instrument = o.instrument
        if 'RXA3' in instrument.upper():
            instrument = 'RxA3'
        duration = o.duration
        instruments_observed[instrument] = instruments_observed.get(instrument, 0.0) + duration
        tau = np.mean([o.wvmtaust, o.wvmtauen])
        tau = float('{:.2F}'.format(tau))
        observed_band = Bands.get_band(tau)
        requested_lower_band = Bands.get_band(o.req_mintau)
        requested_upper_band = Bands.get_band(o.req_maxtau)
        if requested_lower_band <= observed_band <= requested_upper_band:
            band = observed_band
        elif observed_band < requested_lower_band:
            band = requested_lower_band
        elif observed_band > requested_upper_band:
            band = requested_upper_band
        else:
            band = None
        if band not in allocated_bands:
            # Check if its right on boundary:
            orig_band = band
            extra = None
            if band > max(allocated_bands):
                if tau <= 1.05* Bands.get_fromname(max(allocated_bands))[1]:
                    extra = None
                    band = max(allocated_bands)
                else:
                    extra = True
                    band = max(allocated_bands)
            elif band < min(allocated_bands):
                if tau >= 0.95 * Bands.get_fromname(min(allocated_bands))[0]:
                    extra = None
                    band = min(allocated_bands)
                else:
                    extra = True
                    band = min(allocated_bands)

            if extra:
                print('Warning: observed {}s in Band {} ({}), but was allocated band/s {}'.format(
                    duration, orig_band, tau, allocated_bands))
                print('Observation had requested tau of {}:{}'.format(o.req_mintau, o.req_maxtau))
                print('Observation charged to band {}'.format(band))
                extra_time[band] = extra_time.get(band, 0.0) + duration

        bands_observed[band] = bands_observed.get(band, 0.0) + duration
    return instruments_observed, bands_observed, extra_time
