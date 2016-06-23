#!/usr/bin/env python
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

"""
Usage: lp-web.py [-h]
       lp-web.py [-d] [--port=<port>]

Starts a web page for monitoring the progress of the JCMT large programs.

Options:
  -h --help  Show this help.
  -d --debug    Run using flasks debug mode. Allows reloading of code.
  -p, --port=<port>  Select port to run on (defaults to 5000)


"""

from __future__ import absolute_import, division, print_function

from docopt import docopt


import web


# Parse arguments
arguments = docopt(__doc__, help=True)
if arguments['--debug'] is True:
    host = '127.0.0.1'
    debug = True
else:
    host='0.0.0.0'
    debug = None

if '--port' in arguments:
    port = int(arguments['--port'])

else:
    port=None

# Start flask app.
app = web.web_pages()
app.run(host=host, debug=debug, port=port)
