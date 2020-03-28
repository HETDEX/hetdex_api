# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 09:38:06 2017

@author: gregz
"""

import sys
import datetime
import logging

import argparse as ap
from datetime import datetime as dt


def setup_parser():
    ''' BRIEF DESCRIPTION '''
    parser = ap.ArgumentParser(add_help=True)

    parser.add_argument("-sd", "--start_date",
                        help='''Start Date, e.g., 20170321, YYYYMMDD''',
                        type=str, default=None)

    parser.add_argument("-ed", "--end_date",
                        help='''Start Date, e.g., 20170326, YYYYMMDD''',
                        type=str, default=None)

    parser.add_argument("-dl", "--date_length",
                        help='''Days after/before start/end date, e.g., 10''',
                        type=int, default=None)

    parser.add_argument("-r", "--rootdir",
                        help='''Root Directory for Date''',
                        type=str, default='/work/03946/hetdex/maverick')

    parser.add_argument("-in", "--instrument",
                        help='''Instrument, e.g., virus''',
                        type=str, default='virus')

    return parser


def setup_basic_parser():
    ''' BRIEF DESCRIPTION '''
    parser = ap.ArgumentParser(add_help=True)

    parser.add_argument("-d", "--date",
                        help='''Date, e.g., 20170321, YYYYMMDD''',
                        type=str, default=None)

    parser.add_argument("-o", "--observation",
                        help='''Observation number, "00000007" or "7"''',
                        type=str, default=None)

    parser.add_argument("-e", "--exposure_number",
                        help='''Exposure number, 10''',
                        type=int, default=None)

    parser.add_argument("-r", "--rootdir",
                        help='''Root Directory for Reductions''',
                        type=str, default='/work/03946/hetdex/maverick')

    parser.add_argument("-in", "--instrument",
                        help='''Instrument, e.g., lrs2''',
                        type=str, default='lrs2')

    parser.add_argument("-i", "--ifuslot",
                        help='''Ifuslot, e.g., 066''',
                        type=str, default='066')

    parser.add_argument("-s", "--side",
                        help='''Instrument Side, e.g., L''',
                        type=str, default='L')

    return parser


def setup_logging(logname='input_utils'):
    '''Set up a logger for shuffle with a name ``input_utils``.

    Use a StreamHandler to write to stdout and set the level to DEBUG if
    verbose is set from the command line
    '''
    log = logging.getLogger('input_utils')
    if not len(log.handlers):
        fmt = '[%(levelname)s - %(asctime)s] %(message)s'
        fmt = logging.Formatter(fmt)

        level = logging.INFO

        handler = logging.StreamHandler()
        handler.setFormatter(fmt)
        handler.setLevel(level)

        log = logging.getLogger('input_utils')
        log.setLevel(logging.DEBUG)
        log.addHandler(handler)
    return log


def set_daterange(args):
    dateatt = ['start_date', 'end_date']
    if args.date_length is None:
        if args.start_date is None:
            args.log.error('You must include two of the following: '
                           '"start_date", "end_date", or "date_length"')
            sys.exit(1)
        if args.end_date is None:
            args.log.error('You must include two of the following: '
                           '"start_date", "end_date", or "date_length"')
            sys.exit(1)
        dates = {}
        for da in dateatt:
            dates[da] = dt(int(getattr(args, da)[:4]),
                           int(getattr(args, da)[4:6]),
                           int(getattr(args, da)[6:]))

        args.daterange = [datetime.date.fromordinal(i)
                          for i in range(dates[dateatt[0]].toordinal(),
                                         dates[dateatt[1]].toordinal())]
    else:
        if args.start_date is not None and args.end_date is not None:
            args.log.warning('Using "start_date" and "date_length", '
                             'however, you specified "end_date" as well '
                             'which will not be used.')
            args.end_date = None
        if args.start_date is not None:
            base = dt(int(args.start_date[:4]),
                      int(args.start_date[4:6]),
                      int(args.start_date[6:]))
            args.daterange = [base + datetime.timedelta(days=x)
                              for x in range(0, args.date_length)]

        if args.end_date is not None:
            base = dt(int(args.end_date[:4]),
                      int(args.end_date[4:6]),
                      int(args.end_date[6:]))
            args.daterange = [base - datetime.timedelta(days=x)
                              for x in range(0, args.date_length)]

    return args
