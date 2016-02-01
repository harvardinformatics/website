#!/usr/bin/env python
# -*- coding: utf-8 -*- #
from __future__ import unicode_literals
import os
from datetime import datetime
from dateutil.tz import tzlocal

PLUGIN_PATHS = ['../pelican-plugins/']
PLUGINS = ['tag_cloud.tag_cloud',
           'interlinks',
           'pelican-toc',
           'tipue_search',
           'rmd_reader',
           'pin_to_top',  
           ]
KNITR_OPTS_CHUNK = {
    'fig.path': 'images/', 
}
MD_EXTENSIONS = (['codehilite','fenced_code'])
TOC = {
    'TOC_HEADERS' : '^h[1-3]',  # What headers should be included in the generated toc
                                # Expected format is a regular expression
    'TOC_RUN'     : 'true'      # Default value for toc generation, if it does not evaluate
                                # to 'true' no toc will be generated
}
AUTHOR = u'Aaron Kitzmiller'
SITENAME = u'Harvard FAS Informatics'
SITEURL = 'http://localhost:8000'
TAGS_URL = 'tags.html'
# PATH = 'content2'
BANNER = True
BANNER_ALL_PAGES = False
DISPLAY_PAGES_ON_MENU = True
FAVICON = 'favicon.ico'

TIMEZONE = 'America/New_York'

DEFAULT_LANG = u'en'

PYGMENTS_STYLE = 'github'

# Feed generation is usually not desired when developing
FEED_ALL_ATOM = None
FEED_ALL_RSS = 'rss.xml'
CATEGORY_FEED_ATOM = None
TRANSLATION_FEED_ATOM = None
AUTHOR_FEED_ATOM = None
AUTHOR_FEED_RSS = None

# Sidebar stuff, Link images are first
LINKS = (('/images/harvard.jpg','Harvard University', 'http://www.harvard.edu'),
         ('/images/harvard.jpg','Harvard FAS', 'http://www.fas.harvard.edu'),
         ('/images/odybot.png','FAS Research Computing', 'http://rc.fas.harvard.edu'),)

SHOW_ARTICLE_AUTHOR = True
# Social widget
# SOCIAL = (('You can add links in your config file', '#'),
#           ('Another social link', '#'),)

DISPLAY_RECENT_POSTS_ON_SIDEBAR = False
DISPLAY_TAGS_ON_SIDEBAR = True
DISPLAY_TAGS_INLINE = True
TAG_CLOUD_STEPS = 4
TAG_CLOUD_MAX_ITEMS = 100

DISPLAY_CATEGORIES_ON_SIDEBAR = False
DISPLAY_CATEGORIES_ON_MENU = True

DEFAULT_PAGINATION = 10

# Uncomment following line if you want document-relative URLs when developing
RELATIVE_URLS = True

STATIC_PATHS = ['images','extra']
EXTRA_PATH_METADATA = {
    'extra/favicon.ico': {'path': 'favicon.ico'},
}

# Common URLs
INTERLINKS = {
    'account_request': 'https://account.rc.fas.harvard.edu/request/',
    'password_reset': 'https://account.rc.fas.harvard.edu/password_reset/',
    'openauth': 'https://software.rc.fas.harvard.edu/oa/',
    'revoke': 'https://software.rc.fas.harvard.edu/oa/revoke/',
    'module_list': 'https://portal.rc.fas.harvard.edu/apps/modules',
    'slurm': 'http://slurm.schedmd.com/',
    'rc_site': 'https://rc.fas.harvard.edu/',
    'informatics_site': 'http://informatics.fas.harvard.edu/',
    'rchelp' : 'https://portal.rc.fas.harvard.edu/rcrt/submit_ticket',
    'lustre' : 'http://wiki.lustre.org/index.php/Main_Page',
    'rcvpn' : 'https://vpn.rc.fas.harvard.edu',
    'access-and-login' : 'https://rc.fas.harvard.edu/docs/access-and-login.html',
}

DIRECT_TEMPLATES = ['search','index','archives','tags']

CUTOFF_DATE = datetime(2090, 1, 1, 0, 0, 0, tzinfo=tzlocal())