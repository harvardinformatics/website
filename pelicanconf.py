#!/usr/bin/env python
# -*- coding: utf-8 -*- #
from __future__ import unicode_literals
import os

PLUGIN_PATHS = ['../pelican-plugins/']
PLUGINS = ['tag_cloud.tag_cloud','interlinks','extract_toc','tipue_search']
MD_EXTENSIONS = (['toc(permalink=true)','codehilite'])

AUTHOR = u'Aaron Kitzmiller'
SITENAME = u'Harvard FAS Informatics'
SITEURL = ''
TAGS_URL = 'tags'
PATH = 'content'
BANNER = True
BANNER_ALL_PAGES = False
DISPLAY_PAGES_ON_MENU = True
FAVICON = 'favicon.ico'

TIMEZONE = 'America/New_York'

DEFAULT_LANG = u'en'

PYGMENTS_STYLE = 'github'

# Feed generation is usually not desired when developing
FEED_ALL_ATOM = None
CATEGORY_FEED_ATOM = None
TRANSLATION_FEED_ATOM = None
AUTHOR_FEED_ATOM = None
AUTHOR_FEED_RSS = None

# Sidebar stuff, Link images are first
LINKS = (('https://lh5.googleusercontent.com/-5HgubPx3b20/AAAAAAAAAAI/AAAAAAAAE1U/HXqd-gt9HwI/s0-c-k-no-ns/photo.jpg','Harvard University', 'http://www.harvard.edu'),
         ('https://lh5.googleusercontent.com/-5HgubPx3b20/AAAAAAAAAAI/AAAAAAAAE1U/HXqd-gt9HwI/s0-c-k-no-ns/photo.jpg','Harvard FAS', 'http://www.fas.harvard.edu'),
         ('/images/odyssey_branding.png','FAS Research Computing', 'http://rc.fas.harvard.edu'),)

# Social widget
# SOCIAL = (('You can add links in your config file', '#'),
#           ('Another social link', '#'),)

DISPLAY_RECENT_POSTS_ON_SIDEBAR = False
DISPLAY_TAGS_ON_SIDEBAR = True
DISPLAY_TAGS_INLINE = True
TAG_CLOUD_STEPS = 4
TAG_CLOUD_MAX_ITEMS = 100

DISPLAY_CATEGORIES_ON_SIDEBAR = False

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

