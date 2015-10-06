#!/usr/bin/env python
# -*- coding: utf-8 -*- #
from __future__ import unicode_literals
import os

PLUGINS = ['tag_cloud.tag_cloud',]

AUTHOR = u'Aaron Kitzmiller'
SITENAME = u'Harvard FAS Informatics'
SITEURL = ''
TAGS_URL = 'tags'
PATH = 'content'
BANNER = True
BANNER_ALL_PAGES = False
DISPLAY_PAGES_ON_MENU = True

TIMEZONE = 'America/New_York'

DEFAULT_LANG = u'en'

# Feed generation is usually not desired when developing
FEED_ALL_ATOM = None
CATEGORY_FEED_ATOM = None
TRANSLATION_FEED_ATOM = None
AUTHOR_FEED_ATOM = None
AUTHOR_FEED_RSS = None

# Sidebar stuff
LINKS = (('Harvard University', 'http://www.harvard.edu'),
         ('Harvard FAS', 'http://www.fas.harvard.edu'),
         ('FAS Research Computing', 'http://rc.fas.harvard.edu'),)

# Social widget
# SOCIAL = (('You can add links in your config file', '#'),
#           ('Another social link', '#'),)

DISPLAY_RECENT_POSTS_ON_SIDEBAR = True
DISPLAY_TAGS_ON_SIDEBAR = True
DISPLAY_TAGS_INLINE = True
TAG_CLOUD_STEPS = 4
TAG_CLOUD_MAX_ITEMS = 100

DEFAULT_PAGINATION = 10

# Uncomment following line if you want document-relative URLs when developing
RELATIVE_URLS = True

STATIC_PATHS = ['images']
