#!/usr/bin/env python
# -*- coding: utf-8 -*- #
from datetime import datetime
from dateutil.tz import tzlocal

MARKDOWN = {
    'extensions' : ['markdown.extensions.codehilite','markdown.extensions.fenced_code','markdown.extensions.meta','markdown.extensions.footnotes'],
    'extension_configs' : {
        'markdown.extensions.codehilite' : {
            'use_pygments' : True,
            'css_class': 'highlight',
            'linenums' : False,
            'guess_lang' : False,
        },
        "markdown.extensions.toc": {"title": "Table of Contents"},
    }
}

THEME = 'informatics-theme'

AUTHOR = u'Aaron Kitzmiller'
SITENAME = u'Harvard FAS Informatics'
SITEURL = 'https://informatics.fas.harvard.edu'
TAGS_URL = 'tags.html'
PATH = 'content'
OUTPUT_RETENTION = [".gitignore"]
BANNER = True
BANNER_ALL_PAGES = False
DISPLAY_PAGES_ON_MENU = True

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
MENUITEMS = (
    ('About', '/pages/about.html'),
)
SHOW_ARTICLE_AUTHOR = True
SHOW_DATE_MODIFIED = True
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

ARTICLE_EXCLUDES = ['workshops']
PAGE_EXCLUDES = ['workshops']
#STATIC_PATHS = ['images']
EXTRA_PATH_METADATA = {
    'extra/favicon.ico': {'path': 'favicon.ico'},
    'workshops': {'path': 'workshops'},
}

DIRECT_TEMPLATES = ['search','index','archives','tags']

CUTOFF_DATE = datetime(2090, 1, 1, 0, 0, 0, tzinfo=tzlocal())

IGNORE_FILES = ['includes', 'venv']
JINJA2CONTENT_TEMPLATES = ['includes']
