# coding: utf8

import os, fnmatch
from datetime import datetime
from pytz import timezone
import re
from lxml import etree
import argparse

opt_defaults = {'include': ".*\.html|.*\.pdf",
                'exclude': "google.*\.html",
                'timezone': 'Europe/Paris'}


parser = argparse.ArgumentParser(description="Create an XML sitemap file from a source html folder.")
parser.add_argument("--url", help="The website address.")
parser.add_argument("--path", help="The root folder to explore.")
parser.add_argument("--include", help="The include regexp pattern.", default=opt_defaults['include'])
parser.add_argument("--exclude", help="The exclude regexp pattern.", default=opt_defaults['exclude'])
parser.add_argument("--timezone", help="The time zone.", default=opt_defaults['timezone'])
parser.add_argument("--priority", help="The page priority.", type=float, default=.8)
parser.add_argument("--changefreq", help="How often the page is likely to change.", default="monthly",
                    choices=["always", "hourly", "daily", "weekly", "monthly", "yearly", "never"])
parser.add_argument("output", help="The output file name.")

opts = parser.parse_args()

defaults = {'loc': '', 'lastmod': '', 'priority': "{:.1f}".format(opts.priority), 'changefreq': opts.changefreq}
results = []

for root, dirs, files in os.walk(opts.path):
    for fn in files:
        if re.match(opts.include, fn) and not re.match(opts.exclude, fn):
            fullfn = os.path.join(root, fn)
            loc = os.path.join(opts.url, fn)
            dt = datetime.fromtimestamp(os.stat(fullfn).st_ctime, tz=timezone(opts.timezone))
            tzd = dt.strftime("%z")
            tzd = tzd[:-2] + ":" + tzd[-2:]
            lastmod = dt.strftime("%Y-%m-%dT%H:%M") + tzd

            data = defaults.copy()
            data.update(loc=loc, lastmod=lastmod)

            results.append(data)


# convert to xml
class NameSpace:
    xsi = "http://www.w3.org/2001/XMLSchema-instance"


root = etree.Element("urlset", xmlns="http://www.sitemaps.org/schemas/sitemap/0.9", nsmap={'xsi': NameSpace.xsi})
root.attrib[etree.QName(NameSpace.xsi, 'schemaLocation')]= ("http://www.sitemaps.org/schemas/sitemap/0.9 "
                                                            "http://www.sitemaps.org/schemas/sitemap/0.9/sitemap.xsd")


comment = etree.Comment("Auto-generated on {}".format(datetime.now().ctime()))
root.insert(1, comment)

for data in results:
    url_el = etree.SubElement(root, "url")
    loc_el = etree.SubElement(url_el, "loc")
    lastmod_el = etree.SubElement(url_el, "lastmod")
    changefreq_el = etree.SubElement(url_el, "changefreq")
    priority_el = etree.SubElement(url_el, "priority")
    loc_el.text = data['loc']
    lastmod_el.text = data['lastmod']
    priority_el.text = data['priority']
    changefreq_el.text = data['changefreq']

with open(opts.output, 'w') as fd:
    fd.write(etree.tostring(root, encoding='UTF-8', xml_declaration=True, pretty_print=True).decode('utf8'))
