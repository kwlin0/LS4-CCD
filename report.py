#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys

nargs = len(sys.argv)

if nargs != 2:
    sys.exit("No CCD ID entered")
else:
    CCD_ID = sys.argv[1]

from pylatex import Document, PageStyle, Head, Section, Subsection, Command, Figure, Package, Math
from pylatex.utils import bold, italic, NoEscape


def fillout_report(doc, CCD_ID):
    """Add a section, a subsection and some text to the document.

    :param doc: the document
    :type doc: :class:`pylatex.document.Document` instance
    """
    doc.preamble.append(Command('title',  NoEscape(r'\color{red} LS4 CCD '+CCD_ID+' Testing Report')))
    doc.preamble.append(Command('author',  'Kenneth Lin'))
    doc.preamble.append(Command('date', NoEscape(r'\today')))
    doc.append(NoEscape(r'\maketitle'))

    with doc.create(Section('Flats')):
        doc.append('High-light PTC flats were taken up to an exposure time of 450s at 900 nm with monochromator fixed at filter 1 and grating 2.')
        with doc.create(Figure(position='h!')) as ptc0:
                ptc0.add_image('s3-246_Relative_QE_2')
                ptc0.add_caption('High-light PTC')

    with doc.create(Section('Quantum Efficiency')):
        doc.append('Flats were taken from 400 to 1100 nm in 25 nm increments with 15s exposure time. Images are dark-subtracted.')
        
        with doc.create(Figure(position='h!')) as qe0:
                qe0.add_image('s3-246_Relative_QE_2')
                qe0.add_caption('QE curve 1')
    
    with doc.create(Section('Bad pixel estimation')):
        doc.append("A series of darks are taken to estimate hot and bad pixels by fitting pixels on the active surface to a Gaussian and masking pixels \
        over ")
        doc.append(NoEscape(r"$5\sigma$ "))
        doc.append("from the mean. Images are overscan-subtracted.")


if __name__ == '__main__':
    
    geometry_options = {"tmargin": "1in", "lmargin": "1in"}
    doc = Document(geometry_options=geometry_options)
    doc.packages.append(Package('color'))
    header = PageStyle("header")
    
    with header.create(Head("R")):
        header.append(CCD_ID + ' Testing Report')
    
    doc.preamble.append(header)
    doc.change_document_style("header")
    
    fillout_report(doc, CCD_ID)

    doc.generate_pdf(CCD_ID + '_report', clean_tex=True, compiler='/global/u1/k/kwlin/texlive/2021/bin/x86_64-linux/pdflatex')
    tex = doc.dumps()
    print(CCD_ID + '_report.pdf written.')