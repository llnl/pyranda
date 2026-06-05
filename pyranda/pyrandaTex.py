################################################################################
# Copyright (c) 2018, Lawrence Livemore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
#
# LLNL-CODE-749864
# This file is part of pyranda
# For details about use and distribution, please read: pyranda/LICENSE
#
# Written by: Britton J. Olson, olson45@llnl.gov
################################################################################
import os
from .pyrandaUtils import runCMD, version, icopyright, fortran3d, Timeout
import random
import getpass

DIR = os.path.dirname(__file__)
BASE_TEMPLATE = os.path.join(DIR, "latex_template.tex")


class TexFigure:
    def __init__(self, caption, fig, wkdir, size):
        self.caption = caption
        self.fig = fig
        self.wkdir = wkdir
        self.fname = None
        self.size = str(size)

    def render(self):

        import matplotlib.pyplot as plt

        # Make the imagefile
        cwd = os.getcwd()
        os.chdir(self.wkdir)

        fileName = str(int(random.random() * 1.0e6))
        fileName += ".pdf"
        plt.figure(self.fig.number)
        plt.savefig(fileName, format="pdf")
        self.fname = fileName

        os.chdir(cwd)

        figure_tex = (
            r"""
        \begin{figure}[!h]
          \centering
          \includegraphics[width=$$SIZE$$\textwidth]{$$FILE$$}
          \caption{$$CAPTION$$}
        \end{figure}
        """.replace("$$CAPTION$$", self.caption)
            .replace("$$FILE$$", fileName)
            .replace("$$SIZE$$", self.size)
        )

        return figure_tex


class TexSection:
    def __init__(self, title, wkdir, level=0):
        self.title = title
        self.body = ""
        self.level = level
        self.figures = []
        self.wkdir = wkdir

    def render(self):

        sub = "sub" * self.level
        title = r"\%ssection{%s}" % (sub, self.title)

        body = self.body

        # Add figure to body
        for fig in self.figures:
            body += fig.render()

        return title + body

    def addFigure(self, caption, size=0.9):

        import matplotlib.pyplot as plt

        fig = plt.gcf()
        self.figures.append(TexFigure(caption, fig, self.wkdir, size))


class pyrandaTex:
    def __init__(self, pySim):

        self.pySim = pySim

        self.tMap = self.pySim.get_tMap()
        self.template = BASE_TEMPLATE

        self.user = getpass.getuser()

        self.sections = []
        self.packages = ""
        self.pdfFile = None
        self.wkdir = self.pySim.name + "_latex"
        if not os.path.exists(self.wkdir):
            os.makedirs(self.wkdir)

    def renderEqu(self, kind="PDE"):

        equations = []

        if (kind == "PDE") or (kind == "PDE"):
            for eq in self.pySim.equations:
                if eq.kind == kind:
                    equations.append(eq.eqstr)
        elif kind == "IC":
            equations = self.pySim.initial_conditions
        elif kind == "ALG":
            for eq in self.pySim.equations:
                if eq.kind == kind:
                    equations.append(eq.eqstr)

        # Loop over the equations of motion
        #   and return its latex string
        myEOM = ""
        for eq in equations:
            myEOM += fortran3d(eq, self.tMap, latex=True)
            myEOM += r"\\"

        section_eom = r"""
        \begin{align}
           %s
        \end{align}
        
        """ % (myEOM)
        return section_eom

    def addSection(self, title, level=0):

        isection = TexSection(title, self.wkdir, level=level)
        self.sections.append(isection)
        return isection

    def makeTex(self):

        # Main driver to produce latex doc

        # 1 - Load template
        lid = open(self.template)
        tex = lid.read()

        # Fill in form
        form = {}
        form["$$AUTHOR$$"] = self.user
        form["$$PACKAGES$$"] = self.packages
        form["$$TITLE$$"] = "Pyranda - %s" % self.pySim.name
        form["$$PYRANDA$$"] = version()
        form["$$COPY$$"] = icopyright()

        # Section
        AllSections = ""
        for section in self.sections:
            AllSections += section.render()
        form["$$SECTIONS$$"] = AllSections

        # Apply template
        for key in form:
            tex = tex.replace(key, form[key])

        return tex

    def makeDoc(self):

        tex = self.makeTex()

        # 1 - tmp dir
        directory = self.wkdir

        # 2 - Write the tex file
        texname = "pyranda.tex"
        texFile = os.path.join(directory, texname)
        pid = open(texFile, "w")
        pid.writelines(tex)
        pid.close()

        # 3 - Call tex render-er
        cwd = os.getcwd()
        os.chdir(directory)
        texExe = "pdflatex"
        try:
            with Timeout(10):
                _ = runCMD([texExe, texname])
                os.chdir(cwd)
                self.pdfFile = texFile.replace(".tex", ".pdf")

        except Timeout.Timeout:
            os.chdir(cwd)
            print("PDF render timed out")

    def showPDF(self):

        if self.pdfFile:
            try:
                _ = runCMD(["xpdf", self.pdfFile])
                return
            except Exception:
                print("Error: 'xpdf' not found")

            try:
                _ = runCMD(["open", self.pdfFile])
                return
            except Exception:
                print("Error: 'open' not found")
