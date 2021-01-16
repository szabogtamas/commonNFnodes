import os
import pandas as pd

from .nodePrototypes import nextflowCmdProcess


class pdfFromLatex(nextflowCmdProcess):
    """
    Convert the LaTeX format report into pdf.
    """

    def channel_pretreat(self):
        return [
            [
                "Channel",
                "fromPath(params.reference_file)",
                "set{reference_file}",
            ],
        ]

    def directives(self):
        return {
            "publishDir": "'../', mode: 'copy'",
        }

    def customize_features(self):
        self.inputs = [
            "file reportex",
            "file 'references.bib' from reference_file",
        ]
        self.outputs = ['file("*.pdf")']
        self.command = "pdflatex -shell-escape -interaction nonstopmode -halt-on-error -file-line-error $reportex\n            "
        self.command += "biber ${reportex.baseName}\n            "
        self.command += "pdflatex -shell-escape -interaction nonstopmode -halt-on-error -file-line-error $reportex\n"
        return None
