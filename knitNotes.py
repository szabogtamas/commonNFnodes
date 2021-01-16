from .nodePrototypes import nextflowCmdProcess


class knitRmd(nextflowCmdProcess):
    "Knit Rmd notebooks to both pdf and html."

    def customize_features(self):
        self.inputs = [
            "file notebook from notebooks",
        ]
        self.outputs = [
            "file '*.pdf' into notebook_pdfs",
            "file '*.html' into notebook_htmls",
        ]
        self.command = """Rscript -e \"library(dplyr); rmarkdown::render('${notebook}', 'html_document'); rmarkdown::render('${notebook}', 'pdf_document')\""""
        return None
