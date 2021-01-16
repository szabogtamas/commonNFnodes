import os, jinja2

from .nodePrototypes import nextflowProcess
from .nodePrototypes import introSpect

chs = introSpect.flowNodes.createChannelSpecification
intro_folder = os.path.dirname(os.path.realpath(introSpect.__file__))


class notebookEdger(nextflowProcess):
    "Render Rmd notebooks to facilitate editing of specific figures."

    def dependencies(self):
        return {
            "imports": ["import os", "import jinja2"],
            "helpers": [],
            "inhouse_packages": [self.intro_folder],
        }

    def directives(self):
        return {"publishDir": "'../notebooks', mode: 'copy'"}

    def channel_specifications(self):
        return {
            "edger_notebook_pars": chs(
                "tuple",
                (
                    "edger_nb_template",
                    "edger_script_path",
                    "file(edger_counts)",
                    "edger_condition_order",
                ),
                (
                    "edger_nb_template",
                    "edger_script_path",
                    "edger_counts",
                    "edger_condition_order",
                ),
            ),
            "edger_notebook": chs("file", "'edger_DE.Rmd'", "outFile"),
        }

    def process(
        self,
        edger_nb_template,
        edger_script_path,
        edger_counts,
        edger_condition_order,
        edger_condition_colors,
        nb_title=None,
        nb_project_dir=None,
    ):
        if nb_title in [None, "None", ""]:
            nb_title = "Calling DE genes from RNASeq with EdgeR"
        if nb_project_dir in [None, "None", ""]:
            nb_project_dir = "dirname(getwd())"
        nb_kws = dict(
            nb_title=nb_title,
            nb_project_dir=nb_project_dir,
            edger_script_path="pipeline/bin/" + os.path.basename(edger_script_path),
            edger_counts="tables/" + os.path.basename(edger_counts),
            edger_condition_order=edger_condition_order,
            edger_condition_colors=edger_condition_colors,
        )
        with open(edger_nb_template) as f:
            nb = jinja2.Template(f.read())
        return nb.render(nb_kws)

    def customize_features(self):
        node_params = dict(
            intro_folder=intro_folder,
        )
        node_params.update(self.node_params)
        self.intro_folder = node_params["intro_folder"]
        return
