import os, jinja2

from .nodePrototypes import nextflowProcess
from .nodePrototypes import introSpect

chs = introSpect.flowNodes.createChannelSpecification
intro_folder = os.path.dirname(os.path.realpath(introSpect.__file__))


class notebookTemplater(nextflowProcess):
    "Render Rmd notebooks to facilitate editing of specific figures."

    def dependencies(self):
        return {
            "imports": ["import os", "import jinja2"],
            "helpers": [],
            "inhouse_packages": [self.intro_folder],
        }

    def directives(self):
        return {"publishDir": "'../notebooks', mode: 'copy'"}

    def customize_features(self):
        node_params = dict(
            intro_folder=intro_folder,
        )
        node_params.update(self.node_params)
        self.intro_folder = node_params["intro_folder"]
        return


class notebookEdger(notebookTemplater):
    "Render Rmd notebooks to facilitate editing of specific figures."

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


class notebookProgeny(notebookTemplater):
    "Render Rmd notebooks to facilitate editing of specific figures."

    def channel_specifications(self):
        return {
            "progeny_notebook_pars": chs(
                "tuple",
                (
                    "progeny_nb_template",
                    "progeny_script_path",
                    "file(progeny_counts)",
                    "progeny_condition_order",
                    "species",
                ),
                (
                    "progeny_nb_template",
                    "progeny_script_path",
                    "progeny_counts",
                    "progeny_condition_order",
                    "species",
                ),
            ),
            "progeny_notebook": chs("file", "'progeny_paths.Rmd'", "outFile"),
        }

    def process(
        self,
        progeny_nb_template,
        progeny_script_path,
        progeny_counts,
        progeny_condition_order,
        progeny_condition_colors,
        species,
        nb_title=None,
        nb_project_dir=None,
    ):
        if nb_title in [None, "None", ""]:
            nb_title = "Pathway perturbation scores from PROGENy"
        if nb_project_dir in [None, "None", ""]:
            nb_project_dir = "getwd()"
        nb_kws = dict(
            nb_title=nb_title,
            nb_project_dir=nb_project_dir,
            progeny_script_path="pipeline/bin/" + os.path.basename(progeny_script_path),
            progeny_counts="tables/" + os.path.basename(progeny_counts),
            progeny_condition_order=progeny_condition_order,
            progeny_condition_colors=progeny_condition_colors,
            progeny_species=species,
        )
        with open(progeny_nb_template) as f:
            nb = jinja2.Template(f.read())
        return nb.render(nb_kws)


class notebookGO(notebookTemplater):
    "Render Rmd notebooks to facilitate editing of specific figures."

    def channel_specifications(self):
        return {
            "go_notebook_pars": chs(
                "tuple",
                (
                    "go_nb_template",
                    "go_script_path",
                    "go_hitlists",
                    "order_of_conditions",
                    "go_tag",
                    "go_prefix",
                    "species",
                    "mcat",
                    "msubcat",
                ),
                (
                    "go_nb_template",
                    "go_script_path",
                    "go_hitlists",
                    "go_condition_order",
                    None,
                    None,
                    "go_species",
                    "go_maincat",
                    "go_subcat",
                ),
            ),
            "hitlist_notebook": chs(
                "file",
                '"${go_prefix}${go_tag.replace(' + "'\\\"', ''" + ')}.Rmd"',
                "outFile",
            ),
        }

    def process(
        self,
        go_nb_template,
        go_script_path,
        go_hitlists,
        go_condition_order,
        nb_title=None,
        nb_project_dir=None,
        go_species=None,
        go_maincat=None,
        go_subcat=None,
    ):
        if nb_title in [None, "None", ""]:
            nb_title = "GO terms associated to DE genes RNASeq data"
        if nb_project_dir in [None, "None", ""]:
            nb_project_dir = "dirname(getwd())"
        if go_species in [None, "None", ""]:
            go_species = "Mus musculus"
        if go_maincat in [None, "None", ""]:
            go_maincat = "C5"
        if go_subcat in [None, "None", ""]:
            go_subcat = "BP"
        nb_kws = dict(
            nb_title=nb_title,
            nb_project_dir=nb_project_dir,
            go_script_path="pipeline/bin/" + os.path.basename(go_script_path),
            go_hitlists=go_hitlists,
            go_condition_order=go_condition_order,
            go_species=go_species,
            go_maincat=go_maincat,
            go_subcat=go_subcat,
        )
        with open(go_nb_template) as f:
            nb = jinja2.Template(f.read())
        return nb.render(nb_kws)


class notebookGSEA(notebookTemplater):
    "Render Rmd notebooks to facilitate editing of specific figures."

    def process(
        self,
        gsea_nb_template,
        gsea_script_path,
        gsea_input_tables,
        gsea_condition_order,
        nb_title=None,
        nb_project_dir=None,
        gsea_species=None,
        gsea_maincat=None,
        gsea_subcat=None,
    ):
        if nb_title in [None, "None", ""]:
            nb_title = "GSEA of RNASeq data"
        if nb_project_dir in [None, "None", ""]:
            nb_project_dir = "getwd()"
        if gsea_species in [None, "None", ""]:
            gsea_species = "Mus musculus"
        if gsea_maincat in [None, "None", ""]:
            gsea_maincat = "H"
        if gsea_subcat in [None, "None", ""]:
            gsea_subcat = "NULL"

        if isinstance(gsea_input_tables, str):
            gsea_input_tables = gsea_input_tables.split(",")

        nb_kws = dict(
            nb_title=nb_title,
            nb_project_dir=nb_project_dir,
            gsea_script_path="pipeline/bin/" + os.path.basename(gsea_script_path),
            gsea_input_tables=",".join(
                [
                    x.split(":")[0]
                    + ":"
                    + "tables/"
                    + os.path.basename(x.split(":")[1])
                    for x in gsea_input_tables
                ]
            ),
            gsea_condition_order=gsea_condition_order,
            gsea_species=gsea_species,
            gsea_maincat=gsea_maincat,
            gsea_subcat=gsea_subcat,
        )

        with open(gsea_nb_template) as f:
            nb = jinja2.Template(f.read())
        return nb.render(nb_kws)


class notebookPositional(notebookTemplater):
    "Render Rmd notebooks to facilitate editing of specific figures."

    def channel_specifications(self):
        return {
            "positional_notebook_pars": chs(
                "tuple",
                (
                    "positional_nb_template",
                    "positional_script_path",
                    "positional_tabs",
                    "order_of_conditions",
                    "positional_prefix",
                    "species",
                ),
                (
                    "positional_nb_template",
                    "positional_script_path",
                    "positional_input_tables",
                    "positional_condition_order",
                    None,
                    "gsea_species",
                ),
            ),
            "positional_notebook": chs("file", '"${positional_prefix}.Rmd"', "outFile"),
        }

    def process(
        self,
        positional_nb_template,
        positional_script_path,
        positional_input_tables,
        positional_condition_order,
        nb_title=None,
        nb_project_dir=None,
        gsea_species=None,
    ):
        if nb_title in [None, "None", ""]:
            nb_title = "Overrrepresented chromosomal positions"
        if nb_project_dir in [None, "None", ""]:
            nb_project_dir = "getwd()"
        if gsea_species in [None, "None", ""]:
            gsea_species = "Mus musculus"

        if isinstance(positional_input_tables, str):
            positional_input_tables = positional_input_tables.split(",")

        nb_kws = dict(
            nb_title=nb_title,
            nb_project_dir=nb_project_dir,
            positional_script_path="pipeline/bin/"
            + os.path.basename(positional_script_path),
            positional_input_tables=",".join(
                [
                    x.split(":")[0]
                    + ":"
                    + "tables/"
                    + os.path.basename(x.split(":")[1])
                    for x in positional_input_tables
                ]
            ),
            positional_condition_order=positional_condition_order,
            gsea_species=gsea_species,
        )

        with open(positional_nb_template) as f:
            nb = jinja2.Template(f.read())
        return nb.render(nb_kws)
