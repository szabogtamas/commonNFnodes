from .nodePrototypes import nextflowRscriptProcess
from .nodePrototypes import introSpect

chs = introSpect.flowNodes.createChannelSpecification
rcaller_script = "commandR.r"
gsea_script = "gsea_ontologies.r"

class pathwayActivityProgeny(nextflowRscriptProcess):
    "Run an R script wrapping PROGENy to assess pathway activity."

    def channel_pretreat(self):
        return [
            [
                "Channel",
                "fromPath(params.command_r)",
                "set{command_r}",
            ],
        ]

    def directives(self):
        return {
            "publishDir": """'.', saveAs: { it.contains('.tsv') ? "../tables/$it" : "../figures/$it" }, mode: 'copy'"""
        }

    def dependencies(self):
        return {
            "inhouse_packages": [self.rcaller_script],
        }

    def channel_specifications(self):
        return {
            "command_r": chs("file", "'commandR.r'", "--commandRpath"),
            "count_matrix2": chs("file", "count_matrix", "-i"),
            "conditions": chs("val", "conditions", "-l", derive_from_params=True),
            "order_of_conditions": chs(
                "val",
                "order_of_conditions",
                "--conditionOrder",
                derive_from_params=True,
            ),
            "progeny_prefix": chs(
                "val", "progeny_prefix", "-o", derive_from_params=True
            ),
            "progeny_title": chs("val", "progeny_title", None, derive_from_params=True),
            "species_alias": chs(
                "each", "species", "--expSpecies", derive_from_params=True
            ),
            "progeny_nb_template": chs(
                "each", "progeny_nb_template", None, derive_from_params=True
            ),
            "progeny_script_path": chs(
                "each", "progeny_script_path", None, derive_from_params=True
            ),
            "progeny_notebook_pars": chs(
                "tuple",
                (
                    "progeny_nb_template",
                    "progeny_script_path",
                    "file(count_matrix)",
                    "order_of_conditions",
                    "species",
                ),
                (None, None, None, None, None),
            ),
            "progeny_tables": chs("file", "'*.tsv'", None),
            "progeny_figures": chs("tuple", ("progeny_title", "'*.pdf'"), (None, None)),
        }

    def customize_features(self):
        node_params = dict(rcaller_script=rcaller_script)
        node_params.update(self.node_params)
        self.rcaller_script = node_params["rcaller_script"]
        return


class visualiseOntologiesForHitlists(nextflowRscriptProcess):
    "Look for biological functions possibly associated to hitlists."

    def channel_pretreat(self):
        return [
            [
                "hitlists",
                "flatten()",
                "map{ it.text }",
                "map{ it.split('\\t') }",
                "map{ [it[0].replace(' ', '_'), it[1]] }",
                "set{hitlist}",
            ],
            [
                "Channel",
                "fromPath(params.command_r)",
                "set{command_r}",
            ],
        ]

    def directives(self):
        return {"publishDir": "'../figures', mode: 'copy'"}

    def dependencies(self):
        return {
            "inhouse_packages": [self.rcaller_script],
        }

    def channel_specifications(self):
        return {
            "command_r": chs("each", "file('commandR.r')", "--commandRpath"),
            "hitlist": chs("tuple", ("go_tag", "hitlist"), ("-o", "-i")),
            "go_prefix": chs("val", "go_prefix", "-p", derive_from_params=True),
            "species": chs("val", "species", "--msig_species", derive_from_params=True),
            "go_mcat": chs(
                "val", "go_mcat", "--msig_category", derive_from_params=True
            ),
            "go_msubcat": chs(
                "val", "go_msubcat", "--msig_subcategory", derive_from_params=True
            ),
            "adjust_method": chs(
                "val", "adjust_method", "--pAdjustMethod", derive_from_params=True
            ),
            "qcutoff": chs("val", "qcutoff", "--qvalueCutoff", derive_from_params=True),
            "go_nb_template": chs(
                "each", "go_nb_template", None, derive_from_params=True
            ),
            "go_title": chs("val", "go_title", None, derive_from_params=True),
            "go_script_path": chs(
                "each", "go_script_path", None, derive_from_params=True
            ),
            "order_of_conditions": chs(
                "each", "order_of_conditions", None, derive_from_params=True
            ),
            "go_notebook_pars": chs(
                "tuple",
                (
                    "go_nb_template",
                    "go_script_path",
                    "hitlist",
                    "order_of_conditions",
                    "go_tag",
                    "go_prefix",
                    "species",
                    "go_mcat",
                    "go_msubcat",
                ),
                (None, None, None, None, None, None, None, None, None),
            ),
            "go_figures": chs(
                "tuple",
                (
                    "go_title",
                    "go_tag",
                    '"${go_prefix}${go_tag.replace(' + "'\\\"', ''" + ')}.pdf"',
                ),
                (None, None, None),
            ),
        }

    def customize_features(self):
        node_params = dict(rcaller_script=rcaller_script)
        node_params.update(self.node_params)
        self.rcaller_script = node_params["rcaller_script"]
        return


class visualiseGseaGeneSets(nextflowRscriptProcess):
    "Show gene sets differentially affected by conditions using GSEA."

    def channel_pretreat(self):
        return [
            [
                "Channel",
                "fromPath(params.command_r)",
                "set{command_r2}",
            ],
            [
                "gsea_tabs",
                "replaceAll('\\n', '')",
                "into{gsea_tabs; positional_tabs; kegg_tabs; de_tabs_ko}",
            ],
        ]

    def directives(self):
        return {"publishDir": "'../figures', mode: 'copy'"}

    def dependencies(self):
        return {
            "inhouse_packages": [self.rcaller_script],
        }

    def channel_specifications(self):
        return {
            "command_r2": chs("each", "file('commandR.r')", "--commandRpath"),
            "gsea_nb_template": chs(
                "each", "gsea_nb_template", None, derive_from_params=True
            ),
            "gsea_tabs": chs("each", "gsea_tabs", "-i"),
            "gsea_prefix": chs("each", "gsea_prefix", "-p", derive_from_params=True),
            "gsea_script_path": chs(
                "each", "gsea_script_path", None, derive_from_params=True
            ),
            "order_of_conditions": chs(
                "each",
                "order_of_conditions",
                "--conditionOrder",
                derive_from_params=True,
            ),
            "species": chs(
                "each", "species", "--msig_species", derive_from_params=True
            ),
            "gsea_sets": chs(
                "tuple",
                ("gsea_tag", "gsea_title", "mcat", "msubcat"),
                ("-o", None, "--msig_category", "--msig_subcategory"),
                derive_from_params=True,
            ),
            "gsea_notebook_pars": chs(
                "tuple",
                (
                    "gsea_nb_template",
                    "gsea_script_path",
                    "gsea_tabs",
                    "order_of_conditions",
                    "gsea_tag",
                    "gsea_prefix",
                    "species",
                    "mcat",
                    "msubcat",
                ),
                (None, None, None, None, None, None, None, None, None),
            ),
            "gsea_figures": chs(
                "tuple",
                (
                    "gsea_title",
                    "gsea_tag",
                    '"${gsea_prefix}${gsea_tag.replace(' + "'\\\"', ''" + ')}.pdf"',
                ),
                (None, None, None),
            ),
        }

    def customize_features(self):
        node_params = dict(rcaller_script=rcaller_script)
        node_params.update(self.node_params)
        self.rcaller_script = node_params["rcaller_script"]
        return


class visualisePositionalOverrep(nextflowRscriptProcess):
    "Show chromosomal positions differentially affected by conditions using GSEA."

    def channel_pretreat(self):
        return [
            [
                "Channel",
                "fromPath(params.command_r)",
                "set{command_r3}",
            ],
            [
                "Channel",
                "fromPath(params.gsea_r)",
                "set{gsea_import}",
            ],
        ]

    def directives(self):
        return {"publishDir": "'../figures', mode: 'copy'"}

    def dependencies(self):
        return {
            "inhouse_packages": [self.rcaller_script, self.gsea_script],
        }

    def channel_specifications(self):
        return {
            "command_r3": chs("each", "file('commandR.r')", "--commandRpath"),
            "gsea_import": chs("each", "file('gsea_ontologies.r')", None),
            "positional_nb_template": chs(
                "each", "positional_nb_template", None, derive_from_params=True
            ),
            "positional_title": chs(
                "each", "positional_title", "--plot_title", derive_from_params=True
            ),
            "positional_tabs": chs("each", "positional_tabs", "-i"),
            "positional_prefix": chs(
                "each", "positional_prefix", "-o", derive_from_params=True
            ),
            "positional_script_path": chs(
                "each", "positional_script_path", None, derive_from_params=True
            ),
            "order_of_conditions": chs(
                "each",
                "order_of_conditions",
                "--conditionOrder",
                derive_from_params=True,
            ),
            "species": chs(
                "each", "species", "--msig_species", derive_from_params=True
            ),
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
                (None, None, None, None, None, None, None, None, None),
            ),
            "positional_figures": chs(
                "tuple",
                ("positional_title", "positional_prefix", '"${positional_prefix}.pdf"'),
                (None, None, None),
            ),
        }

    def customize_features(self):
        node_params = dict(rcaller_script=rcaller_script, gsea_script=gsea_script)
        node_params.update(self.node_params)
        self.rcaller_script = node_params["rcaller_script"]
        self.gsea_script = node_params["gsea_script"]
        return


class drawPathViews(nextflowRscriptProcess):
    "Show gene expression changes on KEGG pathway maps."

    def channel_pretreat(self):
        return [
            [
                "Channel",
                "fromPath(params.command_r)",
                "set{command_r3}",
            ],
            [
                "kegg_tabs",
                "map{it.replace('" + '\\"' + "', '')}",
                "map{it.split(',')}",
                "flatten()",
                "map{it.split(':')}",
                "map{['\"' + it[0] + '\"', '\"' + it[0] + '\"', '\"' + it[1] + '\"']}",
                "set{kegg_de_tabs}",
            ],
        ]

    def directives(self):
        return {"publishDir": "'../figures', mode: 'copy'"}

    def dependencies(self):
        return {
            "inhouse_packages": [self.rcaller_script],
        }

    def channel_specifications(self):
        return {
            "command_r3": chs("each", "file('commandR.r')", "--commandRpath"),
            "kegg_path_nb": chs("each", "kegg_path_nb", None, derive_from_params=True),
            "kegg_de_tabs": chs(
                "tuple",
                ("condition_prefix", "path_plot_title", "kegg_de_tabs"),
                ("-o", "--plot_title", "-i"),
            ),
            "kegg_paths": chs("each", "kegg_paths", "-p", derive_from_params=True),
            "pathview_script": chs(
                "each", "pathview_script", None, derive_from_params=True
            ),
            "species": chs("each", "species", "--species", derive_from_params=True),
            "pathview_notebook_pars": chs(
                "tuple",
                (
                    "kegg_path_nb",
                    "pathview_script",
                    "kegg_de_tabs",
                    "kegg_paths",
                    "condition_prefix",
                    "species",
                ),
                (None, None, None, None, None, None),
            ),
            "pathview_figures": chs(
                "tuple",
                ("kegg_paths", "condition_prefix", "'*.pdf'"),
                (None, None, None),
            ),
        }

    def customize_features(self):
        node_params = dict(rcaller_script=rcaller_script)
        node_params.update(self.node_params)
        self.rcaller_script = node_params["rcaller_script"]
        return
