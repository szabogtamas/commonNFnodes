from .nodePrototypes import nextflowRscriptProcess
from .nodePrototypes import introSpect

chs = introSpect.flowNodes.createChannelSpecification
rcaller_script = "commandR.r"


class getDEwithEdgeR(nextflowRscriptProcess):
    "Run an R script using EdgeR to get DE genes."

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
            "publishDir": """'.', saveAs: { it.contains('.tsv') || it.contains('.xlsx') ? "../tables/$it" : "../figures/$it" }, mode: 'copy'"""
        }

    def dependencies(self):
        return {
            "inhouse_packages": [self.rcaller_script],
        }

    def channel_specifications(self):
        return {
            "command_r": chs("file", "'commandR.r'", "--commandRpath"),
            "count_matrix": chs("file", "count_matrix", "-i"),
            "conditions": chs("val", "conditions", "-l", derive_from_params=True),
            "order_of_conditions": chs(
                "val",
                "order_of_conditions",
                "--conditionOrder",
                derive_from_params=True,
            ),
            "edger_prefix": chs("val", "edger_prefix", "-o", derive_from_params=True),
            "edger_cluster": chs(
                "val", "edger_cluster", "--clusterSamples", derive_from_params=True
            ),
            "edger_labelvolc": chs(
                "val", "edger_labelvolc", "--labelVolcano", derive_from_params=True
            ),
            "edger_title": chs("val", "edger_title", None, derive_from_params=True),
            "edger_nb_template": chs(
                "each", "edger_nb_template", None, derive_from_params=True
            ),
            "edger_script_path": chs(
                "each", "edger_script_path", None, derive_from_params=True
            ),
            "gex_norm_mat": chs("file", "'*normalized_matrix.tsv'", None),
            "edger_notebook_pars": chs(
                "tuple",
                (
                    "edger_nb_template",
                    "edger_script_path",
                    "file(count_matrix)",
                    "order_of_conditions",
                ),
                (None, None, None, None),
            ),
            "de_tables": chs("file", "'*.tsv'", None),
            "de_sheet": chs("file", "'*.xlsx'", None),
            "de_figures": chs("tuple", ("edger_title", "'*.pdf'"), (None, None)),
        }

    def customize_features(self):
        node_params = dict(
            rcaller_script=rcaller_script,
        )
        node_params.update(self.node_params)
        self.rcaller_script = node_params["rcaller_script"]
        return
