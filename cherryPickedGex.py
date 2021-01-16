import os
import pandas as pd
import seaborn as sns

from .nodePrototypes import nextflowProcess
from .nodePrototypes import introSpect

chs = introSpect.flowNodes.createChannelSpecification
intro_folder = os.path.dirname(os.path.realpath(introSpect.__file__))


def meta_parser(sample_meta_file):
    return pd.read_csv(sample_meta_file, sep="\t").set_index("SampleID")


def show_expression_of_KO_genes(
    sample_meta_file: str,
    normalized_matrix: str,
    ko_list: list,
    ko_dict: dict = "",
    percentile: bool = False,
    heat: bool = False,
):
    """
    Show expression of the genes that were actually knocked out.

    Parameters
    ----------
    sample_meta_file
        File containing metadata.
    normalized_matrix
        Normalized gene expression matrix.
    ko_list
        List of genes that were knocked down.
    ko_dict
        Mapping of human friendly names for genes.
    percentile
        If the percentile of gene expression should be shown instead of CPM.
    heat
        If heatmap representation is required.

    Returns
    -------
    A boxplot corresponding distribution of gene expression for each KO gene.
    """

    if ko_dict in [None, "None", ""]:
        ko_dict = dict()
    if isinstance(ko_dict, str):
        kos = ko_dict.split(",")
        ko_dict = dict()
        for e in kos:
            e2 = e.split(":")
            ko_dict[e2[0]] = e2[-1]
    meta_df = meta_parser(sample_meta_file).loc[:, "Label"].to_dict()
    norm_df = pd.read_csv(normalized_matrix, sep="\t")
    pct_dct = {
        "True": True,
        "False": False,
        "true": True,
        "false": False,
        "1": True,
        "0": False,
        1: True,
        0: False,
        None: False,
        "None": False,
        "": False,
    }
    if percentile in pct_dct:
        percentile = pct_dct[percentile]
    if heat in pct_dct:
        heat = pct_dct[heat]
    if percentile:
        norm_df = norm_df.set_index(["GeneID", "Symbol"])
        norm_df = norm_df.rank(axis=0, pct=True, numeric_only=True)
        norm_df = norm_df.reset_index()
    norm_df["Symbol"] = norm_df["Symbol"].apply(
        lambda x: ko_dict[x] if x in ko_dict else x
    )
    ko_df = norm_df.loc[norm_df["Symbol"].isin(ko_list), :]
    fig, ax = plt.subplots(figsize=(7.2, 3.6))
    if heat:
        ko_df = ko_df.drop(columns=["GeneID"]).set_index(["Symbol"])
        ko_df = ko_df + 1
        ko_df = ko_df.apply(np.log2)
        ko_df = ko_df.loc[[ko_dict[x] if x in ko_dict else x for x in ko_list], :]
        ax = sns.heatmap(ko_df, annot=True, ax=ax)
        ax.set_xlabel("")
        if percentile:
            ax.set_title("Percentile of expressed genes")
        else:
            ax.set_title("Log2(CPM) of genes")
    else:
        ko_df = ko_df.melt(id_vars=["GeneID", "Symbol"])
        ko_df["variable"] = ko_df["variable"].map(meta_df)
        if percentile:
            ylabel = "Percentile of gene expression"
        else:
            ylabel = "Gene expression (logCPM)"
        ax = sns.boxplot(
            x="Symbol",
            y="value",
            hue="variable",
            order=[ko_dict[x] if x in ko_dict else x for x in ko_list],
            data=ko_df,
            linewidth=0.5,
            # color="white",
            fliersize=0.5,
            ax=ax,
        )
        ax.set_xticklabels(
            [item.get_text() for item in ax.get_xticklabels()], rotation=30, ha="right"
        )
        ax.set_xlabel("")
        ax.set_ylabel(ylabel)
    return ax


class expressionOfKOs(nextflowProcess):
    "Show expression of the KO genes."

    def channel_pretreat(self):
        return [
            [
                "Channel",
                "from('" + '"' + "' + params.ko_list.join(',') + '" + '"' + "')",
                "set{ko_list}",
            ],
            [
                "Channel",
                "from('" + '"' + "' + params.ko_dict.join(',') + '" + '"' + "')",
                "set{ko_dict}",
            ],
        ]

    def dependencies(self):
        return {
            "imports": [
                "import os",
                "import numpy as np",
                "import pandas as pd",
                "import seaborn as sns",
                "from matplotlib import pyplot as plt",
            ],
            "helpers": [self.meta_parser],
            "inhouse_packages": [self.intro_folder],
        }

    def directives(self):
        return {"publishDir": "'../figures', mode: 'copy'"}

    def channel_specifications(self):
        return {
            "sample_meta_file": chs(
                "each", "'sample_meta.tsv'", "sample_meta_file", derive_from_params=True
            ),
            "gex_norm_mat": chs("each", "file('gex.tsv')", "normalized_matrix"),
            "ko_title": chs("each", "ko_title", "ko_title", derive_from_params=True),
            "ko_list": chs("each", "ko_list", "ko_list"),
            "ko_dict": chs("each", "ko_dict", "ko_dict"),
            "ko_gex_figures": chs(
                "tuple",
                ("ko_sub", "ko_gex_figure", "ko_gex_percentile", "ko_gex_heat"),
                (None, None, "percentile", "heat"),
                derive_from_params=True,
            ),
            "kogex": chs(
                "tuple",
                ("ko_title", "ko_sub", "file(ko_gex_figure)"),
                (None, None, "outFile"),
            ),
        }

    def customize_features(self):
        node_params = dict(
            intro_folder=intro_folder,
            meta_parser=meta_parser,
        )
        node_params.update(self.node_params)
        self.intro_folder = node_params["intro_folder"]
        self.meta_parser = node_params["meta_parser"]
        self.process = show_expression_of_KO_genes
        return


class changeInKOs(nextflowProcess):
    "Show change in expression for KO genes."

    def channel_pretreat(self):
        return [
            [
                "Channel",
                "from('" + '"' + "' + params.ko_list.join(',') + '" + '"' + "')",
                "set{ko_list2}",
            ],
            [
                "Channel",
                "from('" + '"' + "' + params.ko_dict.join(',') + '" + '"' + "')",
                "set{ko_dict2}",
            ],
        ]

    def dependencies(self):
        return {
            "imports": [
                "import os",
                "import pandas as pd",
                "import seaborn as sns",
                "from matplotlib import pyplot as plt",
            ],
            "inhouse_packages": [self.intro_folder],
        }

    def directives(self):
        return {"publishDir": "'../figures', mode: 'copy'"}

    def channel_specifications(self):
        return {
            "de_tabs_ko": chs("each", "de_res_tabs", "de_res_tabs"),
            "ko_title": chs("each", "ko_title", "ko_title", derive_from_params=True),
            "ko_list2": chs("each", "ko_list", "ko_list"),
            "ko_dict2": chs("each", "ko_dict", "ko_dict"),
            "ko_de_figures": chs(
                "tuple",
                ("ko_sub", "ko_de_figure", "ko_de_percentile"),
                (None, None, "percentile"),
                derive_from_params=True,
            ),
            "kode": chs(
                "tuple",
                ("ko_title", "ko_sub", "file(ko_de_figure)"),
                (None, None, "outFile"),
            ),
        }

    def process(
        self,
        de_res_tabs: str,
        ko_list: list,
        *,
        ko_dict: dict = "",
        order_of_conditions: list = None,
        percentile: bool = False
    ):
        """
        Show fold changes for genes that were actually knocked out.

        Parameters
        ----------
        de_res_tabs
            Path to DE tables, with tags for conditions.
        ko_list
            List of genes that were knocked down.
        ko_dict
            Mapping of human friendly names for genes.
        order_of_conditions
            Order of conditions on the plot (and matching colors).
        percentile
            If the percentile of gene expression should be shown instead of CPM.

        Returns
        -------
        A barplot showing log fold change in each condition.
        """

        if ko_dict in [None, "None", ""]:
            ko_dict = dict()
        if isinstance(ko_dict, str):
            kos = ko_dict.split(",")
            ko_dict = dict()
            for e in kos:
                e2 = e.split(":")
                ko_dict[e2[0]] = e2[-1]
        pct_dct = {
            "True": True,
            "False": False,
            "true": True,
            "false": False,
            "1": True,
            "0": False,
            1: True,
            0: False,
            None: False,
            "None": False,
            "": False,
        }
        if percentile in pct_dct:
            percentile = pct_dct[percentile]

        de_tabs = de_res_tabs.split(",")
        dtd = dict()
        for e in de_tabs:
            tag, pth = e.split(":")
            df = pd.read_csv(pth, sep="\t").reset_index()
            df["Symbol"] = df["Symbol"].apply(
                lambda x: ko_dict[x] if x in ko_dict else x
            )
            df = df.set_index("Symbol")
            df = df["logFC"]
            if percentile:
                df = df.rank(pct=True, numeric_only=True)
            dtd[tag] = df.loc[ko_list] - 0.5
        if order_of_conditions is None:
            order_of_conditions = list(dtd.keys())
        ko_df = pd.DataFrame(dtd).reset_index()
        ko_df = ko_df.melt(id_vars=["Symbol"])
        fig, ax = plt.subplots(figsize=(7.2, 3.6))
        if percentile:
            ylabel = "Percentile of Fold Changes"
        else:
            ylabel = "log(Fold Change)"
        ax = sns.barplot(
            x="Symbol",
            y="value",
            hue="variable",
            order=[ko_dict[x] if x in ko_dict else x for x in ko_list],
            data=ko_df,
            linewidth=0.5,
            hue_order=order_of_conditions,
            ax=ax,
        )
        ax.set_xticklabels(
            [item.get_text() for item in ax.get_xticklabels()], rotation=30, ha="right"
        )
        ax.set_title("Fold changes calculated by DE")
        ax.set_xlabel("")
        ax.set_ylabel(ylabel)
        return ax

    def customize_features(self):
        node_params = dict(
            intro_folder=intro_folder,
        )
        node_params.update(self.node_params)
        self.intro_folder = node_params["intro_folder"]
        return
