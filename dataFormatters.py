import os
import pandas as pd

from .nodePrototypes import nextflowProcess
from .nodePrototypes import introSpect

chs = introSpect.flowNodes.createChannelSpecification
intro_folder = os.path.dirname(os.path.realpath(introSpect.__file__))
rerun_preprocessing = True


def safe_quote(s):
    return '"' + s + '"'


def meta_parser(sample_meta_file):
    return pd.read_csv(sample_meta_file, sep="\t").set_index("SampleID")


class countMatrixFormatter(nextflowProcess):
    def channel_pretreat(self):
        return [
            [
                "Channel",
                "from('" + '"' + "' + params.sample_ids.join(',') + '" + '"' + "')",
                "set{sample_ids}",
            ],
        ]

    def dependencies(self):
        return {
            "imports": ["import pandas as pd"],
            "helpers": [self.meta_parser],
            "inhouse_packages": [self.intro_folder],
        }

    def channel_specifications(self):
        return {
            "count_file": chs("file", "count_file", "count_file")
            if self.rerun_preprocessing
            else chs("val", "count_file", "count_file", derive_from_params=True),
            "genedict_file": chs("val", "genedict_file", "genedict_file"),
            "sample_meta_file": chs("val", "sample_meta_file", "sample_meta_file"),
            "sample_ids": chs("val", "sample_ids", "sample_ids"),
            "count_matrix": chs("file", "'count_matrix.tsv'", "outFile"),
            "count_matrix2": chs("file", "'*.tsv'", None),
        }

    def process(
        self,
        count_file: str,
        genedict_file: str,
        sample_meta_file: str,
        sample_ids: list,
    ) -> pd.DataFrame:
        """
        Remove the first few columns of the count matrix not containing counts, but add gene symbols.

        Parameters
        ----------
        count_file
            Count matrix as it comes from FeatureCounts output.
        genedict_file
            A mapping between ENEMBL IDs and HGNC gene symbols.
        sample_meta_file
            Sample metadata, especially the experimental condition.
        sample_ids
            Sample IDs in the same order as their labels

        Returns
        -------
        A dataframe with unique ID in the first column, human-friendly second and counts in the rest.
        """

        geneDict = (
            pd.read_csv(genedict_file, sep="\t")
            .set_index("Gene name")
            .to_dict()["Gene stable ID"]
        )
        meta_df = meta_parser(sample_meta_file)

        df = pd.read_csv(count_file, sep="\t", skiprows=1)
        df["Geneid2"] = df["Geneid"].apply(
            lambda x: geneDict[x] if x in geneDict else x
        )
        df = df.rename(columns={"Geneid2": "Genesym"})
        df = df.set_index("Genesym")
        df = df.loc[:, ["Geneid"] + [x + "_Aligned.out.sam" for x in sample_ids]]
        df.rename(
            columns={
                x: meta_df.loc[x.replace("_Aligned.out.sam", "")]["Sample"]
                for x in set(df.columns) - set(["Geneid"])
            },
            inplace=True,
        )

        return df

    def customize_features(self):
        node_params = dict(
            intro_folder=intro_folder,
            rerun_preprocessing=rerun_preprocessing,
            meta_parser=meta_parser,
        )
        node_params.update(self.node_params)
        self.intro_folder = node_params["intro_folder"]
        self.rerun_preprocessing = node_params["rerun_preprocessing"]
        self.meta_parser = node_params["meta_parser"]
        return


class hitListCreator(nextflowProcess):
    def channel_pretreat(self):
        return [
            [
                "de_tables",
                "collect()",
                "into{de_tabs; de_tabs2; de_tabs3}",
            ],
        ]

    def dependencies(self):
        return {
            "imports": ["import os", "import pandas as pd"],
            "helpers": [safe_quote],
            "inhouse_packages": [self.intro_folder],
        }

    def channel_specifications(self):
        return {
            "de_tabs": chs("file", "de_tabs", "*de_tables"),
            "hitlists": chs("file", "'*_hitlist.tsv'", None),
        }

    def process(
        self,
        edger_prefix: str,
        *de_tables,
    ) -> list:
        """
        Take topTag DE tables and create hitlists by selecting significantly UP or DOWN regulated genes.

        Parameters
        ----------
        edger_prefix
            Prefix of files added by EdgeR. This will be removed from labels.
        de_tables
            List of paths to de_tables created in the previous step.

        Returns
        -------
        A list of hitlists that will be further analyzed for GO terms.
        """

        hls = []
        bm = []
        bms = []
        de_tabs = []
        for d in de_tables:
            d = os.path.basename(d)
            if d.find(edger_prefix) > -1:
                if d != edger_prefix + "normalized_matrix.tsv":
                    de_tabs.append(d)
        for tbl in de_tabs:
            tbn = (
                tbl.replace(edger_prefix, "").replace(".tsv", "").replace("___", "_")
            )  # , "/")
            df = pd.read_csv(tbl, sep="\t")
            upreg = ",".join(
                df.loc[
                    (df["logFC"] > 1) & (df["FDR"] < 0.05) & df["Symbol"].notna(),
                    "Symbol",
                ].astype(str)
            )
            downreg = ",".join(
                df.loc[
                    (df["logFC"] < -1) & (df["FDR"] < 0.05) & df["Symbol"].notna(),
                    "Symbol",
                ].astype(str)
            )
            hl = tbn + "_UP," + upreg + ":" + tbn + "_DOWN," + downreg
            hls.append([tbn, hl])
            bm.append(hl)
            bms.append(",".join([tbn, upreg, downreg]))
        hls.append(["overview_by_regulation", ":".join(bm)])
        hls.append(["overview_by_condition", ":".join(bms)])
        for i, h in enumerate(hls):
            h = [safe_quote(x) for x in h]
            with open("n" + str(i) + "_hitlist.tsv", "w") as f:
                f.write("\t".join(h))
        return None


class gseaTabCreator(nextflowProcess):
    def dependencies(self):
        return {
            "imports": ["import os", "import pandas as pd"],
            "helpers": [safe_quote],
            "inhouse_packages": [self.intro_folder],
        }

    def channel_specifications(self):
        return {
            "de_tabs2": chs("file", "de_tabs", "*de_tables"),
            "gsea_tabs": chs("stdout", "", None),
        }

    def process(
        self,
        edger_prefix: str,
        *de_tables,
    ) -> list:
        """
        Take topTag DE tables and label them for GSEA.

        Parameters
        ----------
        edger_prefix
            Prefix of files added by EdgeR. This will be removed from labels.
        de_tables
            List of paths to de_tables created in the previous step.

        Returns
        -------
        An expression with condition labels and table paths to be fed into GSEA script.
        """

        de_tabs = []
        gsea_in = ""
        for d in de_tables:
            d = os.path.basename(d)
            if d.find(edger_prefix) > -1:
                if d != edger_prefix + "normalized_matrix.tsv":
                    de_tabs.append(d)
        for tbl in de_tabs:
            tbn = tbl.replace(edger_prefix, "").replace(".tsv", "").replace("___", "/")
            gsea_in += tbn + ":" + os.path.abspath(tbl) + ","
        return '"' + gsea_in[:-1] + '"'
