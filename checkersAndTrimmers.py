from .nodePrototypes import nextflowCmdProcess

wait_to_fetch = False


class rnaBulkTrimmomaticSE(nextflowCmdProcess):
    "Trim Single End reads with Trimmomatic."

    def channel_pretreat(self):
        if self.wait_to_fetch:
            p = dict()
        else:
            p = [
                [
                    "Channel",
                    "from([params.raw_fastqs, params.sample_ids].transpose())",
                    "map{[prams.input_folder + '/' + it[0], it[1]]}",
                    "set{insamples}",
                ],
            ]
        return p

    def customize_features(self):
        node_params = dict(
            wait_to_fetch=wait_to_fetch,
        )
        node_params.update(self.node_params)
        self.wait_to_fetch = node_params["wait_to_fetch"]
        self.inputs = [
            "val adapterFileIllumina from params.adapterFileIllumina",
            "val manycpu from params.manycpu",
            "val indir from params.input_folder",
            "tuple raw, sample from insamples",
        ]
        self.outputs = ['tuple sample, "${sample}_trimed.fastq" into trimmed_fastqs']
        self.command = "trimmomatic SE $raw\\\n            "
        self.command += "${sample}_trimed.fastq\\\n            "
        self.command += "ILLUMINACLIP:${adapterFileIllumina}:2:30:10 LEADING:3 TRAILING:3 MINLEN:36\\\n            "
        self.command += "-threads $manycpu\n"

        return None


class rnaBulkTrimmomaticPE(nextflowCmdProcess):
    "Trim Paired End reads with Trimmomatic."

    def channel_pretreat(self):
        if self.wait_to_fetch:
            p = dict()
        else:
            p = [
                [
                    "Channel",
                    "from([params.forward_fastqs, params.reverse_fastqs, params.sample_ids].transpose())",
                    "map{[prams.input_folder + '/' + it[0], prams.input_folder + '/' + it[1], it[2]]}",
                    "set{insamples}",
                ],
            ]
        return p

    def customize_features(self):
        node_params = dict(
            wait_to_fetch=wait_to_fetch,
        )
        node_params.update(self.node_params)
        self.wait_to_fetch = node_params["wait_to_fetch"]
        self.inputs = [
            "val adapterFileIllumina from params.adapterFileIllumina",
            "val manycpu from params.manycpu",
            "val indir from params.input_folder",
            "tuple forward, reverse, sample from insamples",
        ]
        self.outputs = [
            'tuple sample, "${sample}_trim_1.fastq", "${sample}_trim_2.fastq" into trimmed_fastqs'
        ]
        self.command = (
            "trimmomatic PE ${indir}/${forward} ${indir}/${reverse}\\\n            "
        )
        self.command += (
            "${sample}_trim_1.fastq ${sample}_forward_unpaired.fastq\\\n            "
        )
        self.command += (
            "${sample}_trim_2.fastq ${sample}_reverse_unpaired.fastq\\\n            "
        )
        self.command += "ILLUMINACLIP:$adapterFileIllumina:2:30:10:8:keepBothReads LEADING:3 TRAILING:3 MINLEN:36\\\n            "
        self.command += "-threads $manycpu\n"
        return None
