from . import nextflowCmdProcess


class rnaBulkTrimmomaticSE(nextflowCmdProcess):
    "Trim Single End reads with Trimmomatic."

    def channel_pretreat(self):
        return [
            [
                "Channel",
                "from([params.raw_fastqs, params.sample_ids].transpose())",
                "set{insamples}",
            ],
        ]

    def customize_features(self):
        self.inputs = [
            "val adapterFileIllumina from params.adapterFileIllumina",
            "val manycpu from params.manycpu",
            "val indir from params.input_folder",
            "tuple raw, sample from insamples",
        ]
        self.outputs = ['tuple sample, "${sample}_trimed.fastq" into trimmed_fastqs']
        self.command = "trimmomatic SE ${indir}/${raw}\\\n            "
        self.command += "${sample}_trimed.fastq\\\n            "
        self.command += "ILLUMINACLIP:${adapterFileIllumina}:2:30:10 LEADING:3 TRAILING:3 MINLEN:36\\\n            "
        self.command += "-threads $manycpu\n"

        return None


class rnaBulkTrimmomaticPE(nextflowCmdProcess):
    "Trim Paired End reads with Trimmomatic."

    def channel_pretreat(self):
        return [
            [
                "Channel",
                "from([params.raw_fastqs, params.sample_ids].transpose())",
                "set{insamples}",
            ],
        ]

    def customize_features(self):
        self.inputs = [
            "val adapterFileIllumina from params.adapterFileIllumina",
            "val manycpu from params.manycpu",
            "val indir from params.input_folder",
            "tuple raw, sample from insamples",
        ]
        self.outputs = [
            'tuple sample, "${sample}_trim_1.fastq", "${sample}_trim_2.fastq" into trimmedFastq'
        ]
        self.command = "trimmomatic PE ${indir}/${sample}_1.fastq ${indir}/${sample}_2.fastq\\\n            "
        self.command += (
            "${sample}_trim_1.fastq ${sample}_forward_unpaired.fastq\\\n            "
        )
        self.command += (
            "${sample}_trim_2.fastq ${sample}_reverse_unpaired.fastq\\\n            "
        )
        self.command += "ILLUMINACLIP:$adapterFileIllumina:2:30:10:8:keepBothReads LEADING:3 TRAILING:3 MINLEN:36\\\n            "
        self.command += "-threads $manycpu\n"
        return None
