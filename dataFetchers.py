from .nodePrototypes import nextflowCmdProcess


class downloadPairedSRA(nextflowCmdProcess):
    "Download archive sequencing data from SRA database."

    def directives(self):
        return {"publishDir": "params.input_folder"}

    def customize_features(self):
        self.inputs = [
            "tuple sample from params.sample_ids",
        ]
        self.outputs = [
            'tuple "${sample}_1.fastq", "${sample}_2.fastq", sample into insamples'
        ]
        self.command = "fastq-dump -I --split-files $sample"
        return None
