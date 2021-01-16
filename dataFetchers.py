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
        self.command = "fasterq-dump -I -t ./tmp --split-files $sample"
        return None


class localFetchNunzip(nextflowCmdProcess):
    "Copy files from a local folder and unzip them in a work dir."

    def directives(self):
        return {"publishDir": "params.input_folder"}

    def customize_features(self):
        self.inputs = [
            "file arch sample from params.archives",
        ]
        self.outputs = ["file '*.fastq' into insamples"]
        self.command = "unzip $arch"
        return None
