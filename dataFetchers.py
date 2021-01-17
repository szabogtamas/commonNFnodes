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
        self.command = "fastq-dump -I -O ./tmp -X 500 --split-files --disable-multithreading $sample\n            "
        self.command += "mv tmp/${sample}_1.fastq ${sample}_1.fastq\n            "
        self.command += "mv tmp/${sample}_2.fastq ${sample}_2.fastq"
        return None


class downloadPairedSRAfast(nextflowCmdProcess):
    "Download archive sequencing data from SRA database."

    def directives(self):
        return {"publishDir": "params.input_folder", "label": "'manycpu'"}

    def customize_features(self):
        self.inputs = [
            "val manycpu from params.manycpu",
            "tuple sample from params.sample_ids",
        ]
        self.outputs = [
            'tuple "${sample}_1.fastq", "${sample}_2.fastq", sample into insamples'
        ]
        self.command = "fasterq-dump $sample --split-files --threads $manycpu -t ./tmp\n            "
        self.command += "mv tmp/${sample}_1.fastq ${sample}_1.fastq\n            "
        self.command += "mv tmp/${sample}_2.fastq ${sample}_2.fastq"
        return None


class downloadPairedSRAbuiltin(nextflowCmdProcess):
    "Download archive sequencing data from SRA database."

    def directives(self):
        return {"publishDir": "params.input_folder"}

    def channel_pretreat(self):
        return [
            [
                "Channel",
                "fromSRA(params.sample_ids)",
                "set{sra_samples}",
            ],
        ]

    def customize_features(self):
        self.inputs = [
            "val indir from params.input_folder",
            "tuple sample, fq from sra_samples",
        ]
        self.outputs = [
            'tuple "${fq[0]}", "${fq[1]}", sample into insamples' "stdout sra_succes"
        ]
        self.command += "mv ${fq[0]} ${indir}/${fq[0]}\n            "
        self.command += "mv ${fq[1]} ${indir}/${fq[1]}\n            "
        self.command = "echo $sample"


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
