from . import nextflowCmdProcess


class alignSingleEndReadsWithStar(nextflowCmdProcess):
    "Use STAR to align single end reads with general settings."

    def customize_features(self):
        self.inputs = [
            "val manycpu from params.manycpu",
            "val genomedir from params.genomedir",
            'tuple sample, "${sample}_trimed.fastq" from trimmed_fastqs',
        ]
        self.outputs = ['file "${sample}_Aligned.out.sam" into aligned']
        self.command = "STAR --runMode alignReads\\\n            "
        self.command += "--quantMode GeneCounts\\\n            "
        self.command += "--runThreadN $manycpu\\\n            "
        self.command += "--genomeDir $genomedir\\\n            "
        self.command += "--outFileNamePrefix ${sample}_\\\n            "
        self.command += "--readFilesIn ${sample}_trimed.fastq\\\n"

        self.manualDoc = self.__doc__
        return None


class alignPairedEndReadsWithStar(nextflowCmdProcess):
    "Use STAR to align paired end reads with general settings."

    def customize_features(self):
        self.inputs = [
            "val manycpu from params.manycpu",
            "val genomedir from params.genomedir",
            'tuple sample, "${sample}_trim_1.fastq", "${sample}_trim_2.fastq" into trimmedFastqs',
        ]
        self.outputs = ['file "${sample}_Aligned.out.sam" into aligned']
        self.command = "STAR --runMode alignReads\\\n            "
        self.command += "--quantMode GeneCounts\\\n            "
        self.command += "--runThreadN $manycpu\\\n            "
        self.command += "--genomeDir $genomedir\\\n            "
        self.command += "--outFileNamePrefix ${sample}_\\\n            "
        self.command += "readFilesIn ${sample}_trim_1.fastq ${sample}_trim_2.fastq\\\n"

        self.manualDoc = self.__doc__
        return None
