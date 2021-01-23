from .nodePrototypes import nextflowCmdProcess


class alignPairedEndReadsWithStar(nextflowCmdProcess):
    "Use STAR to align paired end reads with general settings."

    def customize_features(self):
        self.inputs = [
            "val manycpu from params.manycpu",
            "val genomedir from params.genomedir",
            'tuple sample, "${sample}_trim_1.fastq", "${sample}_trim_2.fastq" from trimmed_fastqs',
        ]
        self.outputs = ['file "${sample}_Aligned.out.sam" into aligned']
        self.command = "STAR --runMode alignReads\\\n            "
        self.command += "--quantMode GeneCounts\\\n            "
        self.command += "--runThreadN $manycpu\\\n            "
        self.command += "--genomeDir $genomedir\\\n            "
        self.command += "--outFileNamePrefix ${sample}_\\\n            "
        self.command += (
            "--readFilesIn ${sample}_trim_1.fastq ${sample}_trim_2.fastq\\\n"
        )
        return None


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
        return None


class alignSingleEndReadsWithBWA(nextflowCmdProcess):
    "Use BWA to align single end reads with general settings."

    def customize_features(self):
        self.inputs = [
            "val manycpu from params.manycpu",
            "val genomedir from params.genomedir",
            "val genomeindex from params.genomeindex",
            'tuple sample, "${sample}_trimed.fastq" from trimmed_fastqs',
        ]
        self.outputs = ['file "${sample}_Aligned.out.sam" into aligned']
        self.command = "bwa aln\\\n            "
        self.command += "--o 0\\\n            "
        self.command += "-t $manycpu\\\n            "
        self.command += "$genomeindex\\\n            "
        self.command += "${sample}_trimed.fastq\\\n            "
        self.command += "> ${sample}.sai\n            "
        self.command += "bwa samse\\\n            "
        self.command += "$genomeindex\\\n            "
        self.command += "${sample}.sai\n            "
        self.command += "${sample}_trimed.fastq\\\n            "
        self.command += "> ${sample}_\\\n            "
        return None


class alignSingleEndSmallRNAReadsWithBWA(nextflowCmdProcess):
    "Use BWA to align single end reads with settings optimized for miRNA."

    def customize_features(self):
        self.inputs = [
            "val manycpu from params.manycpu",
            "val genomedir from params.genomedir",
            "val genomeindex from params.genomeindex",
            'tuple sample, "${sample}_trimed.fastq" from trimmed_fastqs',
        ]
        self.outputs = ['file "${sample}_Aligned.out.sam" into aligned']
        self.command = "bwa aln\\\n            "
        self.command += "-t $manycpu\\\n            "
        self.command += "$genomeindex\\\n            "
        self.command += "${sample}_trimed.fastq\\\n            "
        self.command += "> ${sample}.sai\n            "
        self.command += "bwa samse\\\n            "
        self.command += "$genomeindex\\\n            "
        self.command += "${sample}.sai\n            "
        self.command += "${sample}_trimed.fastq\\\n            "
        self.command += "> ${sample}_\\\n            "
        return None
