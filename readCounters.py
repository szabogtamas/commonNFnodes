from .nodePrototypes import nextflowCmdProcess


class countWithFeatureCounts(nextflowCmdProcess):
    "Count reads with FeatureCounts, using general settings."

    def directives(self):
        return {"publishDir": "'../tables', mode: 'copy'", "label": "'manycpu'"}

    def customize_features(self):
        self.inputs = [
            "val manycpu from params.manycpu",
            "val count_file from params.count_file",
            "val genomeannotation from params.genomeannotation",
            "file alignedSams from aligned.collect()",
        ]
        self.outputs = ["file 'counts.tsv' into count_file"]
        self.command = "featureCounts -T $manycpu\\\n            "
        self.command += "-a $genomeannotation\\\n            "
        self.command += "-o $count_file\\\n            "
        self.command += "${alignedSams.join(' ')}\\\n"
        return None


class countSymbolsWithFeatureCounts(nextflowCmdProcess):
    "Count reads with FeatureCounts, using general settings."

    def directives(self):
        return {"publishDir": "'../tables', mode: 'copy'", "label": "'manycpu'"}

    def customize_features(self):
        self.inputs = [
            "val manycpu from params.manycpu",
            "val count_file from params.count_file",
            "val genomeannotation from params.genomeannotation",
            "file alignedSams from aligned.collect()",
        ]
        self.outputs = ["file 'counts.tsv' into count_file"]
        self.command = "featureCounts -T $manycpu\\\n            "
        self.command += "-a $genomeannotation\\\n            "
        self.command += "-o $count_file\\\n            "
        self.command += "-g gene_name\\\n            "
        self.command += "${alignedSams.join(' ')}\\\n"
        return None


class countWithHTSeqCounts(nextflowCmdProcess):
    "Count reads with HTSeqCounts, using general settings."

    def directives(self):
        return {"publishDir": "'../tables', mode: 'copy'", "label": "'manycpu'"}

    def customize_features(self):
        self.inputs = [
            "val manycpu from params.manycpu",
            "val count_file from params.count_file",
            "val genomeannotation from params.genomeannotation",
            "file alignedSams from aligned.collect()",
        ]
        self.outputs = ["file 'counts.tsv' into count_file"]
        self.command = "htseq-count\\\n            "
        self.command += "-o $count_file\\\n            "
        self.command += "${alignedSams.join(' ')}            "
        self.command += "$genomeannotation\\\n"
        return None


class countWithSalmon(nextflowCmdProcess):
    "Count reads with Salmon, using general settings."

    def directives(self):
        return {"publishDir": "'../tables', mode: 'copy'", "label": "'manycpu'"}

    def customize_features(self):
        self.inputs = [
            "val manycpu from params.manycpu",
            "val genomeindex from params.genomeindex" # e.g. /home/szabo/myScratch/SeqmiRNA/indices/salmon
            'tuple sample, "${sample}_trimed.fastq" from trimmed_fastqs',
        ]
            "val manycpu from params.manycpu",
            "val count_file from params.count_file",
            "file alignedSams from aligned.collect()",
        ]
        self.outputs = ["file '${alignedSams}_counts.tsv' into count_file"]
        self.command = "salmon quant\\\n            "
        self.command += "-p $manycpu -l A\\\n            "
        self.command += "--validateMappings\\\n            "
        self.command += "-r ${alignedSams.join(' ')}\\\n            "
        self.command += "-o $count_file\\\n"
        return None

    
class countWithSalmonAligned(nextflowCmdProcess):
    "Count reads with Salmon, using general settings and alignment mode."

    def directives(self):
        return {"publishDir": "'../tables', mode: 'copy'", "label": "'manycpu'"}

    def customize_features(self):
        self.inputs = [
            "val manycpu from params.manycpu",
            "val genomeindex from params.genomeindex" # e.g. salmon.fa
            "file alignedSams from aligned.collect()",
        ]
            "val manycpu from params.manycpu",
            "val count_file from params.count_file",
            "file alignedSams from aligned.collect()",
        ]
        self.outputs = ["file 'counts.tsv' into count_file"]
        self.command = "salmon quant\\\n            "
        self.command += "-t $genomeindex -l A\\\n            "
        self.command += "-p $manycpu -l A\\\n            "
        self.command += "-a ${sample}_trimed.fastq\\\n            "
        self.command += "-o $count_file\\\n"
        return None
