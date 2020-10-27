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
