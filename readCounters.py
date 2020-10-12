from . import nextflowCmdProcess


class countWithFeatureCounts(nextflowCmdProcess):
    "Count reads with FeatureCounts, using general settings."

    def directives(self):
        return {"publishDir": "'../tables', mode: 'copy'", "label": "'manycpu'"}

    def customize_features(self):
        self.inputs = [
            "val manycpu from params.manycpu",
            "val genomeannotation from params.genomeannotation",
            "file alignedSams from aligned.collect()",
        ]
        self.outputs = ["file 'counts.tsv'"]
        self.command = "featureCounts -T $manycpu\\\n            "
        self.command += "-a $genomeannotation\\\n            "
        self.command += "-o counts.tsv\\\n            "
        self.command += "${alignedSams.join(' ')}\\\n"
        return None
