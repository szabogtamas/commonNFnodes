class rnaBulkTrimmomaticSE(nextflowProcess):

    def customize_features(self):
        self.manualDoc = "Trimm reads with Trimmomatic.\n"
        self.inputs = [
            "file raw from params.raw_fastqs",
        ]
        self.outputs = [
            "tuple sraId, "${sraId}_trim_1.fastq", "${sraId}_trim_2.fastq" into trimmedFastq"
        ]
        self.command = "trimmomatic SE ${indir}/${sraId}_1.fastq ${indir}/${sraId}_2.fastq \\ \n            "
        self.command += "${sraId}_trim_1.fastq tr_u1.fastq \\ \n            "
        self.command += "${sraId}_trim_2.fastq tr_u2.fastq \\ \n            "
        self.command += "ILLUMINACLIP:$adapterFileIllumina:2:30:10:8:keepBothReads LEADING:3 TRAILING:3 MINLEN:36 \\\n            "
        
        return None
    
    def compile_command(self):
        return self.command

class rnaBulkTrimmomaticPE(nextflowProcess):

    def customize_features(self):
        self.manualDoc = "Trimm reads with Trimmomatic.\n"
        self.inputs = [
            "val indir from params.indir",
            "val adapterFileIllumina from params.adapterFileIllumina"
        ]
        self.outputs = [
            "tuple sraId, "${sraId}_trim_1.fastq", "${sraId}_trim_2.fastq" into trimmedFastq"
        ]
        self.command = "trimmomatic PE ${indir}/${sraId}_1.fastq ${indir}/${sraId}_2.fastq \\ \n            "
        self.command += "${sraId}_trim_1.fastq tr_u1.fastq \\ \n            "
        self.command += "${sraId}_trim_2.fastq tr_u2.fastq \\ \n            "
        self.command += "ILLUMINACLIP:$adapterFileIllumina:2:30:10:8:keepBothReads LEADING:3 TRAILING:3 MINLEN:36 \\\n            "
        
        return None
    
    def compile_command(self):
        return self.command