import os, sys, subprocess

homefolder = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.dirname(homefolder))

import introSpect


class nextflowCmdProcess(introSpect.flowNodes.nextflowProcess):
    def directives(self):
        return {"label": "'manycpu'"}

    def compile_command(self):
        return self.command


class nextflowRscriptProcess(nextflowCmdProcess):
    def compile_process(self, dr):
        self.manualDoc = self.__doc__
        self.params = {v[2]: v[1] for k, v in self.channel_specifications().items()}
        self.cmdpars = None
        self.cmdouts = dict()
        for r_script in self.dependencies["inhouse_packages"]:
            subprocess.call(["cp", r_script, dr + "/bin/" + os.path.basename(r_script)])
        subprocess.call(
            ["cp", self.command, dr + "/bin/" + os.path.basename(self.command)]
        )

    def compile_command(self):
        return (
            "Rscript bin/"
            + os.path.basename(self.command)
            + "\\\n            "
            + "\\\n            ".join(self.flags)
            + "\n"
        )
