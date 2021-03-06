import os, sys, subprocess

homefolder = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.dirname(homefolder))

import introSpect

nextflowProcess = introSpect.flowNodes.nextflowProcess

generalClusterProfile = """
profiles {
    standard {
        process.executor = 'local'
    }
    cluster {
        process.executor = 'sge'
        process.penv = 'smp'
        process.clusterOptions = { '-V -S /bin/bash -q all.q@apollo-*' }
        process.cpus = 1
    }
}
"""


class nextflowCmdProcess(nextflowProcess):
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
        command = self.command
        deps = self.dependencies()
        if "inhouse_packages" in deps and self.command_locally:
            for r_script in deps["inhouse_packages"]:
                subprocess.call(["cp", r_script, dr + "/bin/" + os.path.basename(r_script)])
            subprocess.call(["cp", command, dr + "/bin/" + os.path.basename(command)])
            self.command = dr + "/bin/" + os.path.basename(command)

    def compile_command(self):
        return (
            "Rscript "
            + self.command
            + "\\\n            "
            + "\\\n            ".join(
                ["" if flag is None else flag for flag in self.flags]
            )
            + "\n"
        )
