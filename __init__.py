import os

__doc__ = """
Extending the FlowNodes script, this package contains Python objects that can be used as
drop-in nodes when building a NextFlow pipeline from within Python. Typical use case is
to assemble the nodes in Jupyter notebook, to facilitate record keeping.
"""

homefolder = os.path.abspath(os.path.dirname(__file__))

for module in os.listdir(homefolder):
    if module == "__init__.py" or module[-3:] != ".py":
        continue
    __import__(module[:-3], globals(), locals(), level=1)

del module

nextflowCmdProcess = nodePrototypes.nextflowCmdProcess
nextflowRscriptProcess = nodePrototypes.nextflowRscriptProcess
