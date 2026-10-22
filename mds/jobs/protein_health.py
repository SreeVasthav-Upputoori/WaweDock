"""RDKit-independent worker for read-only protein health screening."""
from pathlib import Path
from mds.chem.health import check_protein
from mds.jobs import Job


class ProteinHealthJob(Job):
    def __init__(self, source: Path):
        self.source = Path(source)
        self.name = f'Protein health: {self.source.name}'

    def run(self, ctx):
        ctx.check_cancelled()
        return check_protein(self.source, progress=ctx.progress, check_cancelled=ctx.check_cancelled)
