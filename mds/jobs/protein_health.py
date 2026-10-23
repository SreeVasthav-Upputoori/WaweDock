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


class ProteinRepairPlanJob(Job):
    def __init__(self, source, selected_residues=None, ph=7.4, sequence_references=()):
        self.source,self.selected_residues,self.ph=Path(source),selected_residues,ph
        self.sequence_references=tuple(sequence_references)
        self.name='Preview protein repairs'

    def run(self,ctx):
        from mds.chem.repair import propose_repairs
        return propose_repairs(self.source,selected_residues=self.selected_residues,ph=self.ph,sequence_references=self.sequence_references,
                               progress=ctx.progress,check_cancelled=ctx.check_cancelled)


class ProteinRepairJob(Job):
    def __init__(self,plan,selected,output):
        self.plan,self.selected,self.output=plan,list(selected),Path(output)
        self.name='Build protein repair candidate'

    def run(self,ctx):
        from mds.chem.repair_process import run_repair_process
        return run_repair_process(self.plan,self.selected,self.output,ctx)
