"""Sequence-backed, explicitly selected short-loop models; never inferred from numbering alone."""
from __future__ import annotations
from collections import defaultdict
import io
import math
from pathlib import Path
import random
import gemmi
from mds.chem.health import SIDECHAINS, VDW
from mds.chem.edit import atom_records

MAX_LOOP = 8
SEEDS = (42, 1337, 2026)


def _chains(st):
    names=[c.name for c in st[0]]
    return {c.name:c for c in st[0]} if len(set(names))==len(names) else {}


def _positions(chain, sequence):
    """Accept authoritative label_seq, or a unique SEQRES/author-number alignment."""
    residues=[r for r in chain if r.name in SIDECHAINS]
    if not residues or any(r.seqid.icode.strip() for r in residues): return None
    labels=[r.label_seq for r in residues]
    if all(n is not None and 1<=n<=len(sequence) and sequence[n-1]==r.name for r,n in zip(residues,labels)):
        if labels==sorted(set(labels)): return {r.seqid.num:n for r,n in zip(residues,labels)}
    first=residues[0].seqid.num
    span=residues[-1].seqid.num-first+1
    if span<1 or span>len(sequence): return None
    if len({r.seqid.num for r in residues})!=len(residues): return None
    offsets=[offset for offset in range(len(sequence)-span+1) if all(
        0<=offset+r.seqid.num-first<len(sequence) and sequence[offset+r.seqid.num-first]==r.name for r in residues)]
    if len(offsets)!=1: return None
    return {r.seqid.num:offsets[0]+r.seqid.num-first+1 for r in residues}


def propose_loops(source, st, segment_starts, references=(), selected_residues=None, *, check_cancelled=lambda: None):
    from mds.chem.repair import RepairChange, fingerprint
    source=Path(source)
    candidates=[source]+[Path(p) for p in references if Path(p)!=source]
    chains=_chains(st)
    if not chains: return [], ['Loop modelling unavailable: duplicate author chain IDs are ambiguous.']
    evidence=defaultdict(list)
    notes=[]
    for path in candidates:
        check_cancelled()
        try:
            ref=gemmi.read_structure(str(path)); refs=_chains(ref)
            if not refs: continue
            for ci,ch in enumerate(st[0]):
                if ch.name not in refs: continue
                rc=refs[ch.name]
                ref.setup_entities()
                entity=ref.get_entity_of(rc.get_polymer())
                sequence=list(entity.full_sequence) if entity else []
                if not sequence or any(n not in SIDECHAINS for n in sequence): continue
                mapping=_positions(rc,sequence)
                if mapping is None: continue
                ref_res={r.seqid.num:r for r in rc if r.name in SIDECHAINS}
                observed=[(ri,r) for ri,r in enumerate(ch) if r.name in SIDECHAINS]
                if not observed: continue
                # Automatic recovery from project originals needs matching source anchors,
                # not just a similar sequence or a coincidentally identical file name.
                if any(r.seqid.icode.strip() or r.seqid.num not in mapping or r.seqid.num not in ref_res or
                       ref_res[r.seqid.num].name!=r.name or any(
                           r.find_atom(n,'*') is None or ref_res[r.seqid.num].find_atom(n,'*') is None or
                           r.find_atom(n,'*').pos.dist(ref_res[r.seqid.num].find_atom(n,'*').pos)>0.002
                           for n in ('N','CA','C')) for _,r in observed): continue
                for (li,left),(ri,right) in zip(observed,observed[1:]):
                    first,last=mapping[left.seqid.num]+1,mapping[right.seqid.num]-1
                    count=last-first+1
                    if count<=0: continue
                    label=f'{ch.name or "(blank)"}:{left.seqid}-{right.seqid}'
                    if count>MAX_LOOP:
                        notes.append(f'{label}: {count} missing residues exceed the local short-loop limit ({MAX_LOOP}); specialist modelling required.');continue
                    # Do not cross explicit termini, insertion codes, components or ambiguous numbering.
                    if ri!=li+1 or (ci,ri) in segment_starts or (left.subchain and right.subchain and left.subchain!=right.subchain): continue
                    if right.seqid.num-left.seqid.num!=count+1: continue
                    if any(a.altloc not in ('\0',' ') for r in (left,right) for a in r):
                        notes.append(f'{label}: resolve anchor alternate conformations before loop modelling.');continue
                    distance=left.find_atom('C','*').pos.dist(right.find_atom('N','*').pos)
                    if not math.isfinite(distance) or distance>3.8*(count+1):
                        notes.append(f'{label}: endpoints exceed the short-loop reach limit.');continue
                    spec=dict(chain_index=ci,left_index=li,right_index=ri,chain=ch.name,
                        left_number=left.seqid.num,right_number=right.seqid.num,
                        numbers=list(range(left.seqid.num+1,right.seqid.num)),
                        names=sequence[first-1:last],label_positions=list(range(first,last+1)),
                        sequence=sequence,sequence_positions={str(k):v for k,v in mapping.items()},
                        evidence=str(path),evidence_sha256=fingerprint(path))
                    if selected_residues is not None and not ({(0,ci,li),(0,ci,ri)} & set(selected_residues)): continue
                    evidence[(ci,li,ri)].append(spec)
        except (OSError,RuntimeError,ValueError): continue
    changes=[]
    for (ci,li,ri),specs in sorted(evidence.items()):
        signatures={(tuple(s['names']),tuple(s['numbers']),tuple(s['label_positions']),tuple(s['sequence'])) for s in specs}
        if len(signatures)!=1:
            notes.append(f'Chain {ci}, gap at residue index {ri}: conflicting sequence evidence; no loop proposed.');continue
        spec=specs[0]
        changes.append(RepairChange(f'loop:{ci}:{li}:{ri}','loop',ci,ri,
            f"{spec['chain']}:{spec['numbers'][0]}-{spec['numbers'][-1]}: model {len(spec['names'])} missing residues ({'-'.join(spec['names'])}); computational, review required",
            recommended=False,loop=spec))
    if not changes: notes.append('No eligible sequence-backed internal short loops. Numbering gaps alone never authorise reconstruction.')
    else: notes.append('Loop modelling is opt-in and separate from atom completion. PDBFixer builds missing residues; local Amber14 refinement keeps all existing coordinates fixed. No experimental conformation or docking accuracy is implied.')
    return changes, list(dict.fromkeys(notes))


def _snapshot(st, specs, segment_starts):
    """Reserve backend residue IDs for insertions; bridge only the approved gaps."""
    work=gemmi.Structure(); model=gemmi.Model('1'); mapping={}
    bridge={(s['chain_index'],s['right_index']) for s in specs}
    for ci,ch in enumerate(st[0]):
        segment=0; wc=gemmi.Chain(f'W{ci}_{segment}');previous=None
        for ri,res in enumerate(ch):
            if len(wc) and ((ci,ri) in segment_starts or ((ci,ri) not in bridge and (
                previous.find_atom('C','*') is None or res.find_atom('N','*') is None or
                previous.find_atom('C','*').pos.dist(res.find_atom('N','*').pos)>1.8 or
                (previous.subchain and res.subchain and previous.subchain!=res.subchain)))):
                model.add_chain(wc);segment+=1;wc=gemmi.Chain(f'W{ci}_{segment}')
            clone=res.clone();clone.seqid=gemmi.SeqId((ri+1)*100,' ');clone.subchain=wc.name
            wc.add_residue(clone);mapping[(wc.name,str(clone.seqid.num))]=(ci,ri);previous=res
        if len(wc): model.add_chain(wc)
    work.add_model(model);work.setup_entities()
    return work,mapping


def _model_initial(st, specs, segment_starts, seed):
    import openmm
    from openmm import unit
    from pdbfixer import PDBFixer
    from mds.chem.repair import _retain_disulfides
    platform=openmm.Platform.getPlatformByName('CPU');platform.setPropertyDefaultValue('Threads','2')
    work,mapping=_snapshot(st,specs,segment_starts)
    fixer=PDBFixer(pdbxfile=io.StringIO(work.make_mmcif_document().as_string()),platform=platform)
    _retain_disulfides(fixer,mapping,st)
    fixer.missingResidues={};new_mapping={}
    for chain in fixer.topology.chains():
        for index,res in enumerate(chain.residues()):
            key=mapping[(chain.id,res.id)]
            for si,s in enumerate(specs):
                if key==(s['chain_index'],s['right_index']):
                    fixer.missingResidues[(chain.index,index)]=s['names']
                    for j in range(len(s['names'])):
                        # PDBFixer wraps newly generated author IDs modulo 10000.
                        new_key=(chain.id,str((int(res.id)-len(s['names'])+j)%10000))
                        if new_key in new_mapping or new_key in mapping:
                            raise ValueError('Ambiguous backend insertion IDs; model smaller chain segments.')
                        new_mapping[new_key]=(si,j)
    fixer.findMissingAtoms();fixer.missingAtoms={};fixer.missingTerminals={}
    fixer.addMissingAtoms(seed=seed)
    loops=defaultdict(list)
    for atom,pos in zip(fixer.topology.atoms(),fixer.positions.value_in_unit(unit.angstrom)):
        key=(atom.residue.chain.id,atom.residue.id)
        if key in new_mapping:
            si,j=new_mapping[key];loops[(si,j)].append((atom.name,atom.element.symbol,tuple(pos)))
        elif key not in mapping or st[0][mapping[key][0]][mapping[key][1]].find_atom(atom.name,'*') is None:
            raise ValueError(f'Loop backend attempted an unapproved non-loop addition: {key}, {atom.name}.')
    return loops


def _insert(st, specs, loops):
    serial=max([a.serial for m in st for ch in m for r in ch for a in r]+[n for k,v in st.conect_map.items() for n in (k,*v)]+[0])
    for si,s in sorted(enumerate(specs),key=lambda pair:(pair[1]['chain_index'],pair[1]['right_index']),reverse=True):
        ch=st[0][s['chain_index']]
        for j in reversed(range(len(s['names']))):
            r=gemmi.Residue();r.name=s['names'][j];r.seqid=gemmi.SeqId(s['numbers'][j],' ')
            r.het_flag='A';r.subchain=ch[s['left_index']].subchain;r.entity_id=ch[s['left_index']].entity_id;r.label_seq=s['label_positions'][j]
            expected=set(('N CA C O '+SIDECHAINS[r.name]).split())
            values=loops.get((si,j),[])
            if {a[0] for a in values}!=expected: raise ValueError('Backend returned an incomplete or unexpected loop residue.')
            for name,element,xyz in values:
                serial+=1;a=gemmi.Atom();a.name=name;a.element=gemmi.Element(element);a.pos=gemmi.Position(*xyz)
                a.serial=serial;a.occ=1;a.b_iso=0;r.add_atom(a)
            ch.add_residue(r,pos=s['right_index'])
    return st


def _stereocentres(residue):
    if residue.name=='GLY': return []
    centres=[('N','CA','C','CB')]
    if residue.name=='THR': centres.append(('CA','CB','OG1','CG2'))
    if residue.name=='ILE': centres.append(('CA','CB','CG1','CG2'))
    return centres


def _initialise_chirality(residue):
    """Correct generated stereocentres before refinement; never touch input atoms."""
    import numpy as np
    for index,names in enumerate(_stereocentres(residue)):
        atoms=[residue.find_atom(n,'*') for n in names]
        p=[np.array(tuple(a.pos)) for a in atoms]
        normal=np.cross(p[0]-p[1],p[2]-p[1]);length=np.linalg.norm(normal)
        if length<1e-8: raise ValueError('Degenerate generated loop stereocentre.')
        normal/=length
        if np.dot(normal,p[3]-p[1])>0: continue
        target=p[3]-2*np.dot(normal,p[3]-p[1])*normal
        if index:  # CG2 is terminal; preserve the alpha centre and the other branch.
            atoms[3].pos=gemmi.Position(*target);continue
        old=p[3]-p[1];new=target-p[1]
        old/=np.linalg.norm(old);new/=np.linalg.norm(new)
        axis=np.cross(old,new);sine=np.linalg.norm(axis);cosine=np.dot(old,new)
        if sine<1e-8: raise ValueError('Ambiguous generated side-chain orientation.')
        x,y,z=axis
        cross=np.array([[0,-z,y],[z,0,-x],[-y,x,0]])
        rotation=np.eye(3)+cross+cross@cross*((1-cosine)/(sine*sine))
        for a in residue:
            if a.name not in ('N','CA','C','O','OXT'):
                a.pos=gemmi.Position(*(p[1]+rotation@(np.array(tuple(a.pos))-p[1])))


def _refine(st, spec, seed, check_cancelled, ph=7.4):
    """Amber14 patch dynamics/minimisation; only loop atoms can move.

    Temporary completion/protonation of the two flanking residues is used solely
    to parameterise the isolated patch and is never merged into the output.
    """
    import openmm
    from openmm import unit
    from openmm.app import ForceField,CutoffNonPeriodic
    from pdbfixer import PDBFixer
    from mds.chem.repair import _hydrogen_lock
    ch=st[0][spec['chain_index']]
    ri=next(i for i,r in enumerate(ch) if r.seqid.num==spec['left_number'] and not r.seqid.icode.strip())
    n=len(spec['names']); actual=list(ch)[ri:ri+n+2]
    for r in actual[1:-1]: _initialise_chirality(r)
    patch=gemmi.Structure();m=gemmi.Model('1');c=gemmi.Chain('P')
    for i,r in enumerate(actual):
        clone=r.clone();clone.seqid=gemmi.SeqId(i+1,' ');clone.subchain='P'
        for ai in reversed(range(len(clone))):
            if clone[ai].is_hydrogen(): del clone[ai]
        c.add_residue(clone)
    m.add_chain(c);patch.add_model(m);patch.setup_entities()
    platform=openmm.Platform.getPlatformByName('CPU');platform.setPropertyDefaultValue('Threads','2')
    f=PDBFixer(pdbxfile=io.StringIO(patch.make_mmcif_document().as_string()),platform=platform)
    f.missingResidues={};f.findMissingAtoms();f.addMissingAtoms(seed=seed)
    ff=ForceField('amber14/protein.ff14SB.xml','amber14/tip3p.xml')
    with _hydrogen_lock:
        state=random.getstate()
        try: random.seed(seed);f.addMissingHydrogens(ph,forcefield=ff)
        finally: random.setstate(state)
    system=ff.createSystem(f.topology,nonbondedMethod=CutoffNonPeriodic,nonbondedCutoff=0.5*unit.nanometer,constraints=None,removeCMMotion=False)
    atoms=list(f.topology.atoms()); movable={a.index for a in atoms if 1<int(a.residue.id)<n+2}
    heavy={a.index for a in atoms if a.index in movable and a.element.symbol!='H'}
    for a in atoms:
        if a.index not in movable: system.setParticleMass(a.index,0)
    # Force-field torsions alone can leave a reflected initial stereocentre in
    # its local minimum. Explicit L-configuration restraints prevent inversion.
    chirality=openmm.CustomCompoundBondForce(4,
        'chirality_k*min(0,volume-0.0015)^2; volume='
        '(x1-x2)*((y3-y2)*(z4-z2)-(z3-z2)*(y4-y2))'
        '-(y1-y2)*((x3-x2)*(z4-z2)-(z3-z2)*(x4-x2))'
        '+(z1-z2)*((x3-x2)*(y4-y2)-(y3-y2)*(x4-x2))')
    chirality.addGlobalParameter('chirality_k',1e9)
    for residue in f.topology.residues():
        if not 1<int(residue.id)<n+2: continue
        lookup={a.name:a.index for a in residue.atoms()}
        for names in _stereocentres(residue): chirality.addBond([lookup[name] for name in names],[])
    system.addForce(chirality)
    # Explicit trans-peptide restraints prevent generic-template cis/nonplanar
    # minima. The independent output validator still checks every junction.
    peptide=openmm.CustomTorsionForce('omega_k*(1+cos(theta))')
    peptide.addGlobalParameter('omega_k',1000)
    residues=list(f.topology.residues())
    for left,right in zip(residues,residues[1:]):
        a={a.name:a.index for a in left.atoms()};b={a.name:a.index for a in right.atoms()}
        peptide.addTorsion(a['CA'],a['C'],b['N'],b['CA'],[])
    system.addForce(peptide)
    backbone_angles=openmm.HarmonicAngleForce()
    for left,right in zip(residues,residues[1:]):
        a={a.name:a.index for a in left.atoms()};b={a.name:a.index for a in right.atoms()}
        backbone_angles.addAngle(a['CA'],a['C'],b['N'],math.radians(116.6),5000)
        backbone_angles.addAngle(a['C'],b['N'],b['CA'],math.radians(121.7),5000)
    for residue in residues:
        if not 1<int(residue.id)<n+2: continue
        a={a.name:a.index for a in residue.atoms()}
        backbone_angles.addAngle(a['N'],a['CA'],a['C'],math.radians(111.2),5000)
        backbone_angles.addAngle(a['CA'],a['C'],a['O'],math.radians(120.8),5000)
    system.addForce(backbone_angles)
    xyz=list(f.positions.value_in_unit(unit.nanometer))
    original_loop=[a.pos for r in actual[1:-1] for a in r]
    # Fixed environmental atoms retain waters, cofactors and ligands around the loop.
    env=[]
    addresses={(ch.name,r.seqid.num,r.seqid.icode) for r in actual}
    for source_ch in st[0]:
        for r in source_ch:
            if (source_ch.name,r.seqid.num,r.seqid.icode) in addresses: continue
            for a in r:
                if a.occ>0 and not a.is_hydrogen() and a.element.vdw_r>0 and any(a.pos.dist(p)<8.0 for p in original_loop): env.append(a)
    repulsion=openmm.CustomNonbondedForce('0.5*k*step(radius1+radius2-buffer-r)*(radius1+radius2-buffer-r)^2')
    repulsion.addGlobalParameter('k',10000);repulsion.addGlobalParameter('buffer',0.04)
    repulsion.addPerParticleParameter('radius')
    for a in atoms: repulsion.addParticle([VDW.get(a.element.symbol,1.1)*0.1])
    nb=next(force for force in system.getForces() if isinstance(force,openmm.NonbondedForce))
    env_indices=[]
    for a in env:
        index=system.addParticle(0);env_indices.append(index)
        nb.addParticle(0,0.1,0);repulsion.addParticle([VDW.get(a.element.name,a.element.vdw_r)*0.1]);xyz.append(openmm.Vec3(a.pos.x*0.1,a.pos.y*0.1,a.pos.z*0.1))
    if env_indices:
        # OpenMM's CPU kernels share neighbour lists and require matching exclusions.
        for i in range(nb.getNumExceptions()):
            a,b,*_=nb.getExceptionParameters(i);repulsion.addExclusion(a,b)
        repulsion.addInteractionGroup(heavy,set(env_indices));repulsion.setNonbondedMethod(openmm.CustomNonbondedForce.CutoffNonPeriodic);repulsion.setCutoffDistance(0.5)
        system.addForce(repulsion)
    integrator=openmm.LangevinMiddleIntegrator(300*unit.kelvin,5/unit.picosecond,0.001*unit.picoseconds)
    integrator.setRandomNumberSeed(seed)
    context=openmm.Context(system,integrator,platform);context.setPositions(xyz*unit.nanometer)
    openmm.LocalEnergyMinimizer.minimize(context,maxIterations=1500)
    context.setVelocitiesToTemperature(300*unit.kelvin,seed)
    # Several local conformations rather than a single template placement.
    for temperature in (600,400,300):
        check_cancelled();integrator.setTemperature(temperature*unit.kelvin);integrator.step(1500)
    openmm.LocalEnergyMinimizer.minimize(context,tolerance=5,maxIterations=3000)
    state=context.getState(getPositions=True,getEnergy=True)
    result=state.getPositions(asNumpy=True).value_in_unit(unit.angstrom)
    energy=state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    for a in atoms:
        if a.index in heavy:
            target=actual[int(a.residue.id)-1].find_atom(a.name,'*')
            if target is None: raise ValueError('Unexpected loop refinement atom.')
            target.pos=gemmi.Position(*result[a.index])
    del context,integrator
    return float(energy)


def _topology(st, specs):
    from mds.chem.health import SIDE_BONDS
    bonds=[];added=set();regions=[]
    for s in specs:
        ch=st[0][s['chain_index']];ci=s['chain_index']
        start=next(i for i,r in enumerate(ch) if r.seqid.num==s['left_number'] and not r.seqid.icode.strip())
        stop=start+len(s['names'])+1
        for ri,r in enumerate(ch):
            if r.name not in SIDECHAINS: continue
            names={a.name for a in r}
            pairs='N-CA CA-C C-O CA-CB C-OXT '+SIDE_BONDS[r.name]
            for pair in pairs.split():
                a,b=pair.split('-')
                if a in names and b in names: bonds.append(((ci,ri,a),(ci,ri,b)))
            if start<ri<stop: added.update((ci,ri,a.name) for a in r)
        for ri in range(start,stop): bonds.append(((ci,ri,'C'),(ci,ri+1,'N')))
        regions.append((ci,start,stop))
    # Include normal peptide bonds elsewhere for clash exclusions at the anchor boundaries.
    for ci,ch in enumerate(st[0]):
        for ri in range(len(ch)-1):
            a=ch[ri].find_atom('C','*');b=ch[ri+1].find_atom('N','*')
            if a and b and 1.1<=a.pos.dist(b.pos)<=1.8: bonds.append(((ci,ri,'C'),(ci,ri+1,'N')))
    return added,bonds,regions


def validate_loops(st,specs):
    from mds.chem.repair import _validate_additions
    added,bonds,regions=_topology(st,specs)
    radii={a.element.name:float(a.element.vdw_r) for m in st for ch in m for r in ch for a in r if a.element.vdw_r>0}
    radii.update(VDW);radii['H']=1.1
    blockers=_validate_additions(st,added,bonds,radii=radii);metrics=[]
    metal_atoms=[(ch.name,r.seqid,a) for ch in st[0] for r in ch for a in r if a.occ>0 and a.element.is_metal]
    for ci,ri,name in sorted(added):
        atom=st[0][ci][ri].find_atom(name,'*')
        for chain,number,metal in metal_atoms:
            if atom.pos.dist(metal.pos)<3.0:
                blockers.append(f'Loop approaches metal {metal.element.name} {chain}:{number}; specialist coordination review required.')
    for ci,start,stop in regions:
        ch=st[0][ci]
        for ri in range(start,stop):
            left,right=ch[ri],ch[ri+1]
            a=left.find_atom('C','*');b=right.find_atom('N','*')
            d=a.pos.dist(b.pos);omega=math.degrees(gemmi.calculate_dihedral(left.find_atom('CA','*').pos,a.pos,b.pos,right.find_atom('CA','*').pos))
            tag=f'{ch.name}:{left.seqid}-{right.seqid}'
            if not 1.20<=d<=1.45: blockers.append(f'Loop peptide closure outside 1.20-1.45 A: {tag}, {d:.3f} A.')
            if not math.isfinite(omega) or not (abs(omega)>=150 or (right.name=='PRO' and abs(omega)<=30)):
                blockers.append(f'Nonplanar loop peptide: {tag}, omega {omega:.1f} degrees.')
            junction_angles=[math.degrees(gemmi.calculate_angle(left.find_atom('CA','*').pos,a.pos,b.pos)),
                             math.degrees(gemmi.calculate_angle(a.pos,b.pos,right.find_atom('CA','*').pos))]
            if any(not math.isfinite(v) or not 110<=v<=135 for v in junction_angles):
                blockers.append(f'Loop junction bond angle abnormal: {tag}, {junction_angles[0]:.1f}/{junction_angles[1]:.1f} degrees.')
            metrics.append(dict(region=tag,peptide_CN_A=d,omega_degrees=omega,junction_angles_degrees=junction_angles))
        for ri in range(start+1,stop):
            r=ch[ri];n,ca,c=[r.find_atom(name,'*') for name in ('N','CA','C')]
            oxygen=r.find_atom('O','*')
            for aa,bb,low,high in ((n,ca,1.35,1.60),(ca,c,1.40,1.65),(c,oxygen,1.15,1.35)):
                length=aa.pos.dist(bb.pos)
                if not low<=length<=high:
                    blockers.append(f'Loop backbone bond abnormal: {ch.name}:{r.seqid} {aa.name}-{bb.name}, {length:.3f} A.')
            carbonyl_angle=math.degrees(gemmi.calculate_angle(ca.pos,c.pos,oxygen.pos))
            if not 110<=carbonyl_angle<=135: blockers.append(f'Loop carbonyl bond angle abnormal: {ch.name}:{r.seqid}, {carbonyl_angle:.1f} degrees.')
            angle=math.degrees(gemmi.calculate_angle(n.pos,ca.pos,c.pos))
            if not 100<=angle<=125: blockers.append(f'Loop N-CA-C angle abnormal: {ch.name}:{r.seqid}, {angle:.1f} degrees.')
            if r.name!='GLY':
                cb=r.find_atom('CB','*')
                volume=(n.pos-ca.pos).cross(c.pos-ca.pos).dot(cb.pos-ca.pos)
                if not math.isfinite(volume) or volume<=0.5:
                    blockers.append(f'Loop L-amino-acid chirality abnormal: {ch.name}:{r.seqid}.')
            for names in _stereocentres(r)[1:]:
                p=[r.find_atom(name,'*').pos for name in names]
                if (p[0]-p[1]).cross(p[2]-p[1]).dot(p[3]-p[1])<=0.5:
                    blockers.append(f'Loop side-chain stereochemistry abnormal: {ch.name}:{r.seqid}.')
            phi,psi=gemmi.calculate_phi_psi(ch[ri-1],r,ch[ri+1])
            metrics.append(dict(residue=f'{ch.name}:{r.seqid}',phi_degrees=math.degrees(phi),psi_degrees=math.degrees(psi)))
    return list(dict.fromkeys(blockers)),metrics


def build_loop_repairs(plan,changes,selected,output_dir,progress,check_cancelled):
    from uuid import uuid4
    import openmm
    from mds.chem.repair import (RepairResult,fingerprint,_write_structure,reconstructed_atom_indices)
    from mds.chem.health import check_protein
    specs=[c.loop for c in changes]
    if any(not s for s in specs): raise ValueError('Loop proposal has no sequence evidence.')
    source=Path(plan.source)
    for s in specs:
        if fingerprint(s['evidence'])!=s['evidence_sha256']: raise ValueError('Sequence evidence changed; generate a fresh proposal.')
    original=gemmi.read_structure(str(source)); before=atom_records(original)
    def address(a): return (a['model'],a['chain'],a['resi'],a['icode'],a['resn'],a['atom'],a['altloc'])
    original_addresses={address(a) for a in before}
    if len(original_addresses)!=len(before): raise ValueError('Ambiguous atom addresses; loop modelling unavailable.')
    inherited_addresses={address(before[i]) for i in reconstructed_atom_indices(source,before)}
    trials=[];best=None
    for attempt,seed in enumerate(SEEDS):
        check_cancelled();progress(5+attempt*27,f'Modelling short loops: candidate {attempt+1}/{len(SEEDS)}')
        try:
            loops=_model_initial(original,specs,plan.segment_starts,seed)
            st=_insert(original.clone(),specs,loops);energies=[]
            for s in specs:
                check_cancelled();energies.append(_refine(st,s,seed,check_cancelled,plan.ph))
            blockers,metrics=validate_loops(st,specs)
            if not all(math.isfinite(v) for v in energies): raise ValueError('Non-finite refinement energy.')
            trials.append(dict(seed=seed,blockers=blockers,patch_energy_kJ_mol=energies,geometry=metrics))
            score=(len(blockers),sum(energies))
            if best is None or score<best[0]: best=(score,st,seed)
        except (ValueError,RuntimeError,openmm.OpenMMException) as exc:
            trials.append(dict(seed=seed,error=str(exc)))
    if best is None: raise ValueError('No loop candidate could be modelled: '+'; '.join(t.get('error','') for t in trials))
    _,st,seed=best
    # Retain full sequence evidence on output, including when cleanup previously stripped it.
    st.setup_entities()
    for s in specs:
        ch=st[0][s['chain_index']];entity=st.get_entity_of(ch.get_polymer())
        if entity: entity.full_sequence=s['sequence']
        for r in ch:
            if str(r.seqid.num) in s['sequence_positions']: r.label_seq=s['sequence_positions'][str(r.seqid.num)]
    records=atom_records(st)
    added=[a for a in records if address(a) not in original_addresses]
    inherited=[a for a in records if address(a) in inherited_addresses]
    directory=Path(output_dir);directory.mkdir(parents=True,exist_ok=True)
    suffix='.cif' if source.suffix.lower() in ('.cif','.mmcif') else '.pdb'
    candidate=directory/f'protein_{uuid4().hex[:12]}_loops{suffix}'
    _write_structure(st,candidate,added+inherited,plan.source_sha256,source)
    saved=gemmi.read_structure(str(candidate));saved_records=atom_records(saved)
    surviving=[a for a in saved_records if address(a) in original_addresses]
    if len(surviving)!=len(before) or any(address(a)!=address(b) or any(abs(a[k]-b[k])>1e-5 for k in ('x','y','z')) for a,b in zip(before,surviving)):
        raise ValueError('Loop candidate changed existing coordinates/identities; rejected.')
    blockers,metrics=validate_loops(saved,specs)
    after=check_protein(candidate,progress=lambda p,m:progress(87+p*.12,m),check_cancelled=check_cancelled)
    added_addresses={address(a) for a in added}
    for issue in after.issues:
        if issue.code in ('bond_geometry','angle_geometry','invalid_coordinates','backbone_connectivity') and any(address(a) in added_addresses for a in issue.atoms):
            blockers.append('Post-write loop geometry finding: '+issue.location+' / '+issue.code)
    if fingerprint(source)!=plan.source_sha256 or any(fingerprint(s['evidence'])!=s['evidence_sha256'] for s in specs):
        raise ValueError('Source or sequence evidence changed during modelling; candidate rejected.')
    result=RepairResult(plan,list(selected),str(candidate),after,added,[],blockers,[
        'Computational short-loop hypotheses, not recovered experimental coordinates. Review before use; the native conformation is unknown.',
        'No statistical Ramachandran/rotamer, pKa, charge, metal-coordination or docking-readiness validation. Phi/psi values are reported for manual review.',
        'Peptides are modelled with trans-omega restraints, including proline; cis-proline alternatives are not sampled.',
        'Only loop heavy atoms are retained from refinement. Temporary patch hydrogens/terminal atoms/anchor completion are discarded; prepare protonation and charges separately.',
    ],fingerprint(candidate),inherited,[],loops=specs,modelling=dict(backend='PDBFixer + local Amber14/OpenMM',seeds=list(SEEDS),selected_seed=seed,trials=trials,geometry=metrics,
        forcefield='amber14/protein.ff14SB.xml + amber14/tip3p.xml',temporary_patch_ph=plan.ph,temperature_schedule_K=[600,400,300],steps_per_temperature=1500,time_step_ps=0.001,
        environment_radius_A=8,environment_overlap_buffer_A=0.4,environment_repulsion_k_kJ_mol_nm2=10000,L_chirality_restraint_min_A3=1.5,chirality_gate_min_A3=0.5,chirality_restraint_k=1e9,peptide_omega_assumption="trans",peptide_omega_restraint_k_kJ_mol=1000,backbone_angle_restraint_k_kJ_mol_rad2=5000,backbone_angle_targets_degrees=[116.6,121.7,111.2,120.8],existing_coordinates_fixed=True))
    result.export_report(candidate.with_suffix(candidate.suffix+'.repair.json'));progress(100,'Loop candidate ready for review' if result.can_accept else 'Loop candidates blocked by geometry checks')
    return result

