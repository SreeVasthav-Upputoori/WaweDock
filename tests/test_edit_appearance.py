"""Real atom deletion/replacement must retain view preferences and camera."""
from pathlib import Path
import json
import gemmi
import pytest
from conftest import wait_until
from mds.chem.edit import atom_records,edit_atoms
from mds.render.mol3d import Mol3DViewer

DATA=Path(__file__).parent/'data'


@pytest.fixture
def viewer(qapp):
    widget=Mol3DViewer(); widget.resize(600,450); widget.show()
    assert wait_until(lambda:widget._ready,20)
    yield widget
    widget.close()


def js(viewer,expression):
    result=[]; viewer.page().runJavaScript(expression,result.append)
    assert wait_until(lambda:bool(result),10)
    return result[0]


def load(viewer,name,path):
    structure=gemmi.read_structure(str(path))
    viewer.add_object(name,structure.make_pdb_string(),'pdb','protein',atom_records(structure))
    return structure


def state(viewer,name):
    return json.loads(js(viewer,"""JSON.stringify((() => {
      const o=objects[NAME];
      return {rep:o.rep,color:o.color,customColor:o.customColor||null,
        components:o.components,ligandColor:o.ligandColor||null,labels:o.labels,spectrumRange:o.spectrumRange,
        background:backgroundColor,camera:viewer.getView()};
    })())""".replace('NAME',json.dumps(name))))


def test_new_view_defaults_black_and_spectrum(viewer):
    assert js(viewer,'backgroundColor')=='#000000'
    assert viewer.page().backgroundColor().name()=='#000000'
    load(viewer,'protein',DATA/'3PTB.pdb')
    assert js(viewer,"objects.protein.color")=='spectrum'
    assert js(viewer,"objects.protein.model.selectedAtoms({atom:'CA',hetflag:false})[0].style.cartoon.colorscheme.gradient.gradient")=='ROYGB'
    assert js(viewer,"objects.protein.model.selectedAtoms({atom:'CA',hetflag:false})[0].style.cartoon.colorscheme.gradient.valueToHex(objects.protein.spectrumRange[0])")==0x0000ff
    assert js(viewer,"objects.protein.model.selectedAtoms({atom:'CA',hetflag:false})[0].style.cartoon.colorscheme.gradient.valueToHex(objects.protein.spectrumRange[1])")==0xff0000


@pytest.mark.parametrize('scheme,custom',[('spectrum',None),('residue',None),('chain','#bada55')])
def test_hetatom_deletion_undo_redo_keep_appearance(viewer,tmp_path,scheme,custom):
    source=DATA/'3PTB.pdb'
    original=load(viewer,'source',source)
    records=atom_records(original)
    selected=[r for r in records if r['resn'] in {'BEN','HOH','CA'}]
    edited=tmp_path/'edited.pdb'
    removed=edit_atoms(source,edited,selected)
    assert removed>0
    edited_structure=gemmi.read_structure(str(edited))
    js(viewer,"setBackgroundColor('#203040'); setRepresentation('source','ribbon'); "
       +f"setColorScheme('source',{json.dumps(scheme)}); "
       +(f"setCustomColor('source',{json.dumps(custom)}); " if custom else '')
       +"setLigandColor('source','#00bfff',true); setComponentVisible('source','waters',false); true")
    before=state(viewer,'source')
    js(viewer,"pickMolecularSelection('source',objects.source.model.selectedAtoms({atom:'CA',hetflag:false})[0],{button:0},true); true")
    def replace(old,name,structure):
        viewer.replace_structure(old,name,structure.make_pdb_string(),'pdb',atom_records(structure))
        assert state(viewer,name)==before
        assert js(viewer,'Object.keys(objects).length')==1
        assert js(viewer,f"objects[{json.dumps(name)}].model.selectedAtoms({{}}).length")==len(atom_records(structure))
        assert js(viewer,f"objects[{json.dumps(name)}].model.selectedAtoms({{atom:'CA',hetflag:false}}).every(a=>a.style.cartoon.color!=='#ffff00')")
    replace('source','edited',edited_structure)
    replace('edited','source',original)
    replace('source','edited',edited_structure)
    assert viewer.page().backgroundColor().name()=='#203040'


def test_visible_prepared_protein_style_carries_to_editable_source(viewer,tmp_path):
    structure=load(viewer,'source',DATA/'3PTB.pdb')
    load(viewer,'prepared',DATA/'3PTB.pdb')
    js(viewer,"setBackgroundColor('#ffffff'); setRepresentation('prepared','sticks'); setCustomColor('prepared','#12ab34'); true")
    before=state(viewer,'prepared')
    viewer.copy_appearance('prepared','source')
    viewer.remove_object('prepared')
    viewer.replace_structure('source','edited',structure.make_pdb_string(),'pdb',atom_records(structure))
    assert state(viewer,'edited')==before


def test_deleting_terminal_residue_does_not_rescale_rainbow(viewer,tmp_path):
    source=DATA/'3PTB.pdb'
    structure=load(viewer,'source',source)
    records=atom_records(structure)
    protein=[r for r in records if r['resn'] not in {'BEN','HOH','CA'}]
    last=max(r['resi'] for r in protein)
    edited=tmp_path/'shortened.pdb'
    edit_atoms(source,edited,[r for r in protein if r['resi']==last])
    color=js(viewer,"objects.source.model.selectedAtoms({atom:'CA',hetflag:false})[0].style.cartoon.colorscheme.gradient.valueToHex(100)")
    before=state(viewer,'source')
    changed=gemmi.read_structure(str(edited))
    viewer.replace_structure('source','edited',changed.make_pdb_string(),'pdb',atom_records(changed))
    assert state(viewer,'edited')==before
    assert js(viewer,"objects.edited.model.selectedAtoms({atom:'CA',hetflag:false})[0].style.cartoon.colorscheme.gradient.valueToHex(100)")==color
