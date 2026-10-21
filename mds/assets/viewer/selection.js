"use strict";
// Selections contain canonical atom indices, never PDB serial numbers. The Python
// bridge supplies identities from the same Gemmi structure used for file edits.
const selection = {};
// Read-only overlays can be highlighted, but never sent to Python as editable atoms.
const visualSelection = {};
const container = document.getElementById("v");
let selectionMode = "residue", selectionShapes = [], lasso = null;
let suppressPickUntil = 0;

function bindAtomIdentities(name, identities) {
  const o = objects[name], atoms = o.model.selectedAtoms({});
  atoms.forEach(a => { a._wdObject = name; });
  if (!identities) return; // Prepared receptors and docked poses are read-only overlays.
  const valid = identities.length === atoms.length && atoms.every((a, i) => {
    const r = identities[i];
    return a.atom === r.atom && Math.abs(a.x-r.x) < 0.002 &&
      Math.abs(a.y-r.y) < 0.002 && Math.abs(a.z-r.z) < 0.002;
  });
  if (!valid) {
    emit("jserror", "Atom mapping failed for " + name + "; editing is disabled to protect the structure.");
    return;
  }
  atoms.forEach((a, i) => { a._wdIdentity = identities[i]; });
  o.editable = true;
}
function editableAtoms(name) {
  const o = objects[name];
  return o && o.editable && o.visible ? o.model.selectedAtoms({}).filter(a => componentVisible(o, a)) : [];
}
function residueKey(a) {
  const r = a._wdIdentity;
  return JSON.stringify([r.model, r.chain_index, r.residue_index]);
}
function groupKey(atom) {
  const r = atom._wdIdentity;
  if (selectionMode === "atom") return r.atom_index;
  if (selectionMode === "chain" || (selectionMode === "component" && !atom.hetflag)) {
    return JSON.stringify([r.model, r.chain]);
  }
  return residueKey(atom);
}
function atomGroups(atoms) {
  const groups = new Map();
  atoms.forEach(a => {
    const key = groupKey(a);
    if (!groups.has(key)) groups.set(key, []);
    groups.get(key).push(a);
  });
  return groups;
}
function groupAtoms(name, atom) {
  const key = groupKey(atom);
  return editableAtoms(name).filter(a => groupKey(a) === key);
}
function setSelectionMode(mode) {
  if (["atom", "residue", "chain", "component"].includes(mode)) selectionMode = mode;
}
function clearSelectionHighlights() {
  selectionShapes.forEach(s => viewer.removeShape(s));
  selectionShapes = [];
}
function reapplySelectionHighlights() {
  clearSelectionHighlights();
  const names = new Set([...Object.keys(selection), ...Object.keys(visualSelection)]);
  names.forEach(name => {
    const o = objects[name]; if (!o || !o.visible) return;
    const selected = o.editable ? selection[name] : visualSelection[name];
    if (!selected || !selected.size) return;
    const shape = viewer.addShape({color: "#ffff00", opacity: 0.65});
    selectionShapes.push(shape);
    o.model.selectedAtoms({}).forEach(a => {
      const id = o.editable ? a._wdIdentity.atom_index : a.index;
      if (componentVisible(o, a) && selected.has(id)) {
        // Extend beyond VDW spheres so the yellow selection stays visible.
        shape.addSphere({center: a, radius: a.elem === "H" ? 1.35 : 1.95});
      }
    });
  });
  viewer.render();
}
function visualSelectionCount() {
  let count = 0;
  Object.keys(visualSelection).forEach(name => {
    const o = objects[name]; if (!o || !o.visible) return;
    count += o.model.selectedAtoms({}).filter(a => componentVisible(o, a) && visualSelection[name].has(a.index)).length;
  });
  return count;
}
function emitSelection() {
  const detail = {}; let count = 0;
  Object.keys(selection).forEach(name => {
    const chosen = editableAtoms(name).filter(a => selection[name].has(a._wdIdentity.atom_index));
    if (chosen.length) detail[name] = chosen.map(a => a._wdIdentity);
    count += chosen.length;
  });
  const visualCount = visualSelectionCount();
  const total = count + visualCount;
  emit("selection", {detail: detail, summary: total ? total + " atom(s) selected" +
    (visualCount ? " (" + visualCount + " read-only)" : "") : ""});
}
function commitSelection() { reapplySelectionHighlights(); emitSelection(); }
function clearSelection() {
  Object.keys(visualSelection).forEach(name => delete visualSelection[name]);
  Object.keys(selection).forEach(name => delete selection[name]);
  commitSelection();
}
function pickSelection(name, atom, event) {
  if (Date.now() < suppressPickUntil || lasso) return;
  if (!objects[name].editable) {
    emit("measure", "This is a read-only overlay. Select the imported structure to edit it.");
    return;
  }
  const e = event && (event.originalEvent || event) || {};
  if (e.altKey || e.button === 2) return;
  const group = groupAtoms(name, atom), current = selection[name] || new Set();
  const remove = group.every(a => current.has(a._wdIdentity.atom_index));
  if (!e.shiftKey && !e.ctrlKey && !e.metaKey) {
    Object.keys(selection).forEach(n => delete selection[n]);
    Object.keys(visualSelection).forEach(n => delete visualSelection[n]);
    selection[name] = new Set();
  } else selection[name] = current;
  group.forEach(a => {
    const id = a._wdIdentity.atom_index;
    if ((e.shiftKey || e.ctrlKey || e.metaKey) && remove) selection[name].delete(id);
    else selection[name].add(id);
  });
  commitSelection();
}
function molecularResidueKey(a) {
  return a._wdIdentity ? residueKey(a) : JSON.stringify([a.chain, a.resi, a.icode || "", a.resn]);
}
function pickMolecularSelection(name, atom, event, whole = false) {
  const e = event && (event.originalEvent || event) || {};
  if (lasso || Date.now() < suppressPickUntil || e.altKey || e.button === 2 ||
      (!whole && e.detail >= 2)) return;
  const o = objects[name]; if (!o || !o.visible) return;
  const atoms = o.model.selectedAtoms({}).filter(a => componentVisible(o, a));
  const water = WATERS.includes(atom.resn);
  const proteinResidues = new Set(), backbone = new Map();
  if (whole && !water) atoms.forEach(a => {
    const key = molecularResidueKey(a);
    if (!a.hetflag) proteinResidues.add(key);
    if (!backbone.has(key)) backbone.set(key, new Set());
    if (a.atom !== "CA" || a.elem === "C") backbone.get(key).add(a.atom);
  });
  // Modified amino acids may use HETATM; identify their peptide backbone too.
  backbone.forEach((names, key) => {
    if (["N", "CA", "C"].every(name => names.has(name))) proteinResidues.add(key);
  });
  const group = atoms.filter(a => {
    if (whole && water) return WATERS.includes(a.resn);
    if (whole && proteinResidues.has(molecularResidueKey(atom)))
      return proteinResidues.has(molecularResidueKey(a));
    return molecularResidueKey(a) === molecularResidueKey(atom);
  });
  const target = o.editable ? selection : visualSelection;
  const atomId = a => o.editable ? a._wdIdentity.atom_index : a.index;
  const additive = e.shiftKey || e.ctrlKey || e.metaKey;
  const previous = target[name] || new Set();
  const remove = additive && group.every(a => previous.has(atomId(a)));
  if (!additive) {
    Object.keys(selection).forEach(n => delete selection[n]);
    Object.keys(visualSelection).forEach(n => delete visualSelection[n]);
    target[name] = new Set();
  } else target[name] = previous;
  group.forEach(a => {
    if (remove) target[name].delete(atomId(a)); else target[name].add(atomId(a));
  });
  commitSelection();
}
function expandSelectionToChain() {
  Object.keys(selection).forEach(name => {
    const atoms = editableAtoms(name), chains = new Set();
    atoms.forEach(a => { if (selection[name].has(a._wdIdentity.atom_index))
      chains.add(JSON.stringify([a._wdIdentity.model, a._wdIdentity.chain])); });
    atoms.forEach(a => { if (chains.has(JSON.stringify([a._wdIdentity.model, a._wdIdentity.chain])))
      selection[name].add(a._wdIdentity.atom_index); });
  });
  commitSelection();
}
function invertSelection() {
  Object.keys(objects).forEach(name => {
    const atoms = editableAtoms(name); if (!atoms.length) return;
    const previous = selection[name] || new Set();
    selection[name] = new Set(atoms.filter(a => !previous.has(a._wdIdentity.atom_index))
      .map(a => a._wdIdentity.atom_index));
  });
  commitSelection();
}
function setSelectionFromPython(detail) {
  Object.keys(visualSelection).forEach(name => delete visualSelection[name]);
  Object.keys(selection).forEach(name => delete selection[name]);
  Object.keys(detail).forEach(name => {
    selection[name] = new Set();
    editableAtoms(name).forEach(a => {
      if (detail[name].some(r => r.atom_index !== undefined ? r.atom_index === a._wdIdentity.atom_index :
          r.chain === a.chain && r.resi === a.resi &&
          (r.icode === undefined || r.icode === a._wdIdentity.icode)))
        selection[name].add(a._wdIdentity.atom_index);
    });
  });
  commitSelection();
}
function selectWithinRadius(radius) {
  if (!(radius > 0 && Number.isFinite(radius))) return;
  // A union of 3D spheres around the selected atoms, scoped to each source model.
  Object.keys(selection).forEach(name => {
    const atoms = editableAtoms(name), seeds = atoms.filter(a => selection[name].has(a._wdIdentity.atom_index));
    const cutoff2 = radius * radius;
    const grid = new Map(), expanded = new Set(), groups = atomGroups(atoms);
    const cell = a => [Math.floor(a.x/radius), Math.floor(a.y/radius), Math.floor(a.z/radius)];
    const gridKey = (model, x, y, z) => model+":"+x+":"+y+":"+z;
    seeds.forEach(a => {
      const p = cell(a), key = gridKey(a._wdIdentity.model, ...p);
      if (!grid.has(key)) grid.set(key, []);
      grid.get(key).push(a);
    });
    atoms.forEach(a => {
      const key = groupKey(a); if (expanded.has(key)) return;
      const p = cell(a); let hit = false;
      for (let x=-1; x<=1 && !hit; x++) for (let y=-1; y<=1 && !hit; y++) for (let z=-1; z<=1 && !hit; z++) {
        const near = grid.get(gridKey(a._wdIdentity.model, p[0]+x,p[1]+y,p[2]+z)) || [];
        hit = near.some(b => (a.x-b.x)**2 + (a.y-b.y)**2 + (a.z-b.z)**2 <= cutoff2);
      }
      if (!hit) return;
      expanded.add(key);
      groups.get(key).forEach(hit => selection[name].add(hit._wdIdentity.atom_index));
    });
  });
  commitSelection();
}

// Freehand screen polygon, selecting through the visible structure's depth.
// modelToScreen in the bundled 3Dmol returns PAGE coordinates (including canvas offset).
function pointInPolygon(point, polygon) {
  let inside = false;
  for (let i=0, j=polygon.length-1; i<polygon.length; j=i++) {
    const a=polygon[i], b=polygon[j];
    if ((a.y > point.y) !== (b.y > point.y) &&
        point.x < (b.x-a.x)*(point.y-a.y)/(b.y-a.y)+a.x) inside = !inside;
  }
  return inside;
}
function finishLasso(cancel) {
  if (!lasso) return;
  const active = lasso; lasso = null;
  active.svg.remove(); suppressPickUntil = Date.now()+200;
  if (cancel || active.points.length < 3) return;
  if (!active.additive) {
    Object.keys(selection).forEach(name => delete selection[name]);
    Object.keys(visualSelection).forEach(name => delete visualSelection[name]);
  }
  Object.keys(objects).forEach(name => {
    const atoms = editableAtoms(name), hits = new Set(), groups = atomGroups(atoms);
    atoms.forEach(a => {
      const p = viewer.modelToScreen(a);
      if (p && pointInPolygon(p, active.points)) {
        // Group only once per residue/chain even when many atoms fall in the polygon.
        const key = groupKey(a);
        if (hits.has(key)) return;
        hits.add(key);
        selection[name] = selection[name] || new Set();
        groups.get(key).forEach(hit => selection[name].add(hit._wdIdentity.atom_index));
      }
    });
  });
  commitSelection();
}
container.addEventListener("mousedown", e => {
  if (!e.altKey || e.button !== 0) return;
  e.preventDefault(); e.stopImmediatePropagation();
  const svg = document.createElementNS("http://www.w3.org/2000/svg", "svg");
  svg.style.cssText = "position:fixed;inset:0;width:100%;height:100%;pointer-events:none;z-index:50";
  const polygon = document.createElementNS(svg.namespaceURI, "polygon");
  polygon.setAttribute("fill", "rgba(14,165,233,0.18)");
  polygon.setAttribute("stroke", "#0284c7"); polygon.setAttribute("stroke-width", "2");
  svg.appendChild(polygon); document.body.appendChild(svg);
  lasso = {svg, polygon, points: [{x:e.pageX,y:e.pageY}], additive:e.shiftKey || e.ctrlKey || e.metaKey};
}, true);
document.addEventListener("mousemove", e => {
  if (!lasso) return;
  e.preventDefault(); e.stopImmediatePropagation();
  const last = lasso.points[lasso.points.length-1];
  if (Math.hypot(e.pageX-last.x,e.pageY-last.y) < 2) return;
  lasso.points.push({x:e.pageX,y:e.pageY});
  lasso.polygon.setAttribute("points", lasso.points.map(p =>
    (p.x-window.scrollX)+","+(p.y-window.scrollY)).join(" "));
}, true);
document.addEventListener("mouseup", e => {
  if (!lasso) return;
  e.preventDefault(); e.stopImmediatePropagation(); finishLasso(false);
}, true);
window.addEventListener("blur", () => finishLasso(true));
document.addEventListener("keydown", e => {
  if (e.key === "Escape") { finishLasso(true); clearSelection(); }
  if (e.key === "Delete" || e.key === "Backspace") {
    e.preventDefault(); emit("edit", "delete");
  }
  if ((e.ctrlKey || e.metaKey) && e.key.toLowerCase() === "z") {
    e.preventDefault(); emit("edit", e.shiftKey ? "redo" : "undo");
  }
  if ((e.ctrlKey || e.metaKey) && e.key.toLowerCase() === "y") {
    e.preventDefault(); emit("edit", "redo");
  }
});
// Raycast against the rendered atom geometry rather than choosing the nearest
// screen coordinate: empty space must not act on a molecule behind the menu.
function atomAtEvent(event) {
  const point = viewer.mouseXY(event.pageX, event.pageY);
  const candidates = Object.values(objects).flatMap(o => o.model.selectedAtoms({})
    .filter(a => componentVisible(o, a) && Object.keys(a.style || {}).length));
  const hits = viewer.targetedObjects(point.x, point.y, candidates);
  return hits.length ? hits[0].clickable : null;
}
function contextPayload(event) {
  const atom = atomAtEvent(event);
  const name = atom && atom._wdObject, o = name && objects[name];
  let count = 0;
  Object.keys(selection).forEach(n => { count += editableAtoms(n)
    .filter(a => selection[n].has(a._wdIdentity.atom_index)).length; });
  return {x: event.clientX, y: event.clientY, name: name || null,
    kind: o ? o.kind : null, target: atom ? describe(atom) : null,
    editable: !!(o && o.editable), atomIndex: atom && atom._wdIdentity ? atom._wdIdentity.atom_index : null,
    representation: o ? o.rep : null, color: o ? o.color : null,
    customColor: o ? o.customColor : null, components: o ? o.components : {},
    labels: !!(o && o.labels), selectionCount: count, focusCount: count + visualSelectionCount(), background: backgroundColor,
    objects: Object.keys(objects).map(n => ({name:n, visible:objects[n].visible})),
    selectionMode: selectionMode};
}
function selectContextTarget(name, atomIndex, mode) {
  const atom = editableAtoms(name).find(a => a._wdIdentity.atom_index === atomIndex);
  if (!atom) return;
  const previous = selectionMode; setSelectionMode(mode);
  pickSelection(name, atom, {}); setSelectionMode(previous);
}
container.addEventListener("contextmenu", e => {
  e.preventDefault(); e.stopImmediatePropagation();
  emit("context", contextPayload(e));
}, true);

container.addEventListener("dblclick", e => {
  if (e.button !== 0 || e.altKey || measureMode !== "none" || lasso) return;
  const atom = atomAtEvent(e); if (!atom) return;
  e.preventDefault(); e.stopImmediatePropagation();
  pickMolecularSelection(atom._wdObject, atom, e, true);
}, true);
