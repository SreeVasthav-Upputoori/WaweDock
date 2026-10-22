"""Pose comparison state shared by the desktop UI and standalone viewer audits."""
from __future__ import annotations


class PoseReviewSession:
    pose_name = "pose: current"

    def __init__(self, viewer):
        self.viewer = viewer
        self.enforced_visibility = {}
        self._protein_visibility = {}
        self._reference_visibility = {}
        self.reference = None
        self.has_pose = False

    def clear(self):
        self.enforced_visibility.clear()
        self._protein_visibility.clear()
        self._reference_visibility.clear()
        self.reference = None
        self.has_pose = False

    def apply(self, proteins, references, reference, *, locked=False, hide_protein=False):
        """Apply display-only masks; remember original visibility for later restoration."""
        visible = {**proteins, **references}
        desired = {}
        if not hide_protein:
            desired.update({n:v for n,v in self._protein_visibility.items() if n in proteins})
            self._protein_visibility.clear()
        if not locked or reference not in references:
            desired.update({n:v for n,v in self._reference_visibility.items() if n in references})
            self._reference_visibility.clear()
            self.reference = None
        enforced = {}
        if hide_protein:
            for name, state in proteins.items():
                self._protein_visibility.setdefault(name, state)
                enforced[name] = False
        if locked and reference in references:
            for name, state in references.items():
                self._reference_visibility.setdefault(name, state)
                enforced[name] = name == reference
            self.reference = reference
        desired.update(enforced)
        self.enforced_visibility = enforced
        for name, state in desired.items():
            if visible[name] != state:
                self.viewer.set_visible(name, state)

    def requested_visibility(self, name, visible):
        return self.enforced_visibility.get(name, visible)

    def show_pose(self, data, fmt="sdf", *, keep_camera=True):
        first = not self.has_pose
        self.viewer.add_object(self.pose_name, data, fmt, "ligand",
                               preserve_view=keep_camera, preserve_style=True)
        self.has_pose = True
        if first:
            self.viewer._call("setLigandColor", self.pose_name, "#ff9933", True)

    def clear_pose(self):
        if self.has_pose:
            self.viewer.remove_object(self.pose_name)
        self.has_pose = False

    def fit_view(self):
        names = ([self.reference] if self.reference else [])
        if self.has_pose:
            names.append(self.pose_name)
        if names:
            self.viewer._call("zoomToObjects", names)
