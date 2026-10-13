Possible next steps, if you want any of them:

A split rule that doesn't flip under small perturbations. This would give the triangles the quads' skew insensitivity without the quads' assembly cost. It's also the open question behind keeping the zero tip.
Faster quad assembly. Most of the 9–19× cost is in the touching pairs (shared sides, shared corners and self pairs).
Splitting badly distorted quads into triangles. Quads whose corner Jacobians differ by a large factor (the "slivers" mesh) are integrated less accurately, which shows up as a high column defect. Splitting them would fix that.