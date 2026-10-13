# O3 F0 full-wave carrier convergence

The preregistered paired mirror gate was `<= 3%`. It was not changed after
looking at the FEM results.

| Mesh pair | h bulk, mm | h defect, mm | Free mirror | Void mirror | Paired mirror | F0 |
|---|---:|---:|---:|---:|---:|---|
| independent | 1.5 | 0.35 | 0.629% | 0.887% | 25.783% | fail |
| independent | 1.0 | 0.25 | 0.192% | 0.312% | 10.184% | fail |
| common conforming | 1.5 | 0.35 | 0.877% | 0.887% | 4.695% | fail |
| common conforming | 1.0 | 0.25 | 0.318% | 0.311% | **1.805%** | **pass** |

The independent free/void matrices are individually accurate but their small
difference amplifies unrelated triangulation error. The accepted pair starts
from one conforming triangulation: the free model retains a filled circular
subdomain with aluminium properties, while the void model removes only those
interior elements. No operator, architecture, defect or gate changed during
this remediation.

The accepted response and frozen projections are in
[`fullwave_carrier_commonmesh_h1p0_v0/`](../fullwave_carrier_commonmesh_h1p0_v0/README.md).
The earlier failed directories are retained as the audit trail. F0 is a
carrier physics gate, not a broadband/transient detection result.
