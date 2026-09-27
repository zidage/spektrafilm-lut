"""Grey-ramp probe: scene stops -> APD / ADX16 -> ACES round trip."""
import sys, numpy as np
from spektrafilm_lut_creator.aces_lmt import NegativeADXModel, NegativeADXSpec, aces_to_adx16

film = sys.argv[1] if len(sys.argv) > 1 else "kodak_vision3_250d"
balance = sys.argv[2] if len(sys.argv) > 2 else "grey"
metric = sys.argv[3] if len(sys.argv) > 3 else "print"
m = NegativeADXModel(NegativeADXSpec(film_profile=film, balance=balance, density_metric=metric))
print(film, "APD Dmin", m.pd_dmin.round(3), "ACES grey gain", m.aces_gain.round(3))
stops = np.arange(-10, 12.1, 1.0)
aces = (0.18 * 2.0 ** stops)[:, None] * np.ones(3)
cdd = m.cdd(aces)
adx = m.adx16(aces)
out = m.aces_out(aces)
ideal = aces_to_adx16(aces)
print(" stop |   CDD (R G B) raw      |   ADX16 R G B        | ideal ADX16 (IDT inverse) | out stops R G B")
for s, c, a, i, o in zip(stops, cdd, adx, ideal, out):
    print(f"{s:5.0f} | {c[0]:6.3f} {c[1]:6.3f} {c[2]:6.3f} | {a[0]:6.0f} {a[1]:6.0f} {a[2]:6.0f} | {i[0]:6.0f} {i[1]:6.0f} {i[2]:6.0f} | "
          + " ".join(f"{np.log2(max(v,1e-9)/0.18):6.2f}" for v in o))
